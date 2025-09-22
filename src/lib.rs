use nom::AsBytes;
use noodles::bam;
use noodles::bam::io::Writer;
use noodles::sam::Header;
use noodles::sam::alignment::record::mapping_quality::MappingQuality;
use noodles::bgzf as bgzf;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::iter::zip;
use std::path::Path;
use std::path::PathBuf;
use std::str;
use std::sync::{Arc, Mutex};

mod bam_index;
mod hmm;
mod motif_correction;
mod repeats;
mod io;

use crate::bam_index::check_bai;
use crate::hmm::{Module, Hmm};
use crate::repeats::TandemRepeat;
use crate::io::get_modules;

pub fn run(
    bam_file: &Path, motif_file: &Path, output: String, out_bam: bool,
    params: (bool, u8, Option<char>, bool)
) {
    check_bai(bam_file);

    let header = header(bam_file);
    let out_bam = if out_bam { Some(Arc::new(Mutex::new(init_bam(&output, &header)))) } else { None };
    let out_tsv = Arc::new(Mutex::new(init_tsv(&output)));

    let motif_records = read_motifs(motif_file);
    motif_records.par_iter().for_each(|motif_record| {
        process_motif(motif_record, bam_file, params, out_tsv.clone(), out_bam.clone());
    });

    println!("Annotation finished successfully.");
    // TODO:
    // sort bam
    // create bai index
    //     let filename = args.output.to_string() + ".bam";
    //     check_bai(filename);
}

fn process_motif(
    motif_record: &(Vec<u8>, TandemRepeat, Vec<u8>), 
    bam_file: &Path,
    params: (bool, u8, Option<char>, bool),
    out_tsv: Arc<Mutex<File>>,
    out_bam: Option<Arc<Mutex<Writer<bgzf::Writer<File>>>>>,
) {
    // load bam
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_file)
        .expect("Unable to read the associated index (.bai).");
    let header = reader.read_header().unwrap();

    //  select relevant reads
    let (left_flank, repeat, right_flank) = motif_record;
    let (dedup, q, score, print_quality) = params;
    let region_str = format!("{}:{}-{}", repeat.reference, repeat.start + 1, repeat.end);
    let region = region_str.parse().unwrap();
    let reads: Vec<_> = reader
        .query(&header, &region).unwrap()
        .map(|x| x.expect("Incorrect read."))
        .collect();
    let raw_count = reads.len();
    let reads: Vec<_> = reads.into_iter()
        .filter(|x| !x.sequence().is_empty())
        .filter(|x| !(dedup && x.flags().is_duplicate()))
        .filter(|x| !mapq_less_than(x, q))
        .collect();
    let filt_count = reads.len();
    println!("{region_str}: {filt_count}/{raw_count}");

    //  build HMM
    let modules = get_modules(left_flank, repeat, right_flank);
    let model = Hmm::from(&modules).log();

    let (annotation, annotated_reads) = annotate_reads(reads.into_iter(), model, repeat, score, print_quality);

    // write to files
    out_tsv.lock().unwrap().write_all(annotation.as_bytes()).expect("Cannot write to output file.");
    match out_bam {
        None => {},
        Some(mutex) => {
            let mut writer = mutex.lock().unwrap();
            for record in annotated_reads {
                writer.write_record(&header, &record).expect("Cannot write to out bam.");
            }
        }
    }
}

fn annotate_reads<T>(reads: T, model: Hmm, repeat: &TandemRepeat, score: Option<char>, print_quality: bool)
    -> (String, Vec<bam::Record>)
where
    T: Iterator<Item = bam::Record>,
{
    let mut annotation_str = String::new();
    let mut annotated_reads = Vec::<bam::Record>::new();
    for (i, read) in reads.enumerate() {

        annotated_reads.push(read.clone());

        let seq: Vec<_> = read.sequence().iter().collect();
        let qual: Vec<_> = read.quality_scores().as_ref().iter().map(|x| x + 33).collect();

        let qual_mod =
            if let Some(x) = score { vec![x as u8; seq.len()] } 
            else { qual.clone() };

        let qual_str =
            if print_quality { qual.clone() }
            else { "No quality".bytes().collect() };

        let (likelihood, annotation) = model.log_predict(&seq, &qual_mod);

        let (new_annot, reconstructed_read) = model.realign(&annotation, &seq);
        let reconstructed_reference = model.reconstruct_sequence(&new_annot);
        let mods = model.reconstruct_mod_ids(&new_annot);

        // b"name\tmotif\tread_sn\tread_id\tmate_order\tread\treference\tmodules\tquality\tlog_likelihood\n"
        let name: String = match &repeat.name {
            Some(x) => x.to_string(),
            None    => "None".to_string(),
        };
        let motif = repeat;
        let read_sn = i;
        let read_id = str::from_utf8(read.name().unwrap().as_bytes()).unwrap();
        let mate_order = mate_order(&read);
        let read = str::from_utf8(&reconstructed_read).unwrap();
        let reference = str::from_utf8(&reconstructed_reference).unwrap();
        let modules = str::from_utf8(&mods).unwrap();
        let quality = str::from_utf8(&qual_str).unwrap();
        let log_likelihood = likelihood;

        let mlen = mods.len();
        let mut left_bg = 0;
        while left_bg < mlen && mods[left_bg] == b'-' { left_bg += 1; }
        let mut right_bg = 0;
        while right_bg < mlen && mods[(mlen - 1) - right_bg] == b'-' { right_bg += 1; }

        let mismatches_str = generate_mismatches(&reconstructed_read, &reconstructed_reference);
        let n_deletions = mismatches_str.bytes().filter(|x| *x == b'D').count();
        let n_insertions = mismatches_str.bytes().filter(|x| *x == b'I').count();
        let n_mismatches = mismatches_str.bytes().filter(|x| *x == b'M').count();

        let n_modules = repeat.copy_number.len() + 2;

        let mut module_bases: Vec<u8> = Vec::with_capacity(n_modules);
        let mut module_repetitions: Vec<u8> = Vec::with_capacity(n_modules);
        let mut module_sequences: Vec<String> = Vec::with_capacity(n_modules);
        for i in 0..n_modules {
            let mb = get_module_bases(&mods, i);
            let mr = get_module_repetitions(mb, &repeat.copy_unit, i);
            let ms = get_module_sequences(&mods, i);
            module_bases.push(mb);
            module_repetitions.push(mr);
            module_sequences.push(ms);
        }

        let module_bases = module_bases.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        let module_repetitions = module_repetitions.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        let module_sequences = module_sequences.join(",");

        let line = format!("\
            {name}\t{motif}\t\
            {read_sn}\t{read_id}\t{mate_order}\t{quality}\t{log_likelihood}\t\
            {read}\t\
            {reference}\t\
            {n_modules}\t{left_bg}\t{module_bases}\t{right_bg}\t{module_repetitions}\t{module_sequences}\t\
            {modules}\t\
            {n_deletions}\t{n_insertions}\t{n_mismatches}\t\
            {mismatches_str}\n\
            "
        );
        // let line = format!("\
        //     {name}\t{motif}\t{read_id}\n\
        //     {read}\n\
        //     {reference}\n\
        //     {modules}\n\
        //     {mismatches_str}\n\
        //     {read_sn}\t{mate_order}\t{quality}\t{log_likelihood}\n\
        //     {n_modules}\t{left_bg}\t{module_bases}\t{right_bg}\t{module_repetitions}\t{module_sequences}\n\
        //     {n_deletions}\t{n_insertions}\t{n_mismatches}\n\
        //     "
        // );
        annotation_str.push_str(&line);
    }
    return (annotation_str, annotated_reads);
}

fn get_module_sequences(_mods: &[u8], _idx: usize) -> String {
    let ms = "".to_string();  // TODO: finish this
    return ms
}

fn get_module_repetitions(mb: u8, copy_units: &[Vec<u8>], idx: usize) -> u8 {
    if mb == 0 { return 0; }
    if idx == 0 { return 1; }
    if idx == copy_units.len() + 1 { return 1; }
    if idx > copy_units.len() + 1 { panic!("This should never happen."); }
    let copy_len: u8 = copy_units[idx - 1].len().try_into().unwrap();
    return mb / copy_len;
}

fn get_module_bases(mods: &[u8], idx: usize) -> u8 {
    const ASCII_ZERO: usize = 48;
    let idx: u8 = (idx + ASCII_ZERO).try_into().unwrap();
    let count = mods.iter().filter(|&&x| x == idx).count();
    return count.try_into().unwrap();
}

fn header<P: AsRef<Path>>(bam_filename: P) -> Header {
    let file = File::open(bam_filename).unwrap();
    let header = bam::io::Reader::new(file).read_header().unwrap();
    return header;
}

fn init_tsv(filename: &str) -> File {
    let mut out = File::create(filename).expect("Cannot open file for writing.");
    let line = b"\
    name\tmotif\t\
    read_sn\tread_id\tmate_order\tquality\tlog_likelihood\t\
    read\t\
    reference\t\
    n_modules\tleft_bg\tmodule_bases\tright_bg\tmodule_repetitions\tmodule_sequences\t\
    modules\t\
    n_deletions\tn_insertions\tn_mismatches\t\
    mismatches_str\n\
    ";

    out.write_all(line).expect("Cannot write to output file.");
    return out;
}

fn init_bam(tsv_file: &str, header: &Header) -> Writer<bgzf::Writer<File>> {
    let mut filename = PathBuf::from(tsv_file);
    filename.set_extension("bam");
    let new_bam = File::create(filename).expect("Cannot open file for writing.");
    let mut writer = bam::io::Writer::new(new_bam);
    writer.write_header(header).unwrap();
    return writer;
}

fn mapq_less_than(rec: &bam::Record, x: u8) -> bool {
    let x = MappingQuality::new(x)
        .expect("Mapq is from 0 to 254. 255 is reserved for None.");
    let Some(q) = rec.mapping_quality() else { return false };
    return q < x;
}

fn read_motifs(filename: &Path) -> Vec<(Vec<u8>, TandemRepeat, Vec<u8>)> {
    let file = File::open(filename).expect("Cannot find nomenclature file.");
    let reader = BufReader::new(file);

    // let crash = |_| panic!("line {}: Nomenclature {} malformatted. Accepted format is <chr>:g.<start>_<end><sequence>[repetitions].", i+1, split[1])
    // assert!(split.len() == 4,
    // "Malformatted line, expected format is <name>\\t<left_flank>\\t<hgvs_nomenclature>\\t<right_flank>\\n.");
    // Accepted format is <chr>:g.<start>_<end><sequence>[repetitions].

    let mut result = Vec::new();
    for line in reader.lines().skip(1) {
        let line = line.expect("Cannot read line from nomenclature file.").trim().to_owned();
        let split: Vec<_> = line.split('\t').collect();

        let name = split[0].to_owned();
        let left_flank = split[2].as_bytes().to_owned();
        let mut repeat: TandemRepeat = split[1].parse().expect("Malformatted nomenclature found.");
        repeat.name = Some(name);
        let right_flank = split[3].as_bytes().to_owned();

        result.push((left_flank, repeat, right_flank));
    }
    return result;
}

fn generate_mismatches(read: &[u8], reference: &[u8]) -> String {
    let mut result = String::with_capacity(read.len());
    for (x, y) in zip(read, reference) {
        match (x, y) {
            (_,    b'-') => { result.push('_'); }
            (b'_', _   ) => { result.push('D'); }
            (_,    b'_') => { result.push('I'); }
            (x,    y   ) => {
                if x == y { result.push('_'); } else { result.push('M'); }
            }
        }
    }
    return result;
}

fn mate_order(read: &bam::Record) -> String {
    if read.flags().is_first_segment() { "1".to_string() }
    else if read.flags().is_last_segment() { "2".to_string() } 
    else {
        // println!("Read {} does not have pair information.", read.read_name().unwrap());
        "0".to_string()
    }
}

#[test]
fn can_load_bam() {
    use noodles::bam::io::reader;
    let mut reader = reader::Builder
        .build_from_path("./data/test/mini.bam").unwrap();

    reader.read_header().unwrap(); // this is necessary here
    for result in reader.records() {
        let record = result.unwrap();
        println!("{:?}", record);
    }
}

#[test]
fn count_present() {
    let hgvs = File::open("data/test/HGVS.txt").unwrap();
    let reader = BufReader::new(hgvs);

    // let mut present_count = 0;
    let present_count = 0;
    let mut max_count = 0;
    for line in reader.lines() {
        let line = line.unwrap().trim().to_owned();
        let _tr: TandemRepeat = line.parse().unwrap();
        // if is_present(&tr, &references) {
        //     present_count += 1;
        // } else {
            // println!("{}", tr);
            // print_diff(&tr, &references);
            // println!();
        // }
        max_count += 1;
    }
    println!("Present repeats: {}/{}", present_count, max_count);
}

#[test]
fn can_read_and_parse_hgvs_file() {
    let file = "./data/test/mini_HGVS.txt";
    let file = File::open(file).unwrap();
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line.unwrap();
        let line = line.trim();
        println!("{}", line);
        let tr: TandemRepeat = line.parse().unwrap();
        println!("{:?}", tr);
    }
}

#[test]
fn tmp_fn_name() {
    // CTTCCTCCTCCTCATCGGTGGCGGCGGCGGCGGCGTCAGGCCAGTGCCGCGGCTTTCTCTCCGCGCCTGTGTTCGCCGGGACGCATTCGGGGCGGGCGGCGGCGGCGGCAGCGGCGGCCGCGGCGGCGGGCGGGGCCGCCCCCCGCCT
    // CCCFFFFFHHHHGIIJJIHIJIJJJHGHDDDBDDBBDDDDDDDDCCCC5@BDDBBCCCCACDBBDBDBBCCCCCBDDBBD@BDDBDDDDDDDDDD@5.5&)0)&5))0)&)2&55)0&&&&&)&&)&))&&&&&)&&&&)&&50)&&&
    // -----------------------------------------------------------------00000000000000000000000000000011111111111111111111111111111111111111111111122222222
    // GCG[4]GCA[1]GCG[2]GCC[1]GCG[3]G[1]GCG[1]GGGCCGCC[1]
    //
    // I think annotation is wrong
    // Independent of annotation, seq nomenclature is wrong as well
    //
    // SPD     chr2:g.176093059_176093103GCG[15]       CCTGTGTTCGCCGGGACGCATTCGGGGCGG  TCCGGCTTTGCGTACCCCGGGACCTCTGAG
    // result 15

    // HISEQ1:26:HA2RRADXX:1:1203:16720:7919
    let seq:  Vec<u8> = b"CTTCCTCCTCCTCATCGGTGGCGGCGGCGGCGGCGTCAGGCCAGTGCCGCGGCTTTCTCTCCGCGCCTGTGTTCGCCGGGACGCATTCGGGGCGGGCGGCGGCGGCGGCAGCGGCGGCCGCGGCGGCGGGCGGGGCCGCCCCCCGCCT".to_vec();
    let qual: Vec<u8> = b"CCCFFFFFHHHHGIIJJIHIJIJJJHGHDDDBDDBBDDDDDDDDCCCC5@BDDBBCCCCACDBBDBDBBCCCCCBDDBBD@BDDBDDDDDDDDDD@5.5&)0)&5))0)&)2&55)0&&&&&)&&)&))&&&&&)&&&&)&&50)&&&".to_vec();

    // SPD
    let left_flank:  Vec<u8> = b"CCTGTGTTCGCCGGGACGCATTCGGGGCGG".to_vec();
    let right_flank: Vec<u8> = b"TCCGGCTTTGCGTACCCCGGGACCTCTGAG".to_vec();
    let repeat: TandemRepeat = "chr2:g.176093059_176093103GCG[15]".parse().expect("Malformatted nomenclature found.");

    let modules = get_modules(&left_flank, &repeat, &right_flank);
    let model = Hmm::from(&modules).log();
    let (_likelihood, annotation) = model.log_predict(&seq, &qual);

    for x in &annotation { print!("{}", x / 10); }
    println!();
    for x in &annotation { print!("{}", x % 10); }
    println!();
    println!("{}", str::from_utf8(&seq).unwrap());
    println!("{}", str::from_utf8(&qual).unwrap());

    let (partition, mod_ids) = model.partition_to_units(&annotation);

    let exp_split: Vec<&[u8]> = vec![
        b"CTTCCTCCTCCTCATCGGTGGCGGCGGCGGCGGCGTCAGGCCAGTGCCGCGGCTTTCTCTCCGCG",
        b"CCTGTGTTCGCCGGGACGCATTCGGGGCGG",
        b"GCG", b"GCG", b"GCG", b"GCG", b"GCA", b"GCG", b"GCG", b"GCC", b"GCG", b"GCG", b"GCG", b"GGC", b"GGG", b"GCC", b"GCC",
        b"CCCCGCCT"
    ];

    let exp_mod_ids: Vec<_> = vec![
        usize::MAX, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2
    ];

    for (i, p) in partition.into_iter().enumerate() {
        println!("{}", str::from_utf8(&seq[p]).unwrap());
        println!("{}", str::from_utf8(exp_split[i]).unwrap());
        println!("{} {}", mod_ids[i], exp_mod_ids[i])
    }
}

