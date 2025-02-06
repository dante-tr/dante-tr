use nom::AsBytes;
use noodles::bam;
use noodles::bam::io::Writer;
use noodles::sam::Header;
use noodles::sam::alignment::record::mapping_quality::MappingQuality;
use noodles::bgzf as bgzf;
use rayon::prelude::*;
use core::panic;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
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
    let tmp = format!("{}:{}-{}", repeat.reference, repeat.start + 1, repeat.end);
    let region = tmp.parse().unwrap();
    let reads = reader
        .query(&header, &region).unwrap()
        .map(|x| x.expect("Incorrect read."))
        .filter(|x| !x.sequence().is_empty())
        .filter(|x| !(dedup && x.flags().is_duplicate()))
        .filter(|x| !mapq_less_than(x, q));

    //  build HMM
    let modules = get_modules(left_flank, repeat, right_flank);
    let model = Hmm::from(&modules).log();

    let (annotation, annotated_reads) = annotate_reads(reads, model, repeat, score, print_quality);

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

fn header<P: AsRef<Path>>(bam_filename: P) -> Header {
    let file = File::open(bam_filename).unwrap();
    let header = bam::io::Reader::new(file).read_header().unwrap();
    return header;
}

fn init_tsv(filename: &str) -> File {
    let mut out = File::create(filename).expect("Cannot open file for writing.");
    out.write_all(
        b"name\tmotif\tread_sn\tread_id\tmate_order\tread\treference\tmodules\tquality\tlog_likelihood\n"
    ).expect("Cannot write to output file.");
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
            else { vec![] };


        let (likelihood, annotation) = model.log_predict(&seq, &qual_mod);

        let (new_annot, reconstructed_read) = model.realign(&annotation, &seq);
        let reconstructed_reference = model.reconstruct_sequence(&new_annot);
        let mods = model.reconstruct_mod_ids(&new_annot);
        let name: String = match &repeat.name {
            Some(x) => x.to_string(),
            None    => "None".to_string(),
        };

        annotation_str.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            name, repeat, i,
            str::from_utf8(read.name().unwrap().as_bytes()).unwrap(),
            mate_order(&read),
            str::from_utf8(&reconstructed_read).unwrap(),
            str::from_utf8(&reconstructed_reference).unwrap(),
            str::from_utf8(&mods).unwrap(),
            str::from_utf8(&qual_str).unwrap(),
            likelihood
        ));
    }
    return (annotation_str, annotated_reads);
}

fn read_motifs(filename: &Path) -> Vec<(Vec<u8>, TandemRepeat, Vec<u8>)> {
    let file = File::open(filename).expect("Cannot find nomenclature file.");
    let reader = BufReader::new(file);

    let mut result = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.expect("Cannot read line from nomenclature file.").trim().to_owned();
        let split: Vec<_> = line.split('\t').collect();
        assert!(split.len() == 4,
            "Malformatted line, expected format is <name>\\t<left_flank>\\t<hgvs_nomenclature>\\t<right_flank>\\n.");
        let name = split[0].to_owned();
        let left_flank = split[1].as_bytes().to_owned();
        let mut repeat: TandemRepeat = split[2].parse()
            .unwrap_or_else(|_| panic!("\
                line {}: Nomenclature {} malformatted. \
                Accepted format is <chr>:g.<start>_<end><sequence>[repetitions].\
            ", i+1, split[1]));
        let right_flank = split[3].as_bytes().to_owned();
        repeat.name = Some(name);

        result.push((left_flank, repeat, right_flank));
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

