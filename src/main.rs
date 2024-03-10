use clap::Parser;
use noodles::bam;
use noodles::bam::Writer;
use noodles::sam::Header;
use noodles::sam::alignment::Record;
use noodles::sam::record::quality_scores::Score;
use noodles::sam::record::MappingQuality;
use noodles::bgzf as bgzf;
use rayon::prelude::*;
use core::panic;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;
use std::str;
use std::sync::{Arc, Mutex};

mod bam_index;
mod cli;
mod consistency;
mod hmm;
mod motif_correction;
mod repeats;
mod io;

use crate::bam_index::check_bai;
use crate::cli::Args;
use crate::consistency::ensure_consistency;
use crate::hmm::{Module, Hmm};
use crate::motif_correction::correct_repeats;
use crate::repeats::TandemRepeat;
use crate::io::{get_modules, read_reference};

fn main() {
    let args = Args::parse();

    // checks
    check_bai(&args.bam_file);
    let header = header(&args.bam_file);

    let bam_refs = read_bam_refs(&header);
    let references = read_reference(&args.ref_file);
    let repeats = read_motifs(&args.motif_file);

    let (references, repeats) = ensure_consistency(bam_refs, references, repeats);

    let mut repeats = repeats;
    if args.correction { repeats = correct_repeats(&references, &repeats); }
    let repeats = repeats;

    { // scope for out_tsv and out_bam
    let out_tsv = init_tsv(&args.output);
    let out_tsv = Arc::new(Mutex::new(out_tsv));

    let out_bam = if args.out_bam {
        let out_bam = init_bam(&args.output, &header);
        let out_bam = Arc::new(Mutex::new(out_bam));
        Some(out_bam)
    } else {
        None
    };

    repeats.par_iter().for_each(|repeat| {
        // load bam
        let mut reader = bam::indexed_reader::Builder::default()
            .build_from_path(&args.bam_file)
            .expect("Unable to read the associated index (.bai).");
        let header = reader.read_header().unwrap();

        //  build HMM
        let modules = get_modules(repeat, &references, args.flank);
        let model = Hmm::from(&modules).log();

        //  select relevant reads
        let tmp = format!("{}:{}-{}", repeat.reference, repeat.start + 1, repeat.end);
        let region = tmp.parse().unwrap();
        let reads = reader
            .query(&header, &region).unwrap()
            .map(|x| x.expect("Incorrect read."))
            .filter(|x| !(args.dedup && x.flags().is_duplicate()))
            .filter(|x| !mapq_less_than(x, args.q));

        let (annotation, annotated_reads) = annotate_reads(reads, model, repeat, &args);

        // write to files
        out_tsv.lock().unwrap().write_all(annotation.as_bytes()).expect("Cannot write to output file.");
        match &out_bam {
            None => {},
            Some(mutex) => {
                let mut writer = mutex.lock().unwrap();
                for record in annotated_reads {
                    writer.write_record(&header, &record).expect("Cannot write to out bam.");
                }
            }
        }
    });
    } // end scope for out_tsv and out_bam

    println!("Annotation finished successfully.");
    // TODO:
    // sort bam
    // create bai index
    //     let filename = args.output.to_string() + ".bam";
    //     check_bai(filename);
}

fn header<P: AsRef<Path>>(bam_filename: P) -> Header {
    let file = File::open(bam_filename).unwrap();
    let header = bam::Reader::new(file).read_header().unwrap();
    return header;
}

fn init_tsv(filename: &str) -> File {
    let mut out = File::create(filename).expect("Cannot open file for writing.");
    out.write_all(
        b"name\tmotif\tread_sn\tread_id\tmate_order\tread\treference\tmodules\tlog_likelihood\n"
    ).expect("Cannot write to output file.");
    return out;
}

fn init_bam(tsv_file: &str, header: &Header) -> Writer<bgzf::Writer<File>> {
    let mut filename = PathBuf::from(tsv_file);
    filename.set_extension("bam");
    let new_bam = File::create(filename).expect("Cannot open file for writing.");
    let mut writer = bam::Writer::new(new_bam);
    writer.write_header(header).unwrap();
    return writer;
}

fn mapq_less_than(rec: &Record, x: u8) -> bool {
    let x = MappingQuality::new(x)
        .expect("Mapq is from 0 to 254. 255 is reserved for None.");
    let Some(q) = rec.mapping_quality() else { return false };
    return q < x;
}

fn annotate_reads<T>(reads: T, model: Hmm, repeat: &TandemRepeat, args: &Args)
    -> (String, Vec<Record>)
where
    T: Iterator<Item = Record>,
{
    let mut annotation_str = String::new();
    let mut annotated_reads = Vec::<Record>::new();
    for (i, read) in reads.enumerate() {

        annotated_reads.push(read.clone());

        let seq: Vec<_> = read.sequence().as_ref().iter().map(|&x| x.into()).collect();
        let qual = if let Some(x) = args.score {
            vec![x as u8; seq.len()]
        } else {
            read.quality_scores().as_ref().iter().map(|&x| remap(x)).collect()
        };


        let (likelihood, annotation) = model.log_predict(&seq, &qual);

        let (new_annot, reconstructed_read) = model.realign(&annotation, &seq);
        let reconstructed_reference = model.reconstruct_sequence(&new_annot);
        let mods = model.reconstruct_mod_ids(&new_annot);
        let name: String = match &repeat.name {
            Some(x) => x.to_string(),
            None    => "None".to_string(),
        };

        annotation_str.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            name, repeat, i,
            read.read_name().unwrap(),
            mate_order(&read),
            str::from_utf8(&reconstructed_read).unwrap(),
            str::from_utf8(&reconstructed_reference).unwrap(),
            str::from_utf8(&mods).unwrap(),
            likelihood
        ));
    }
    return (annotation_str, annotated_reads);
}

fn read_bam_refs(header: &Header) -> HashMap<String, usize> {
    let mut result = HashMap::new();
    for s in header.reference_sequences().iter() {
        let name = s.0.to_string();
        let length = s.1.length().get();
        result.insert(name.clone(), length);
    }
    return result;
}

fn read_motifs(filename: &str) -> Vec<TandemRepeat> {
    let repeats = if is_named_format(filename) {
        read_nomenclature_with_names(filename)
    } else {
        read_nomenclature(filename)
    };

    return repeats;
}

fn is_named_format(filename: &str) -> bool {
    let file = File::open(filename).expect("Cannot find nomenclature file.");
    let reader = BufReader::new(file);
    let count = reader
        .lines().next().expect("Empty nomenclature file?")
        .expect("Cannot read line from nomenclature file.")
        .trim().to_owned()
        .split('\t').count();

    match count {
        1 => { false },
        2 => { true },
        _ => { panic!("Unexpected number of columns in nomenclature file.") }
    }
}

fn read_nomenclature_with_names(filename: &str) -> Vec<TandemRepeat> {
    let file = File::open(filename).expect("Cannot find nomenclature file.");
    let reader = BufReader::new(file);

    let mut repeats = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.expect("Cannot read line from nomenclature file.").trim().to_owned();
        let split: Vec<_> = line.split('\t').collect();
        assert!(split.len() == 2, "Malformatted line, expected format is <name>\\t<hgvs_nomenclature>\\n.");
        let name = split[0].to_owned();
        let mut repeat: TandemRepeat = split[1].parse()
            .unwrap_or_else(|_| panic!("\
                line {}: Nomenclature {} malformatted. \
                Accepted format is <chr>:g.<start>_<end><sequence>[repetitions].\
            ", i+1, split[1]));
        repeat.name = Some(name);

        repeats.push(repeat);
    }
    return repeats;
}

fn read_nomenclature(filename: &str) -> Vec<TandemRepeat> {
    let mut repeats = Vec::new();

    let file = File::open(filename).expect("Cannot find nomenclature file.");
    let reader = BufReader::new(file);

    for (i, line) in reader.lines().enumerate() {
        let line = line.expect("Cannot read line from nomenclature file.").trim().to_owned();
        let repeat = line.parse()
            .unwrap_or_else(|_| panic!("\
                line {}: Nomenclature {} malformatted.\
                Accepted format is <chr>:g.<start>_<end><sequence>[repetitions].\
            ", i+1, line));

        repeats.push(repeat);
    }
    return repeats;
}

fn mate_order(read: &Record) -> String {
    if read.flags().is_first_segment() { "1".to_string() }
    else if read.flags().is_last_segment() { "2".to_string() } 
    else {
        // println!("Read {} does not have pair information.", read.read_name().unwrap());
        "0".to_string()
    }
}

fn remap(x: Score) -> u8 {
    let c: char = x.into();
    return c as u8;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::iter::zip;
    use hgvs::parser::HgvsVariant;
    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn can_load_bam() {
        let mut reader = bam::indexed_reader::Builder::default()
            .build_from_path("data/test/mini.bam").unwrap();

        let header = reader.read_header().unwrap();

        for result in reader.records(&header) {
            let record = result.unwrap();
            println!("{:?}", record);
        }
    }

    #[test]
    fn can_load_fasta() {
        // let sequences = read_reference("data/chromosomeX.fna");
        // let hgvs = File::open("data/mini_HGVS.txt").unwrap();
        // let reader = BufReader::new(hgvs);

        // let expected = vec![
        //     false, false, true, false, false, false, false, true, true, false
        // ];
        // for (i, line) in reader.lines().enumerate() {
            // let line = line.unwrap();
            // let line = line.trim();
            // let tr: TandemRepeat = line.parse().unwrap();
            // let is_correct = is_present(&tr, &sequences);
            // assert_eq!(is_correct, expected[i]);
        // }
    }

    #[test]
    fn can_read_tsv_nomenclature() {
        let filename = "data/test/nomenclature_hgs_1Q_with_names.tsv";
        let motifs1 = read_motifs(filename);

        let filename = "data/test/nomenclature_hgs_1Q_wo_names.tsv";
        let motifs2 = read_motifs(filename);

        for (m1, m2) in zip(motifs1, motifs2) {
            assert_eq!(m1.reference, m2.reference);
            assert_eq!(m1.start, m2.start);
            assert_eq!(m1.end, m2.end);
            assert_eq!(m1.copy_unit, m2.copy_unit);
            assert_eq!(m1.copy_number, m2.copy_number);
        }
    }

    #[test]
    fn count_present() {
        // let references = read_reference("data/chromosomeX.fna");
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

//     fn print_diff(tr: &TandemRepeat, refs: &HashMap<String, Vec<u8>>) {
//         let n = 10;
//         let rflank = ref_region(refs, &tr.reference, tr.start-n, tr.start).unwrap();
//         let ref_repeat = ref_region(refs, &tr.reference, tr.start, tr.end).unwrap();
//         let lflank = ref_region(refs, &tr.reference, tr.end, tr.end+n).unwrap();
//         println!("{} {} {}", 
//             str::from_utf8(rflank).unwrap(),
//             str::from_utf8(ref_repeat).unwrap(),
//             str::from_utf8(lflank).unwrap()
//         );
//         println!("{} {} {}",
//             " ".repeat(n),
//             str::from_utf8(&tr.sequence()).unwrap(),
//             " ".repeat(n)
//         );
//     }

    #[test]
    fn can_parse_hgvs() {
        let _record = "NM_01234.5:c.456-6_*22A>T";
        let _record = "NC_000017.11:g.43091687del";
        let tmp: HgvsVariant = _record.parse().unwrap();
        println!("{:?}", tmp);

        println!("{}", tmp.accession().value);
    }

    #[test]
    fn can_read_and_parse_hgvs_file() {
        let file = "data/test/mini_HGVS.txt";
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
    fn does_not_overflow() {
        let references = read_reference("data/test/chromosomeX.fna");
        let motif: TandemRepeat = "NC_000023.11:g.284585_284614AC[15]".parse().unwrap();
        let repeats = vec![motif];
        correct_repeats(&references, &repeats);
    }
}
