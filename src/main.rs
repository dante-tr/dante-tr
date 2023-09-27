use noodles::bam as bam;
use noodles::fasta as fasta;
use noodles::sam::record::quality_scores::Score;
use noodles::sam::alignment::Record;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::str;
use std::sync::{Arc, Mutex};
use clap::Parser;

mod cli;
mod consistency;
mod hmm;
mod repeats;
mod motif_correction;

use crate::cli::Args;
use crate::consistency::ensure_consistency;
use crate::motif_correction::correct_repeats;
use crate::hmm::{HMM, Module};
use crate::repeats::TandemRepeat;

fn main() {
    let args = Args::parse();
    let flank = 20;

    let bam_refs = read_bam_refs(&args.bam_file);
    let references = read_reference(&args.ref_file);
    let repeats = read_nomenclature(&args.hgvs_file);

    let (references, repeats) = ensure_consistency(bam_refs, references, repeats);
    let valid_repeats = correct_repeats(&references, &repeats);

    let mut out = File::create(&args.out_file).expect("Cannot open file for writing.");
    out.write_all(b"motif\tread_sn\tread_id\tmate_order\tread\treference\tmodules\tlog_likelihood\n")
        .expect("Cannot write to output file.");

    let out = Arc::new(Mutex::new(out));
    valid_repeats.par_iter().for_each(|repeat| {
        // load bam
        let mut reader = bam::indexed_reader::Builder::default()
            .build_from_path(&args.bam_file)
            .expect("Unable to read the associated index (.bai).");
        let header = reader.read_header().unwrap();

        //  build HMM
        let modules = get_modules(&repeat, &references, flank);
        let model = HMM::from(&modules).log();

        //  select relevant reads
        let tmp = format!("{}:{}-{}", repeat.reference, repeat.start+1, repeat.end);
        let region = tmp.parse().unwrap();
        let reads = reader.query(&header, &region).unwrap();

        let mut buffer = String::new();
        for (i, read) in reads.enumerate() {
            let read: Record = read.expect("Incorrect read.");
            let seq: Vec<_> = read.sequence().as_ref().iter().map(|&x| x.into()).collect();
            let qual: Vec<_> = read.quality_scores().as_ref().iter().map(|&x| remap(x)).collect();
            let (likelihood, annotation) = model.log_predict(&seq, &qual);

            let reconstructed_reference = model.reconstruct_sequence(&annotation);
            let reconstructed_read = model.realign_read(&annotation, &seq); 
            let mods = model.reconstruct_mod_ids(&annotation);

            buffer.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                repeat, i,
                read.read_name().unwrap(),
                mate_order(&read),
                str::from_utf8(&reconstructed_read).unwrap(),
                str::from_utf8(&reconstructed_reference).unwrap(),
                str::from_utf8(&mods).unwrap(),
                likelihood
            ));
        }
        out.lock().unwrap().write_all(buffer.as_bytes())
            .expect("Cannot write to output file.");
    })
}

fn read_bam_refs(filename: &str) -> HashMap<String, usize> {
    let mut result = HashMap::new();

    let file = File::open(filename).unwrap();
    let header = bam::Reader::new(file).read_header().unwrap();

    for s in header.reference_sequences().iter() {
        let name = s.0.to_string();
        let length = s.1.length().get();
        result.insert(name.clone(), length);
    }
    return result;
}

fn read_reference(filename: &str) -> HashMap<String, Vec<u8>> {
    let mut reader = fasta::reader::Builder
        .build_from_path(filename).unwrap();

    let mut result = HashMap::new();
    for record in reader.records() {
        let record = record.unwrap();

        result.insert(
            record.name().to_string(),
            (&record.sequence()[..]).to_vec()
            // ^- Is there a better way to get Vec<u8>
            // Do I need Vec<u8>? Cannot I leave it as Sequence?
        );
    }
    return result;
}

fn read_nomenclature(filename: &str) -> Vec<TandemRepeat> {
    let mut repeats = Vec::new();

    let file = File::open(filename).expect("Cannot find nomenclature file.");
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line
            .expect("Cannot read line from nomenclature file.")
            .trim().to_owned();
        let repeat = line.parse().expect("Cannot parse nomenclature.");
        repeats.push(repeat);
    }
    return repeats;
}

fn get_modules(
    repeat: &TandemRepeat, refs: &HashMap<String, Vec<u8>>, flank_size: usize
) -> Vec<Module> {
    let refseq = refs.get(&repeat.reference).unwrap(); // safe due to nomenclature check
    assert!(repeat.start >= flank_size,
        "Cannot create left flank of size {flank_size} for repeat {repeat}.");
    let left_flank = &refseq[(repeat.start-flank_size)..repeat.start];
    assert!(repeat.end + flank_size <= refseq.len(),
        "Cannot create right flank of size {flank_size} for repeat {repeat}.");
    let right_flank = &refseq[repeat.end..(repeat.end+flank_size)];

    let mut modules = Vec::new();
    modules.push(left_flank.into());
    modules_add_motif(&mut modules, &repeat);
    modules.push(right_flank.into());
    return modules;
}

fn modules_add_motif(modules: &mut Vec<Module>, motif: &TandemRepeat) {
    for i in 0..motif.copy_unit.len() {
        modules.push((&motif.copy_unit[i][..], motif.copy_number[i]).into())
    }
}

fn mate_order(read: &Record) -> String {
    if read.flags().is_first_segment() { "1".to_string() }
    else if read.flags().is_last_segment() { "2".to_string() } 
    else {
        println!("Read {} does not have pair information.", read.read_name().unwrap());
        "0".to_string()
    }
}

fn remap(x: Score) -> u8 {
    let c: char = x.into();
    return c as u8;
}

#[cfg(test)]
mod tests {
    use hgvs::parser::HgvsVariant;
    use std::fs::File;
    use std::io::{prelude::*, BufReader};
    use super::*;

    #[test]
    fn can_load_bam() {
        let mut reader = bam::indexed_reader::Builder::default()
            .build_from_path("data/mini.bam").unwrap();

        let header = reader.read_header().unwrap();

        for result in reader.records(&header) {
            let record = result.unwrap();
            println!("{:?}", record);
        }
    }

    #[test]
    fn can_load_fasta() {
        let sequences = read_reference("data/chromosomeX.fna");
        let hgvs = File::open("data/mini_HGVS.txt").unwrap();
        let reader = BufReader::new(hgvs);

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
    fn count_present() {
        let references = read_reference("data/chromosomeX.fna");
        let hgvs = File::open("data/HGVS.txt").unwrap();
        let reader = BufReader::new(hgvs);

        let mut present_count = 0;
        let mut max_count = 0;
        for line in reader.lines() {
            let line = line.unwrap().trim().to_owned();
            let tr: TandemRepeat = line.parse().unwrap();
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
        let file = "data/mini_HGVS.txt";
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
}

