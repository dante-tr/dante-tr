use clap::Parser;
use noodles::bam as bam;
use noodles::fasta as fasta;
use noodles::sam::record::quality_scores::Score;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::str;
use std::sync::{Arc, Mutex};

use crate::consistency::ensure_consistency;
use crate::hmm::HMM;
use crate::hmm::Module;
use crate::repeats::TandemRepeat as TandemRepeat;

mod consistency;
mod hmm;
mod repeats;

const FLANK: usize = 20;

// Predict short tandem repeat annotation
#[derive(Parser)]
pub struct Args {
    /// Reference in FASTA format
    #[arg(short='f')]
    pub ref_file: String,

    /// Reads mapped to reference in BAM format, index (.bai) has to be present
    #[arg(short='b')]
    pub bam_file: String,

    /// Repeats in HGVS nomenclature, one per line
    #[arg(short='n')]
    pub hgvs_file: String,

    /// Output file in TSV format.
    #[arg(short='o')]
    pub out_file: String,
}

fn main() {
    let args = Args::parse();

    let bam_refs = read_bam_refs(&args.bam_file);
    let references = read_reference(&args.ref_file);
    let repeats = read_nomenclature(&args.hgvs_file);

    let (references, repeats) = ensure_consistency(bam_refs, references, repeats);
    let valid_repeats = correct_repeats(&references, &repeats);

    let mut out = File::create(&args.out_file).expect("Cannot open file for writing.");
    out.write_all(b"motif\tread_sn\tread_id\tread\treference\tmodules\tlog_likelihood\n")
        .expect("Cannot write to output file.");

    let out = Arc::new(Mutex::new(out));
    valid_repeats.par_iter().for_each(|repeat| {
        // load bam
        let mut reader = bam::indexed_reader::Builder::default()
            .build_from_path(&args.bam_file)
            .expect("Unable to read the associated index (.bai).");
        let header = reader.read_header().unwrap();

        //  build HMM
        let modules = get_modules(&repeat, &references, FLANK);
        let model = HMM::from(&modules).log();

        //  select relevant reads
        let tmp = format!("{}:{}-{}", repeat.reference, repeat.start+1, repeat.end);
        let region = tmp.parse().unwrap();
        let reads = reader.query(&header, &region).unwrap();

        let mut buffer = String::new();
        for (i, read) in reads.enumerate() {
            let read = read.expect("Incorrect read.");
            let seq: Vec<_> = read.sequence().as_ref().iter().map(|&x| x.into()).collect();
            let qual: Vec<_> = read.quality_scores().as_ref().iter().map(|&x| remap(x)).collect();
            let (likelihood, annotation) = model.log_predict(&seq, &qual);

            let reconstructed_reference = model.reconstruct_sequence(&annotation);
            let reconstructed_read = model.realign_read(&annotation, &seq); 
            let mods = model.reconstruct_mod_ids(&annotation);

            buffer.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                repeat, i,
                read.read_name().unwrap(),
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

fn correct_repeats(refs: &HashMap<String, Vec<u8>>, repeats: &Vec<TandemRepeat>) -> Vec<TandemRepeat> {
    let mut valid_repeats = Vec::new();
    for motif in repeats.iter() {
        if is_present(&motif, &refs) {
            valid_repeats.push(motif.clone());
        } else {
            // eprintln!("Motif {} is not present. Correcting...", motif);
            let from = motif.start - FLANK;
            let to = motif.end + FLANK;
            let seq = ref_region(refs, &motif.reference, from, to)
                .expect("Unable to get reference region.");

            // let corrected_motif = correct_motif(&seq, &motif, FLANK);
            println!(
                "{}\n{}\n{}\n",
                motif,
                str::from_utf8(&seq).unwrap(),
                str::from_utf8(&motif.view(from, to)).unwrap(),
                // str::from_utf8(&corrected_motif.view(from, to)).unwrap(),
            );

            // valid_repeats.push(corrected_motif);
        }
    }
    return valid_repeats;
}

fn correct_motif(seq: &[u8], repeat: &TandemRepeat, flank: usize) -> TandemRepeat {
    let qual = b"~".repeat(seq.len());

    let mut modules = Vec::new();
    modules_add_motif(&mut modules, &repeat); // TODO: make this function?

    let model = HMM::from(&modules).log();
    let (_, annotation) = model.log_predict(&seq, &qual);

    let suggested_repeat = {
        let mut new_repeat = repeat.clone();
        let start = match annotation.iter().position(|&x| x != 0) {
            None => {
                eprintln!("Unable to match with reference.");
                return new_repeat;
            },
            Some(x) => { x }
        };
        let start = repeat.start - flank + start;
        new_repeat.start = start;
        new_repeat.end = start + repeat.sequence().len();
        new_repeat
    };

    let mut orig_motif = b"-".repeat(flank);
    orig_motif.extend_from_slice(&repeat.sequence());
    orig_motif.extend_from_slice(&b"-".repeat(flank));

    return suggested_repeat;
}

fn fn3(model: &HMM, annotation: &[usize]) -> (usize, usize) {
    let start: usize = 0;
    let end: usize = 6; //TODO: model.get_end();

    let mut m_start = usize::MIN;
    let mut m_end = usize::MAX;
    for (i, &state) in annotation.iter().enumerate() {
        if state == start { m_start = i; }
        if state == end && m_end == usize::MAX { m_end = i; }
    }
    return (m_start + 1, m_end);
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

fn ref_region<'a>(
    refseq: &'a HashMap<String, Vec<u8>>, id: &str, start: usize, end: usize
) -> Option<&'a[u8]> {
    let seq = match refseq.get(id) {
        None => { return None; },
        Some(x) => { x },
    };
    return Some(&seq[start..end]);
}

fn is_present(tr: &TandemRepeat, seq: &HashMap<String, Vec<u8>>) -> bool {
    let ref_repeat = match ref_region(seq, &tr.reference, tr.start, tr.end) {
        None => { return false; },
        Some(x) => { x },
    };
    let hgvs_repeat = &tr.sequence();
    if ref_repeat != hgvs_repeat {
        return false;
    }
    return true;
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

        let expected = vec![
            false, false, true, false, false, false, false, true, true, false
        ];
        for (i, line) in reader.lines().enumerate() {
            let line = line.unwrap();
            let line = line.trim();
            let tr: TandemRepeat = line.parse().unwrap();
            let is_correct = is_present(&tr, &sequences);
            assert_eq!(is_correct, expected[i]);
        }
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
            if is_present(&tr, &references) {
                present_count += 1;
            } else {
                // println!("{}", tr);
                // print_diff(&tr, &references);
                // println!();
            }
            max_count += 1;
        }
        println!("Present repeats: {}/{}", present_count, max_count);
    }

    fn print_diff(tr: &TandemRepeat, refs: &HashMap<String, Vec<u8>>) {
        let n = 10;
        let rflank = ref_region(refs, &tr.reference, tr.start-n, tr.start).unwrap();
        let ref_repeat = ref_region(refs, &tr.reference, tr.start, tr.end).unwrap();
        let lflank = ref_region(refs, &tr.reference, tr.end, tr.end+n).unwrap();
        println!("{} {} {}", 
            str::from_utf8(rflank).unwrap(),
            str::from_utf8(ref_repeat).unwrap(),
            str::from_utf8(lflank).unwrap()
        );
        println!("{} {} {}",
            " ".repeat(n),
            str::from_utf8(&tr.sequence()).unwrap(),
            " ".repeat(n)
        );
    }

    #[test]
    fn can_parse_hgvs() {
        let _record = "NM_01234.5:c.456-6_*22A>T";
        let _record = "NC_000017.11:g.43091687del";
        let tmp: HgvsVariant = _record.parse().unwrap();
        println!("{:?}", tmp);

        println!("{}", tmp.accession().value);
    }

    #[test]
    fn can_move_motif() {
        let motif: TandemRepeat = "SEQ1:g.6_15CG[5]".parse().unwrap();
        let flank = 5;
        let from = motif.start - flank;
        let to = motif.end + flank;
        let seq = &b"AAAAAAACGCGCGCGCGAAA"[from..to];

        let expected_motif: TandemRepeat = "SEQ1:g.8_17CG[5]".parse().unwrap();
        let corrected_motif = correct_motif(&seq, &motif, flank);

        println!(
            "{} -> {}\n{}\n{}\n{}",
            motif, corrected_motif,
            str::from_utf8(&seq).unwrap(),
            str::from_utf8(&motif.view(from, to)).unwrap(),
            str::from_utf8(&corrected_motif.view(from, to)).unwrap(),
        );

        assert_eq!(expected_motif, corrected_motif);
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

