use clap::Parser;
use noodles::bam as bam;
use noodles::fasta as fasta;
use noodles::sam::record::quality_scores::Score;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::str;

use crate::hmm::HMM;
use crate::hmm::Module;
use crate::repeats::TandemRepeat;

mod hmm;
mod repeats;

// Predict short tandem repeat annotation
#[derive(Parser, Debug)]
struct Args {
    /// Reference fasta file
    #[arg(short, long)]
    fasta: String,
    /// HGVS nomenclature, one per line
    #[arg(short, long)]
    nomenclature: String,
    /// BAM file, BAI index have to be present
    #[arg(short, long)]
    bam: String,
}

fn main() {
    let args = Args::parse();
    let references = read_reference(&args.fasta);
    let repeats = read_nomenclature(&args.nomenclature);
    let bam_refs = read_bam_refs(&args.bam);

    let mut valid_repeats = Vec::new();
    for repeat in repeats {
        if is_present(&repeat, &references) {
            valid_repeats.push(repeat);
        }
    }

    valid_repeats.par_iter().for_each(|repeat| {
        let mut reader = bam::indexed_reader::Builder::default()
            .build_from_path(&args.bam).unwrap();
        let header = reader.read_header().unwrap();

        let modules = get_modules(&repeat, &references, 20);
        let model = HMM::from(&modules).log();

        let region = format!("{}:{}-{}", repeat.reference, repeat.start+1, repeat.end).parse().unwrap();
        let reads = reader.query(&header, &region).unwrap();

        for read in reads {
            let read = read.expect("Incorrect read.");
            let seq: Vec<_> = read.sequence().as_ref().iter().map(|&x| x.into()).collect();
            let qual: Vec<_> = read.quality_scores().as_ref().iter().map(|&x| remap(x)).collect();
            let (likelihood, annotation) = model.log_predict(&seq, &qual);

            let reconstructed_reference = model.reconstruct_sequence(&annotation);
            let reconstructed_read = model.realign_read(&annotation, &seq); 
            let mods = model.reconstruct_mod_ids(&annotation);

            println!(">{} {} {}\n{}\n{}\n{}", 
                read.read_name().unwrap(), repeat, likelihood,
                str::from_utf8(&reconstructed_read).unwrap(),
                str::from_utf8(&reconstructed_reference).unwrap(),
                str::from_utf8(&mods).unwrap()
            );
        }
    })
}

fn remap(x: Score) -> u8 {
    let c: char = x.into();
    return c as u8;
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
    for i in 0..repeat.copy_unit.len() {
        modules.push((&repeat.copy_unit[i][..], repeat.copy_number[i]).into());
    }
    modules.push(right_flank.into());
    return modules;
}

#[cfg(test)]
mod tests {
    use hgvs::parser::HgvsVariant;
    use std::fs::File;
    use std::io::{prelude::*, BufReader};
    use super::*;

    #[test]
    fn can_get_sequence_id_from_bam() {
        let filename: &str = "data/mini2.bam";
        let file = File::open(filename).unwrap();
        let index = bam::bai::read(filename.to_owned() + ".bai").unwrap();
        let mut reader = bam::IndexedReader::new(file.try_clone().unwrap(), index.clone());

        let header = reader.read_header().unwrap();
        let seqs = header.reference_sequences();
        for s in seqs.iter() {
            let name = s.0.to_string();
            let length = s.1.length().get();
            println!("{} {}", name, length);
        }
    }

    #[test]
    fn can_get_sequence_id_from_fasta() {
        let filename = "data/chromosomeX.fna";
        let file = File::open(filename).unwrap();
        let mut reader = fasta::Reader::new(BufReader::new(file));

        for record in reader.records() {
            let record = record.unwrap();
            let name = record.name().to_string();
            let length = record.sequence().len();
            println!("{} {}", name, length);
        }
    }

    #[test]
    fn test_reference_checking() {
        let references = read_reference("data/chromosomeX.fna");
        let repeats = read_nomenclature("data/mini_HGVS.txt");
        let bam_refs = read_bam_refs("data/mini2.bam");
        println!("{:?}", references.keys());
        println!("{:?}", repeats);
        println!("{:?}", bam_refs);
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
                println!("{}", tr);
                print_diff(&tr, &references);
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

