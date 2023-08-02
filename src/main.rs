use std::str;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::collections::HashMap;
use noodles::fasta as fasta;
use noodles::bam as bam;

use crate::repeats::TandemRepeat as TandemRepeat;
use crate::hmm::HMM;

mod repeats;
mod hmm;

fn main() {
    // read reference
    let references = read_reference("data/chromosomeX.fna");
    // read nomenclature
    let hgvs = File::open("data/mini_HGVS.txt").unwrap();
    let reader = BufReader::new(hgvs);

    // check nomenclature w.r.t. reference
    let mut valid_repeats = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap().trim().to_owned();
        let tr: TandemRepeat = line.parse().unwrap();
        if is_present(&tr, &references) {
            valid_repeats.push(tr);
        }
    }

    // load bam
    let mut reader = bam::indexed_reader::Builder::default()
        .build_from_path("data/mini2.bam").unwrap();
    let header = reader.read_header().unwrap();

    for repeat in valid_repeats {
        //  build HMM
        // let model = HMM::from(get_modules(repeat, reference));
        let model = HMM::default(); // TODO
        //  select relevant reads
        let tmp = format!("{}:{}-{}", repeat.reference, repeat.start+1, repeat.end);
        let region = tmp.parse().unwrap();
        let reads = reader.query(&header, &region).unwrap();

        println!("{}", repeat);
        for read in reads {
            let read = read.expect("Incorrect read.");
            // println!("{:?}", read.unwrap());
            println!("{}", read.sequence());
            println!("{}", read.quality_scores());
            let seq =  b"ACTGCA";   // TODO from read.sequence()
            let qual = b":F::F:";   // TODO from read.quality()
            let (likelihood, annotation) = model.log_predict(seq, qual);
            // postfilter
            // report()
        }
        println!();
        //  report_row()
    }
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

