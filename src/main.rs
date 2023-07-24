use std::str;

mod repeats;
use std::collections::HashMap;
use noodles::fasta as fasta;

use repeats::TandemRepeat as TandemRepeat;

fn main() {
    // read nomenclature
    // read reference
    // check nomenclature w.r.t. reference
    // load bam
    // for all nomenclatures:
    //     build HMM
    //     reads = bam.query()
    //     for read in reads:
    //         prob, annotation = HMM.annotate(read)
    //         postfilter
    //         report()
    //     report_row()
    println!("Bu!");
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

fn check_repeat(tr: &TandemRepeat, refseq: &HashMap<String, Vec<u8>>) -> bool {
    let seq = match refseq.get(&tr.reference) {
        None => { return false; },
        Some(x) => { x },
    };
    let seq_repeat1 = &seq[tr.start..tr.end];
    let seq_repeat2 = &tr.sequence();
    if seq_repeat1 != seq_repeat2 {
        println!("{}", str::from_utf8(&seq[tr.start-10..tr.end+10]).unwrap());
        println!("{}", str::from_utf8(seq_repeat2).unwrap());
        return false;
    }
    println!("{}", str::from_utf8(seq_repeat1).unwrap());
    println!("{}", str::from_utf8(seq_repeat2).unwrap());
    return true;
}

#[cfg(test)]
mod tests {
    use noodles::bam as bam;
    use hgvs::parser::HgvsVariant;
    use std::fs::File;
    use std::io::{prelude::*, BufReader};
    use super::*;

    #[test]
    fn can_load_bam() {
        let mut reader = bam::indexed_reader::Builder::default()
            .build_from_path("data/mini.bam").unwrap();

        let header = reader.read_header().unwrap();

        // let region = "sq0:5-8".parse().unwrap();
        // let query = reader.query(&header, &region).unwrap();

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
            let is_correct = check_repeat(&tr, &sequences);
            assert_eq!(is_correct, expected[i]);
        }
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

