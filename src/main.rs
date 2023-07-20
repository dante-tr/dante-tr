mod repeats;
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

#[cfg(test)]
mod tests {
    use noodles::bam as bam;
    use noodles::fasta as fasta;
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
        let mut reader = fasta::reader::Builder
            .build_from_path("data/chromosomeX.fna").unwrap();

        for result in reader.records() {
            let record = result.unwrap();

            println!("{}\t{}", record.name(), record.sequence().len());
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

