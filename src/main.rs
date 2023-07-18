fn main() {
    println!("Bu!");
}

#[cfg(test)]
mod tests {
    use noodles::bam as bam;
    use noodles::fasta as fasta;
    use std::fs::File;
    use std::io::BufReader;

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
        let mut reader = File::open("data/chromosomeX.fna")
            .map(BufReader::new)
            .map(fasta::Reader::new)
            .unwrap();

        for result in reader.records() {
            let record = result.unwrap();
            println!("{}\t{}", record.name(), record.sequence().len());
        }
    }
}

