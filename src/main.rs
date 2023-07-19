fn main() {
    println!("Bu!");
}

#[cfg(test)]
mod tests {
    use noodles::bam as bam;
    use noodles::fasta as fasta;
    use hgvs::parser::HgvsVariant;

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
        let _record = b"NC_000023.11:g.2789717_2789870ATTTT[30]";
        let _record = "NM_000044.3:g.123_191CAG[25]";
        // let _record = "NM_01234.5:c.456-6_*22A>T";
        // let _record = "NC_000017.11:g.43091687del";
        let tmp: HgvsVariant = _record.parse().unwrap();
        println!("{:?}", tmp);

        println!("{}", tmp.accession().value);
        // println!("{}", tmp.loc_edit().loc);
    }
}

