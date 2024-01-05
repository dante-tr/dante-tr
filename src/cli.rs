use clap::Parser;

// Predict short tandem repeat annotation
#[derive(Parser, Debug)]
pub struct Args {
    /// Reference in FASTA format
    #[arg(short='f')]
    pub ref_file: String,

    /// Reads mapped to reference in BAM format, index (.bai) has to be present
    #[arg(short='b')]
    pub bam_file: String,

    /// Repeats in HGVS nomenclature, one per line or TSV with name and HGVS
    #[arg(short='m')]
    pub motif_file: String,

    /// Correct repeats
    #[arg(short='c')]
    pub correction: bool,

    /// Output file in TSV format.
    #[arg(short='o')]
    pub out_file: String,
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn empty_args_prints_help() {
        let args = Args::try_parse_from(["remastr"].iter()).err().unwrap();
        println!("{}", args); 
    }

    #[test]
    fn prints_help() {
        let args = Args::try_parse_from(["remastr", "-h"].iter()).err().unwrap();
        println!("{}", args);
    }

    #[test]
    fn cli_small_example() {
        let args = Args::try_parse_from([
            "remastr",
            "-f", "data/chromosomeX.fna",
            "-m", "data/mini_HGVS.txt",
            "-b", "data/mini2.bam",
            "-o", "tmp.txt"
        ].iter()).unwrap();
        println!("{:?}", args);
    }

    #[test]
    fn cli_example() {
        let args = Args::try_parse_from([
            "remastr",
            "-f", "data/chromosomeX.fna",
            "-m", "data/nomenclature_hgs_1Q_wo_names.tsv",
            "-b", "data/mini2.bam",
            "-c",
            "-o", "tmp.txt"
        ].iter()).unwrap();
        println!("{:?}", args);
    }
}

