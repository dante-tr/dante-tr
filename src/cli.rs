use clap::Parser;

// Predict short tandem repeat annotation
#[derive(Parser, Debug)]
pub struct Args {
    /// Reference in FASTA format
    #[arg(short='f')]
    pub ref_file: String,

    /// Reads mapped to reference in BAM format
    #[arg(short='b')]
    pub bam_file: String,

    /// Repeats in HGVS nomenclature, one per line or TSV with name and HGVS
    #[arg(short='m')]
    pub motif_file: String,

    /// Output prefix. remaSTR outputs annotations (<OUTPUT>.tsv) and annotated reads (<OUTPUT>.bam).
    /// BAM contains only reads which overlap with motif positions and pass the filters.
    #[arg(short='o', verbatim_doc_comment)]
    pub output: String,

    /// Correct repeats
    #[arg(short='c')]
    pub correction: bool,

    /// Flank size
    #[arg(long="flank", default_value_t=30)]
    pub flank: usize,

    /// Minimum mapping quality to annotate
    #[arg(long="quality", default_value_t=30)]
    pub q: u8,
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
            "-o", "tmp"
        ].iter()).unwrap();
        println!("{:?}", args);
        assert!(!args.correction);
    }

    #[test]
    fn cli_example() {
        let args = Args::try_parse_from([
            "remastr",
            "-f", "data/chromosomeX.fna",
            "-m", "data/nomenclature_hgs_1Q_wo_names.tsv",
            "-b", "data/mini2.bam",
            "-c",
            "-o", "tmp"
        ].iter()).unwrap();
        println!("{:?}", args);
        assert!(args.correction);
    }
}

