use clap::Parser;

// Predict short tandem repeat annotation
#[derive(Parser, Debug, PartialEq, Eq)]
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

    /// Output annotations in tsv format.
    #[arg(short='o', verbatim_doc_comment)]
    pub output: String,

    /// Output annotated reads in BAM. The filename depends on the filename of the output.
    /// For <OUTPUT>.tsv the annotated reads will be writen to <OUTPUT>.bam
    /// BAM contains only reads which overlap with motif positions and pass the filters.
    #[arg(short='a', verbatim_doc_comment)]
    pub out_bam: bool,

    /// Correct repeats using a set of heuristics.
    /// First, the end position of a motif is adjusted in correspondence to the length of the
    /// motif. Then, the motif is checked against reference, and if it does not align, it is
    /// removed.
    #[arg(short='c', verbatim_doc_comment)]
    pub correction: bool,

    /// Filter out reads marked as PCR or optical duplicate (SAM flag 0x400)
    #[arg(short='d')]
    pub dedup: bool,

    /// Flank size
    #[arg(long="flank", default_value_t=30)]
    pub flank: usize,

    /// Minimum mapping quality to annotate
    #[arg(long="quality", default_value_t=30)]
    pub q: u8,
}

#[test]
fn print_noargs() {
    let cmd = "remastr";
    let args = Args::try_parse_from(cmd.split_whitespace()).err().unwrap();
    println!("{}", args); 
}

#[test]
fn print_help() {
    let cmd = "remastr -h";
    let args = Args::try_parse_from(cmd.split_whitespace()).err().unwrap();
    println!("{}", args);
}

#[test]
fn test_cli_minimal() {
    let cmd = "remastr -f reference.fna -b reads.bam -m motifs.tsv -o output.tsv";

    let args = Args::try_parse_from(cmd.split_whitespace()).unwrap();
    let result = Args {
        ref_file  : "reference.fna".to_string(),
        bam_file  : "reads.bam".to_string(),
        motif_file: "motifs.tsv".to_string(),
        output    : "output.tsv".to_string(),
        out_bam   : false,
        correction: false,
        dedup     : false,
        flank     : 30,
        q         : 30
    };
    assert_eq!(args, result);
}

#[test]
fn test_cli_maximal() {
    let cmd = "
        remastr
            -f reference.fna
            -b reads.bam
            -m motifs.tsv
            -o output.tsv
            --flank 40
            --quality 20
            -a -c -d
    ";

    let args = Args::try_parse_from(cmd.split_whitespace()).unwrap();
    let result = Args {
        ref_file  : "reference.fna".to_string(),
        bam_file  : "reads.bam".to_string(),
        motif_file: "motifs.tsv".to_string(),
        output    : "output.tsv".to_string(),
        out_bam   : true,
        correction: true,
        dedup     : true,
        flank     : 40,
        q         : 20
    };
    assert_eq!(args, result);
}
