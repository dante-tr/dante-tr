use clap::Parser;
use std::path::PathBuf;

use remastr::run;

fn main() {
    let args = Args::parse();
    run(
        &PathBuf::from(args.bam_file), &PathBuf::from(args.motif_file), args.output, args.out_bam,
        (args.dedup, args.q, args.score, args.print_quality)
    );
}

// Predict short tandem repeat annotation
#[derive(Parser, Debug, PartialEq, Eq)]
struct Args {
    /// Reads mapped to reference in BAM format
    #[arg(short='b')]
    bam_file: String,

    /// Repeats in HGVS nomenclature
    #[arg(short='m')]
    motif_file: String,

    /// Output annotations in tsv format.
    #[arg(short='o', verbatim_doc_comment)]
    output: String,

    /// Output annotated reads in BAM. The filename depends on the filename of the output.
    /// For <OUTPUT>.tsv the annotated reads will be writen to <OUTPUT>.bam
    /// BAM contains only reads which overlap with motif positions and pass the filters.
    #[arg(short='a', verbatim_doc_comment)]
    out_bam: bool,

    /// Filter out reads marked as PCR or optical duplicate (SAM flag 0x400)
    #[arg(short='d')]
    dedup: bool,

    /// Minimum mapping quality to annotate
    #[arg(long="quality", default_value_t=30)]
    q: u8,

    /// Quality score used for reads
    #[arg(short='s')]
    score: Option<char>,

    /// Print quality scores. Only used for debugging
    #[arg(long="print-quality", action)]
    print_quality: bool,
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
    let cmd = "remastr -b reads.bam -m motifs.tsv -o output.tsv";

    let args = Args::try_parse_from(cmd.split_whitespace()).unwrap();
    let result = Args {
        bam_file  : "reads.bam".to_string(),
        motif_file: "motifs.tsv".to_string(),
        output    : "output.tsv".to_string(),
        out_bam   : false,
        dedup     : false,
        q         : 30,
        score     : None,
        print_quality : false
    };
    assert_eq!(args, result);
}

#[test]
fn test_cli_maximal() {
    let cmd = "
        remastr
            -b reads.bam
            -m motifs.tsv
            -o output.tsv
            --quality 20
            -a -d
    ";

    let args = Args::try_parse_from(cmd.split_whitespace()).unwrap();
    let result = Args {
        bam_file  : "reads.bam".to_string(),
        motif_file: "motifs.tsv".to_string(),
        output    : "output.tsv".to_string(),
        out_bam   : true,
        dedup     : true,
        q         : 20,
        score     : None,
        print_quality : false
    };
    assert_eq!(args, result);
}
