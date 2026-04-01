mod reporting;

use clap::Parser;
use std::path::PathBuf;
use std::fs;

use remastr::run_v2;

// fn main() {
//     let args = Args::parse();
//     run(
//         &PathBuf::from(args.bam_file), &PathBuf::from(args.motif_file), args.output, args.out_bam,
//         (args.dedup, args.q, args.score, args.print_quality)
//     );
// }

fn main() {
    let args = ArgsNew::parse();

    // mkdir -p <output>
    fs::create_dir_all(&args.output).unwrap_or_else(|err| {
        println!("! {:?}", err.kind());
    });

    run_v2(&args.bam_file, &args.motif_file, &args.output, args.out_bam);
    println!("Finished successfully.");
    // let annotations: Vec<String> = vec!["./output/motifs/ALS.annotations.tsv".to_string()];
    // let genotypes: Vec<String>   = vec!["./output/motifs/ALS.genotypes.json".to_string()];
    // reporting::report(annotations, genotypes).unwrap();
}

// Predict short tandem repeat annotation
#[derive(Parser, Debug, PartialEq, Eq)]
struct ArgsNew {
    /// Reads mapped to reference in BAM format
    #[arg(short='b')]
    bam_file: PathBuf,

    /// Table in tsv format containing data about repeat such as HGVS nomenclature, flanks, etc.
    #[arg(short='m')]
    motif_file: PathBuf,

    /// Directory where the results of annotation will be written
    #[arg(short='o')]
    output: PathBuf,

    /// Output annotated reads in BAM.
    /// BAM files contain only reads which overlap with motif positions.
    #[arg(long="output-bams", action, verbatim_doc_comment)]
    out_bam: bool,
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

#[test]
fn test_cli_new() {
    let cmd = "dante_cli -b reads.bam -m motifs.tsv -o output_dir --output-bams";

    let args = ArgsNew::try_parse_from(cmd.split_whitespace()).unwrap();
    let result = ArgsNew {
        bam_file  : PathBuf::from("reads.bam"),
        motif_file: PathBuf::from("motifs.tsv"),
        output    : PathBuf::from("output_dir"),
        out_bam   : true,
    };
    assert_eq!(args, result);
}

