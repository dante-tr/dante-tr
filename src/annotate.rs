use noodles::fasta;
use std::str;
use clap::Parser;

mod hmm;
mod io;
mod repeats;

use crate::hmm::{Hmm, Module};
use crate::io::{get_modules, read_reference};
use crate::repeats::TandemRepeat;

// Predict short tandem repeat annotation
#[derive(Parser, Debug, PartialEq, Eq)]
pub struct Args {
    /// Reference in FASTA format
    #[arg(short='f')]
    pub ref_file: String,

    /// Reads in FASTA format
    #[arg(short='r')]
    pub read_file: String,

    /// HGVS nomenclature
    #[arg(short='n')]
    pub nomenclature: String,

    /// Flank size
    #[arg(long="flank", default_value_t=30)]
    pub flank: usize,
}


fn main() {
    let args = Args::parse();

    // let reference = "data/real/grch38_decoy.fa";
    let reference = args.ref_file;
    // let reads = "data/2024-03_test/wrong_reads_long.fasta";
    // let reads = "data/2024-03_test/wrong_reads.fasta";
    let reads = args.read_file;
    // let repeat = "chr3:g.129172580_129172733GCAG[20]ACAG[9]AC[19]";
    let repeat = args.nomenclature;

    let references = read_reference(&reference);
    let repeat: TandemRepeat = repeat.parse().unwrap();

    let modules = get_modules(&repeat, &references, args.flank);
    let model = Hmm::from(&modules).log();

    let mut reader = fasta::reader::Builder.build_from_path(reads).unwrap();
    for record in reader.records() {
        let record = record.unwrap();

        let seq = record.sequence().as_ref().to_vec();
        let qual = vec![b'I'; seq.len()];

        let (likelihood, annotation) = model.log_predict(&seq, &qual);
        let (new_annot, reconstructed_read) = model.realign(&annotation, &seq);

        let reconstructed_reference = model.reconstruct_sequence(&new_annot);
        let mods = model.reconstruct_mod_ids(&new_annot);
        println!(
            "{}\n{}\n{}\n{}",
            likelihood,
            str::from_utf8(&reconstructed_read).unwrap(),
            str::from_utf8(&reconstructed_reference).unwrap(),
            str::from_utf8(&mods).unwrap()
        );
    }
}
