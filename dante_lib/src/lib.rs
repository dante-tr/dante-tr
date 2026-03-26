mod annotation;
mod bam_index;
mod bam_ops;
mod genotyping;
mod hmm;
mod io;
mod motif_correction;
mod phasing;
mod repeats;

use rayon::prelude::*;
use std::path::Path;

pub fn run_v2(bam_file: &Path, motif_file: &Path, output: &Path, out_bam_flag: bool) {
    bam_index::check_bai(bam_file);

    let motif_records = io::read_motifs(motif_file);
    motif_records.par_iter().for_each(|motif_record| {
    // motif_records.iter().for_each(|motif_record| {

        let (left_flank, repeat, right_flank) = motif_record;
        let name = repeat.name.as_ref().unwrap().clone();
        let region_str = format!("{}:{}-{}", repeat.reference, repeat.start + 1, repeat.end);
        let region = region_str.parse().unwrap();

        let mut relevant_reads = bam_ops::RelevantReads::from(bam_file, region);
        if out_bam_flag {
            let h = relevant_reads.header();
            let out_bam_file = output.join(name.to_owned() + ".annotated.bam");
            let mut out_bam = bam_ops::init_bam(&out_bam_file.to_string_lossy(), &h);
            for record in relevant_reads.iter() {
                out_bam.write_record(&h, &record).expect("Cannot write to out bam.");
            }
            // TODO: sort bam + create bai index
        }

        // build HMM and annotate reads - polars alternative
        let modules = io::get_modules(left_flank, repeat, right_flank);
        let model = hmm::Hmm::from(&modules).log();
        let mut annotation_df /*: DataFrame */ = annotation::annotate_reads(relevant_reads.iter(), model, repeat);
        let genotyping_result = genotyping::genotype(&annotation_df, &modules);
        let phasing_results = phasing::phase(&annotation_df, &genotyping_result);

        // write results to tsv
        let out_tsv_file = output.join(name.to_owned() + ".annotations.tsv");
        annotation::print_tsv_file(&mut annotation_df, &out_tsv_file).expect("Failed writing tsv file.");

        let out_tsv_file = output.join(name.to_owned() + ".annotations.dbg.txt");
        annotation::print_dbg_file(&annotation_df, &out_tsv_file).expect("Failed writing dbg file.");

        // let out_tsv_file = output.join(name.to_owned() + ".annotations.tsv");
        // use crate::annotation::parse_tsv_file;
        // let mut tmp_df = parse_tsv_file(&out_tsv_file).expect("Err");
        // let out_tsv_file = output.join(name.to_owned() + ".annotations2.tsv");
        // print_tsv_file(&mut tmp_df, &out_tsv_file).expect("Failed writing tsv file.");

        // // write genotyping result to json
        let out_json_file = output.join(name.to_owned() + ".genotypes.json");
        let json_str = serde_json::to_string(&genotyping_result).expect("");
        io::print_to_file(&json_str, &out_json_file).expect("Failed writing json file.");

        // // write phasing result to json
        let out_json_file = output.join(name.to_owned() + ".phasing.json");
        let json_str = serde_json::to_string(&phasing_results).expect("");
        io::print_to_file(&json_str, &out_json_file).expect("Failed writing json file.");
    });

    println!("Annotation finished successfully.");
}


