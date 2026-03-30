mod annotation;
mod bam_ops;
mod genotyping;
mod hmm;
mod io;
mod motif_correction;
mod repeats;

use rayon::prelude::*;
use std::path::Path;

pub fn run_v2(bam_file: &Path, motif_file: &Path, output: &Path, out_bam_flag: bool) {
    bam_ops::check_bai(bam_file);
    let motif_records = io::read_motifs2(motif_file);
    motif_records.par_iter().for_each(|motif_record| {
        let mut relevant_reads = bam_ops::RelevantReads::from(bam_file, &motif_record.region());
        let modules            = io::get_modules(motif_record);
        let model              = hmm::Hmm::from(&modules).log();
        let mut annotation_df  = annotation::annotate(relevant_reads.iter(), model, motif_record);
        let genotyping_result  = genotyping::genotype(&annotation_df, &modules);
        let phasing_results    = genotyping::phase(&annotation_df, &genotyping_result);

        // write results
        let name = &motif_record.name;
        let out_bam_file  = output.join(name.to_owned() + ".annotated.bam");
        let out_tsv_file  = output.join(name.to_owned() + ".annotations.tsv");
        let out_txt_file  = output.join(name.to_owned() + ".annotations.dbg.txt");
        let out_json_file = output.join(name.to_owned() + ".genotypes.json");

        if out_bam_flag { relevant_reads.write_to_file(&out_bam_file); }
        annotation::print_tsv_file(&mut annotation_df, &out_tsv_file).expect("Failed writing tsv file.");
        annotation::print_dbg_file(&annotation_df, &out_txt_file).expect("Failed writing dbg file.");
        let json_str = serde_json::to_string(&phasing_results).expect("");
        io::print_to_file(&json_str, &out_json_file).expect("Failed writing json file.");
    });
    println!("Finished successfully.");
}
