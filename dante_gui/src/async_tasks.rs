use std::path::PathBuf;
use std::process::Command;

#[cfg(not(feature = "local"))]
use crate::App;

pub async fn run_annotation(motif_file: PathBuf, bam_file: PathBuf, output_file: PathBuf) -> String {
    // optional params
    let out_bam = false;
    let dedup = false;
    let print_quality = false;
    let q = 30;
    let score: Option<char> = None;

    let out = output_file.to_string_lossy().to_string();
    remastr::run(
        &bam_file, &motif_file, out,
        out_bam, (dedup, q, score, print_quality)
    );

    let mut sample = bam_file.clone();
    sample.set_extension("");
    let sample = sample.file_name().unwrap().to_string_lossy();
    return format!("Annotation of sample {} successful.", sample);
}

pub async fn run_genotyping(annotation_file: PathBuf, dante_output_dir: PathBuf) -> String {
    let output_log;

    #[cfg(feature = "local")] {
    let script = "./../../dante/dante_remastr_standalone.py".to_string();
    let mut cmd = Command::new("python");
    let cmd = cmd.arg(&script)
        .arg("--input-tsv").arg(&annotation_file)
        .arg("--output-dir").arg(&dante_output_dir)
        .arg("--verbose");
    println!("{:?}", cmd);
    output_log = cmd.output()
        .expect("failed to run python part of Dante locally");
    }

    #[cfg(not(feature = "local"))] {
    let bin = format!("{}/dante_remastr_standalone", App::DATA_DIR);
    output_log = Command::new(bin)
        .arg("--input-tsv").arg(&annotation_file)
        .arg("--output-dir").arg(&dante_output_dir)
        .arg("--verbose")
        .output()
        .expect("failed to run python part of Dante");
    }

    println!("{:?}", output_log);

    return format!("Genotyping of sample {} successful.", dante_output_dir.file_name().unwrap().to_string_lossy());
}
