use std::io::Write;
use std::fs::{self, File};
use std::path::{Path, PathBuf};

fn main() {
    // let mut log = File::create("build.log").expect("Unable to create file");

    let mut paths: Vec<(PathBuf, PathBuf)> = Vec::new();
    collect_paths(&mut paths, "assets");
    // writeln!(log, "{:#?}", paths).unwrap();

    let out_dir = std::env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("embedded_assets.rs");
    // writeln!(log, "Writing to {}", dest_path.display()).unwrap();

    let mut out_file = File::create(&dest_path).unwrap();
    writeln!(out_file, "const EMBEDDED_FILES: [(&str, &[u8]); {}] = [", paths.len()).unwrap();
    for entry in paths {
        writeln!(
            out_file,
            "    (\"{}\", include_bytes!(\"{}\")),",
            entry.0.display(),
            entry.1.display(),
        ).unwrap();
    }
    writeln!(out_file, "];").unwrap();
}

fn collect_paths(paths: &mut Vec<(PathBuf, PathBuf)>, dir: impl AsRef<Path>) {
    for entry in fs::read_dir(dir).unwrap() {
        let entry = entry.unwrap().path();
        if entry.is_dir() {
            collect_paths(paths, entry);
        } else {
            let first = entry.clone();
            let second = entry.canonicalize().unwrap();
            paths.push((first, second));
        }
    }
}
