use std::{error::Error, fs::File, io::Read};
use std::path::Path;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Write;

use minijinja::{context, Environment};
use serde_json::Value;

pub(crate) fn report(args: &crate::ArgsNew, motifs_dir: &Path) -> Result<(), Box<dyn Error>> {
    let version = env!("CARGO_PKG_VERSION");

    let mut motifs = Vec::new();
    let mut buf = String::new();
    let motif_names: Vec<String> = read_motifs(&args.motif_file).unwrap();
    for motif in motif_names {
        let filename = motifs_dir.join(motif.to_owned() + ".motif.json");
        let mut fp = File::open(filename).unwrap();
        let _ = fp.read_to_string(&mut buf).unwrap();
        let json: Value = serde_json::from_str(&buf).expect("JSON was not well-formatted");
        motifs.push(json);
        buf.clear();
    };

    let data = context! {
        dante_version => version,
        dante_params => context! {
            file_bam => args.bam_file,
            file_motif => args.motif_file,
            max_noms => 5,
            // is_male => true,        // TODO
        },
        motifs => motifs
    };

    let mut env = Environment::new();
    env.add_template("report_static", include_str!("../templates/report_static.html")).unwrap();
    env.add_filter("percent", |value: f64| { format!("{:.0}%", value * 100.0) });
    let tmpl = env.get_template("report_static").unwrap();
    let result = tmpl.render(context! { data => data }).unwrap();
    let output_file = args.output.join("report_static.html");
    let mut file = File::create(output_file).unwrap();
    file.write_all(result.as_bytes()).unwrap();

    Ok(())
}

fn read_motifs(filename: &Path) -> Result<Vec<String>, Box<dyn Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut result = Vec::new();
    for line in reader.lines().skip(1) {
        let line = line?.trim().to_owned();
        let name = *line.split('\t').collect::<Vec<_>>().first().ok_or("0 columns")?;
        let name = name.to_string();
        result.push(name);
    }
    return Ok(result);
}
