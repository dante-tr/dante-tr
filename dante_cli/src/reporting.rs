use std::{error::Error, fs::File, io::Read};

use polars::prelude::*;
use minijinja::{Environment, context};
use serde_json::Value;
use std::io::Write;

pub(crate) fn report(annotations: Vec<String>, genotypes: Vec<String>) -> Result<(), Box<dyn Error>> {
    // read tsv
    let file = File::open(&annotations[0]).unwrap();
    let opts = CsvReadOptions::default().with_parse_options(CsvParseOptions::default().with_separator(b'\t'));
    let df = CsvReader::new(file).with_options(opts).finish().unwrap();

    // read json
    let mut x = File::open(&genotypes[0]).unwrap();
    let mut buf = String::new();
    let y = x.read_to_string(&mut buf).unwrap();
    let json: Value = serde_json::from_str(&buf).expect("JSON was not well-formatted");

    // collect data
    let ctx = context!(name => "John");

    // create alignment reports

    // create main report
    let mut env = Environment::new();
    // let tmp = include_str!("./../../dante_py/templates/report_template2_static.html");
    // let tmp = include_str!("./report_static.html");
    env.add_template("hello", "Hello {{ name }}!").unwrap();
    let tmpl = env.get_template("hello").unwrap();
    let result = tmpl.render(ctx).unwrap();

    // ...profit
    println!("{:?}", df);
    println!("{:?}", json);
    println!("{:?}", y);
    println!("{:?}", annotations);
    println!("{:?}", genotypes);
    println!("{}", result);
    return Ok(());
}

#[test]
fn test_reporting() {
    // cargo test report
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/ALS.annotations.tsv
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/ALS.genotypes.json
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/DM2.annotations.tsv
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/DM2.genotypes.json

    let ctx = context!(
        name => "John"
    );

    let mut env = Environment::new();
    env.add_template("report_static", include_str!("../templates/report_static.html")).unwrap();
    let tmpl = env.get_template("report_static").unwrap();
    let result = tmpl.render(ctx).unwrap();
    let mut file = File::create("report_static.html").unwrap();
    file.write_all(result.as_bytes()).unwrap();
}
