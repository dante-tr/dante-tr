use std::{error::Error, fs::File, io::Read};

use polars::prelude::*;
use minijinja::{Environment, context};
use serde_json::Value;

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
    let tmp = include_str!("./../../dante_py/templates/report_template2_static.html");
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
