use typst_as_lib::typst_kit_options::TypstKitFontOptions;
use typst_as_lib::TypstEngine;
use typst_as_lib::package_resolver::PackageResolver;
use typst_as_lib::package_resolver::FileSystemCache;

use std::path::PathBuf;

use crate::analysis_family::Data as FamilyData;
use crate::App;

use std::fs::File;
use std::io::{self, BufRead, Write, BufReader};
use std::collections::HashMap;
use std::error::Error;
use std::process::Command;


pub(super) fn simple_report(data: &FamilyData) {
    let output_pdf = "typst_report.pdf";
    let typst_template = work_on_me();

    // Here be dragons
    let typst_cache = App::DATA_DIR.to_string() + "/typst_cache";
    let pkg_resolver = PackageResolver::builder()
        .cache(FileSystemCache(PathBuf::from(typst_cache)))
        .build();

    let template = TypstEngine::builder()
        .main_file(typst_template)
        .search_fonts_with(TypstKitFontOptions::default())
        .add_file_resolver(pkg_resolver)
        .with_file_system_resolver(App::DATA_DIR)
        .build();

    let doc = template.compile().output
        .expect("typst::compile() returned an error!");

    let options = Default::default();

    let pdf = typst_pdf::pdf(&doc, &options).expect("Could not generate pdf.");
    std::fs::write(output_pdf, pdf).expect("Could not write pdf.");
}

fn work_on_me() -> String {
    run_pdf().unwrap();
    return std::fs::read_to_string("src/report_data/report_result.typ")
        .expect("Failed to read the generated Typst template");
}

fn run_pdf() -> Result<(), Box<dyn Error>> {
    println!("idem");
    let output_filename = "src/report_data/report_result.typ";
    let mut output_file = File::create(output_filename)?;
    println!("1");
    let metadata_rows = get_metadata_row_count("src/report_data/metadata.tsv")?;
    println!("2");
    let base_height = 70 - 16;
    let row_height = 16;
    let rect_height = base_height + metadata_rows * row_height;
    
    writeln!(
        output_file,
        "#set page(margin: 15mm) // A4 is default\n\
        \n#image(\"logo.png\", width: 180mm)\n\
        #align(right + top)[#text(20pt, strong[RESULTS REPORT])]\n\
        #align(right + top)[*Report ID:* 2025022]\n\
        \n#place(right, dy: -8pt, dx: 9pt, rect(width: 190mm, height: {}pt, fill: rgb(195, 215, 255), radius: 15%))\n\
        #place(right, dy: 20pt, rect(width: 182mm, height: {}pt, fill: rgb(255, 255, 255), radius: 15%))\n\
        #block[\n  #text(14pt, rgb(7, 7, 87), strong[Sample information])]\n\
        \n#set text(size: 10pt)\n\
        #set table(\n  stroke: none,\n  align: (x, y) => (\n    if x > 0 {{ center }}\n    else {{ left }}\n  )\n)",
        rect_height, rect_height - 35
    )?;
    
    process_metadata(&mut output_file)?;

    let content = r#"
#place(right, dy: 10pt, dx: 9pt, rect(width: 190mm, height: 320pt, fill: rgb(195, 215, 255), radius: 5%))
#place(right, dy: 38pt, dx: -2pt, rect(width: 182mm, height: 285pt, fill: rgb(255, 255, 255), radius:5%))
#" "
#block[#text(14pt, rgb(7, 7, 87), strong[Target information])]
#set text(size: 8pt)
"#;

    writeln!(output_file, "{}", content)?;

    process_strset(&mut output_file)?;
    Ok(())
}


fn get_metadata_row_count(filename: &str) -> io::Result<usize> {
    let input_file = File::open(filename)?;
    let reader = BufReader::new(input_file);
    Ok(reader.lines().count() - 1)
}


fn process_metadata(output_file: &mut File) -> io::Result<usize> {
    let input_filename = "src/report_data/metadata.tsv";
    println!("3");
    let input_file = File::open(input_filename)?;
    let reader = BufReader::new(input_file);
    let mut lines = reader.lines();

    let required_headers = vec![
        "No.", "Sample ID", "Sample SI", "Gender", 
        "Patient position in analysis", "Affection status", "Family ID"
    ];

    let mut row_count = 0;

    if let Some(Ok(header)) = lines.next() {
        let headers: Vec<&str> = header.split('\t').collect();
        let mut col_indices: HashMap<&str, usize> = HashMap::new();
        
        for (i, h) in headers.iter().enumerate() {
            if required_headers.contains(h) {
                col_indices.insert(h, i);
            }
        }

        writeln!(output_file, "\n#table(")?;
        writeln!(output_file, "  columns: {},", required_headers.len())?;
        writeln!(output_file, "  table.header(")?;
        for (i, h) in required_headers.iter().enumerate() {
            let color = "text(rgb(5, 5, 126))";
            if i == required_headers.len() - 1 {
                writeln!(output_file, "    {}[*{}*]", color, h)?;
            } else {
                writeln!(output_file, "    {}[*{}*],", color, h)?;
            }
        }
        writeln!(output_file, "  ),")?;

        let mut row_number = 1;
        for line in lines.flatten() {
            let values: Vec<&str> = line.split('\t').collect();
            write!(output_file, "  [{}],", row_number)?;
            for h in &required_headers[1..] {
                let value = col_indices.get(h).and_then(|&i| values.get(i)).unwrap_or(&"");
                write!(output_file, "[{}],", value)?;
            }
            writeln!(output_file)?;
            row_number += 1;
            row_count += 1;
        }

        writeln!(output_file, ")")?;
    }

    Ok(row_count)
}

fn process_strset(output_file: &mut File) -> Result<(), Box<dyn Error>> {
    let file = File::open("src/report_data/STRset.tsv")?;
    println!("4");
    let reader = BufReader::new(file);
    let mut headers: Vec<String> = Vec::new();
    let mut values: HashMap<String, String> = HashMap::new();

    let mut lines = reader.lines();
    if let Some(header_line) = lines.next() {
        headers = header_line?.split('\t').map(String::from).collect();
    } else {
        return Err("Empty file".into());
    }
    
    for line in lines {
        let row: Vec<String> = line?.split('\t').map(String::from).collect();
        if let Some(index) = headers.iter().position(|h| h == "Disease ID") {           // premysliet
            if row.get(index) == Some(&"DM1".to_string()) {
                for (i, value) in row.iter().enumerate() {
                    values.insert(headers[i].clone(), value.clone());
                }
                break;
            }
        }
    }
    
    if values.is_empty() {
        return Err("Disease not found".into());
    }
    
    writeln!(
        output_file,
        "#table(
  columns: 3,
  [#table(
    align: left,
    columns: 2,
    stroke: none,
    text(rgb(5, 5, 126))[*Disease*], [*{}*],
    text(rgb(5, 5, 126))[Disease abbreviation], [{}],
    text(rgb(5, 5, 126))[OMIM ID], [#{}],
    text(rgb(5, 5, 126))[Motif complexity], [{}],
    text(rgb(5, 5, 126))[*Clinically relevant unit (HGVS)*], [*{}*],
    text(rgb(5, 5, 126))[Clinically relevant unit (historical)], [{}],
    text(rgb(5, 5, 126))[Whole motif (HGVS)], [{}],
    text(rgb(5, 5, 126))[Whole motif (historical)], [{}],
    text(rgb(5, 5, 126))[HGVS nomenclature (GRCh38)], [{}],
    text(rgb(5, 5, 126))[Molecular mechanism], [{}],
    text(rgb(5, 5, 126))[Motif - Notes], [{}],
    text(rgb(5, 5, 126))[Citation (references)], [{}],
  )], [#table(
    align: left,
    columns: 2,
    stroke: none,
    text(rgb(5, 5, 126))[Gene], [{}],
    text(rgb(5, 5, 126))[Gene abbreviation], [{}],
    text(rgb(5, 5, 126))[Inheritance], [{}],
    text(rgb(5, 5, 126))[Physiological range], [{}],
    text(rgb(5, 5, 126))[Premutation range], [{}],
    text(rgb(5, 5, 126))[Pathogenic range], [{}],
    text(rgb(5, 5, 126))[Grey-zone range], [{}],
  )], [#table(
    align: left,
    columns: 2,
    stroke: none,
    text(rgb(5, 5, 126))[Chromosome], [{}],
    text(rgb(5, 5, 126))[Gene context], [{}],
    text(rgb(5, 5, 126))[Protein context], [{}],
  )],
)",
        values.get("Disease name").unwrap_or(&"".to_string()),
        values.get("Disease ID").unwrap_or(&"".to_string()),
        values.get("OMIM ID").unwrap_or(&"".to_string()),
        values.get("Motif complexity").unwrap_or(&"".to_string()),
        values.get("Clinically relevant unit (HGVS)").unwrap_or(&"".to_string()),
        values.get("Clinically relevant unit (historical)").unwrap_or(&"".to_string()),
        values.get("Whole motif (HGVS)").unwrap_or(&"".to_string()),
        values.get("Whole motif (historical)").unwrap_or(&"".to_string()),
        values.get("HGVS nomenclature (GRCh38)").unwrap_or(&"".to_string()),
        values.get("Molecular mechanism").unwrap_or(&"".to_string()),
        values.get("Notes").unwrap_or(&"".to_string()),
        values.get("Citation (references)").unwrap_or(&"".to_string()),
        values.get("Gene").unwrap_or(&"".to_string()),
        values.get("Gene abbreviation").unwrap_or(&"".to_string()),
        values.get("Inheritance").unwrap_or(&"".to_string()),
        values.get("Physiological range").unwrap_or(&"".to_string()),
        values.get("Premutation range").unwrap_or(&"".to_string()),
        values.get("Pathogenic range").unwrap_or(&"".to_string()),
        values.get("Grey-zone range").unwrap_or(&"".to_string()),
        values.get("Chromosome").unwrap_or(&"".to_string()),
        values.get("Gene context").unwrap_or(&"".to_string()),
        values.get("Protein context").unwrap_or(&"".to_string()),
    );

    Ok(())
}