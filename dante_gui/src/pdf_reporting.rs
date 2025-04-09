use typst_as_lib::package_resolver::FileSystemCache;
use typst_as_lib::package_resolver::PackageResolver;
use typst_as_lib::typst_kit_options::TypstKitFontOptions;
use typst_as_lib::TypstEngine;

use std::path::PathBuf;

use crate::analysis_family::Data as FamilyData;
use crate::App;

use std::collections::HashMap;
use std::fs::File;
use std::fs;
use std::io::{self, BufRead, BufReader};

pub fn simple_report(data: &FamilyData) {
    let output_pdf = "typst_report.pdf";
    let typst_template = run_pdf("src/report_data/metadata.tsv", "src/report_data/STRset.tsv", &1)
    .expect("Failed to generate Typst content");        //number of samples will be counted from the number of metadata files

    // Here be dragons
    let typst_cache = App::DATA_DIR.to_string() + "/typst_cache";
    let pkg_resolver = PackageResolver::builder().cache(FileSystemCache(PathBuf::from(typst_cache))).build();

    let template = TypstEngine::builder()
        .main_file(typst_template)
        .search_fonts_with(TypstKitFontOptions::default())
        .add_file_resolver(pkg_resolver)
        .with_file_system_resolver(App::DATA_DIR)
        .build();

    let doc = template.compile().output.expect("typst::compile() returned an error!");

    let options = Default::default();

    let pdf = typst_pdf::pdf(&doc, &options).expect("Could not generate pdf.");
    std::fs::write(output_pdf, pdf).expect("Could not write pdf.");
}

fn run_pdf(metadata_path: &str, strset_path: &str, samples: &i32) -> Result<String, io::Error> {
    let base_height = 54;
    let row_height = 16;
    let rect_height = base_height + samples * row_height;
    let rect_height_minus_35 = rect_height - 35;

    let template = fs::read_to_string("src/report_data/first.txt")?;

    let mut content = template
        .replace("{rect_height}", &rect_height.to_string())
        .replace("{rect_height_minus_35}", &rect_height_minus_35.to_string());

    let metadata_content = process_metadata_to_string(metadata_path)?;
    content.push_str(&metadata_content);

    let template2 = fs::read_to_string("src/report_data/second.txt")?;
    content.push_str(&template2);

    content.push_str(&process_strset_to_string(strset_path)?);

    Ok(content)
}

fn process_metadata_to_string(metadata_path: &str) -> io::Result<String> {
    let input_file = File::open(metadata_path)?;
    let reader = BufReader::new(input_file);

    let required_headers = vec![
        "No.",
        "Sample ID",
        "Sample SI",
        "Gender",
        "Patient position in analysis",
        "Affection status",
        "Family ID",
    ];

    let mut output = String::new();

    let mut lines = reader.lines();
    if let Some(Ok(header)) = lines.next() {
        let headers: Vec<&str> = header.split('\t').collect();
        let mut col_indices: HashMap<&str, usize> = HashMap::new();

        for (i, h) in headers.iter().enumerate() {
            if required_headers.contains(h) {
                col_indices.insert(h, i);
            }
        }

        output.push_str("\n#table(\n");
        output.push_str(&format!("  columns: {},\n", required_headers.len()));
        output.push_str("  table.header(\n");
        for (i, h) in required_headers.iter().enumerate() {
            let color = "text(rgb(5, 5, 126))";
            if i == required_headers.len() - 1 {
                output.push_str(&format!("    {}[*{}*]\n", color, h));
            } else {
                output.push_str(&format!("    {}[*{}*],\n", color, h));
            }
        }
        output.push_str("  ),\n");

        // here we will iterate over each metadata file
        let lines: Vec<_> = lines.flatten().collect();
        if let Some(line) = lines.get(1) {
            let row_number = 1;
            let values: Vec<&str> = line.split('\t').collect();
            
            output.push_str("  [");
            output.push_str(&format!("{}", row_number));
            output.push_str("],");

            for h in &required_headers[1..] {
                let value = col_indices
                    .get(h)
                    .and_then(|&i| values.get(i))
                    .unwrap_or(&"");
                output.push_str(&format!("[{}],", value));
            }
            output.push('\n');
        }

        output.push_str(")\n");
    }

    Ok(output)
}

fn process_strset_to_string(strset_path: &str) -> io::Result<String> {
    let file = File::open(strset_path)?;
    let reader = BufReader::new(file);
    let mut headers: Vec<String> = Vec::new();
    let mut values: HashMap<String, String> = HashMap::new();

    let mut lines = reader.lines();
    if let Some(header_line) = lines.next() {
        headers = header_line?.split('\t').map(String::from).collect();
    } else {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Empty file"));
    }

    for line in lines {
        let row: Vec<String> = line?.split('\t').map(String::from).collect();
        if let Some(index) = headers.iter().position(|h| h == "Disease ID") {
            if row.get(index) == Some(&"DM1".to_string()) {
                for (i, value) in row.iter().enumerate() {
                    values.insert(headers[i].clone(), value.clone());
                }
                break;
            }
        }
    }

    if values.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Disease not found"));
    }

    let template = fs::read_to_string("src/report_data/third.txt")?;

    let output = template
        .replace("{DN}", values.get("Disease name").unwrap_or(&"".to_string()))
        .replace("{DI}", values.get("Disease ID").unwrap_or(&"".to_string()))
        .replace("{OI}", values.get("OMIM ID").unwrap_or(&"".to_string()))
        .replace("{MC}", values.get("Motif complexity").unwrap_or(&"".to_string()))
        .replace("{HGVS}", values.get("Clinically relevant unit (HGVS)").unwrap_or(&"".to_string()))
        .replace("{hist}", values.get("Clinically relevant unit (historical)").unwrap_or(&"".to_string()))
        .replace("{WMHGVS}", values.get("Whole motif (HGVS)").unwrap_or(&"".to_string()))
        .replace("{WMhist}", values.get("Whole motif (historical)").unwrap_or(&"".to_string()))
        .replace("{nom}", values.get("HGVS nomenclature (GRCh38)").unwrap_or(&"".to_string()))
        .replace("{MM}", values.get("Molecular mechanism").unwrap_or(&"".to_string()))
        .replace("{N}", values.get("Notes").unwrap_or(&"".to_string()))
        .replace("{C}", values.get("Citation (references)").unwrap_or(&"".to_string()))
        .replace("{G}", values.get("Gene").unwrap_or(&"".to_string()))
        .replace("{GA}", values.get("Gene abbreviation").unwrap_or(&"".to_string()))
        .replace("{I}", values.get("Inheritance").unwrap_or(&"".to_string()))
        .replace("{PhR}", values.get("Physiological range").unwrap_or(&"".to_string()))
        .replace("{PrR}", values.get("Premutation range").unwrap_or(&"".to_string()))
        .replace("{PaR}", values.get("Pathogenic range").unwrap_or(&"".to_string()))
        .replace("{GZR}", values.get("Grey-zone range").unwrap_or(&"".to_string()))
        .replace("{Chr}", values.get("Chromosome").unwrap_or(&"".to_string()))
        .replace("{GC}", values.get("Gene context").unwrap_or(&"".to_string()))
        .replace("{PC}", values.get("Protein context").unwrap_or(&"".to_string()));
    
    Ok(output)
}