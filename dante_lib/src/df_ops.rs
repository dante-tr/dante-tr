use polars::prelude::*;

use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::io::Write;

pub(crate) fn print_tsv_file(df: &DataFrame, p: &Path) -> Result<(), Box<dyn Error>> {
    // The python part ( ./../../dante_py/src_new/constants.py:13:1 ) needs these columns:
    // MOTIF_COLUMN_ID = "name"
    // MOTIF_COLUMN_NAME = "motif"
    // MOTIF_COLUMN_READ_ID = "read_id"
    // MOTIF_COLUMN_MODULES = "modules"
    // MOTIF_COLUMN_N_MODS = "n_modules"
    // MOTIF_COLUMN_MOD_CLASS = "module_classes"
    // MOTIF_COLUMN_MISMATCHES_STR = "mismatches_str"
    // MOTIF_COLUMN_MODULE_REPETITIONS = "module_repetitions"
    // MOTIF_COLUMN_MODULE_NOMENCLATURES = "module_nomenclatures"
    // MOTIF_COLUMN_MODULE_SEQUENCES = "module_sequences"

    let required_cols = [
        "name", "motif", "read_id", "modules", "mismatches_str", "n_modules",
        "module_sequences", "module_nomenclatures", "module_repetitions", "module_classes",
    ];
    let mut tmp = df.select(required_cols)?;
    let file = File::create(p)?;
    CsvWriter::new(file).with_separator(b'\t').finish(&mut tmp)?;
    return Ok(());
}

pub(crate) fn print_dbg_txt_file(df: &DataFrame, p: &Path) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(p)?;
    // use polars::frame::row::Row;
    // let mut row = Row::default();
    for i in 0..df.height() {
        let row = df.get_row(i)?.0;
        let col_names = df.columns();
        for (name, value) in std::iter::zip(col_names, row) {
            writeln!(file, "{}\t{}", name.name(), value)?;
        }
        writeln!(file)?;
    }
    return Ok(());
}

pub(crate) fn print_dbg_tsv_file(df: &DataFrame, p: &Path) -> Result<(), Box<dyn Error>> {
    let file = File::create(p)?;
    CsvWriter::new(file).with_separator(b'\t').finish(&mut df.clone())?;
    return Ok(());
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn construct_df(
    names: Vec<String>,
    motifs: Vec<String>,
    read_sns: Vec<u64>,
    read_ids: Vec<String>,
    mate_orders: Vec<String>,
    qualities: Vec<String>,
    reads: Vec<String>,
    references: Vec<String>,
    moduleses: Vec<String>,
    mismatches_strs: Vec<String>,
    log_likelihoods: Vec<f32>,
    left_bgs: Vec<u64>,
    right_bgs: Vec<u64>,
    n_deletionses: Vec<u64>,
    n_insertionses: Vec<u64>,
    n_mismatcheses: Vec<u64>,
    n_moduleses: Vec<u64>,
    module_baseses: Vec<String>,
    module_repetitionses: Vec<String>,
    module_sequenceses: Vec<String>,
    module_nomenclatureses: Vec<String>,
    module_classeses: Vec<String>
) -> Result<DataFrame, Box<dyn Error>> {
    let result = df![
        "name"                 => names,                    // "ALS"
        "motif"                => motifs,                   // "chr15:g.22786680_22786703GGC[8]"
                                                            // motif modules
        "read_sn"              => read_sns,                 // 0
        "read_id"              => read_ids,                 // "HISEQ1:29:HA2WPADXX:2:2202:2985:13224"
        "mate_order"           => mate_orders,              // "1"
                                                            // TODO: add seq?
        "quality"              => qualities,                // "No quality" TODO: add qual?
        "read"                 => reads,                    // "CCTCTTCCTGCTCCTCCCCCACCCGTCCCCCTCCCCTCCCCCGCCCGCGCCTCCCGGTCACCCCCCATCCCGCCCCGCGGGGCGCGGCGCGCAGGCGCAGGCTCGGAGGGCGGGCGCGGGCGGAATGGGGACTGCAGCTGCGGCAGCG"
        "reference"            => references,               // "---------------------------------------------------------------------------------------------------------------------GGGCGGAATGGGGACTGCAGCTGCGGCAGCG"
        "modules"              => moduleses,                // "---------------------------------------------------------------------------------------------------------------------0000000000000000000000000000001"
        "mismatches_str"       => mismatches_strs,          // "____________________________________________________________________________________________________________________________________________________"
        "log_likelihood"       => log_likelihoods,          // -172.339767
        "left_bg"              => left_bgs,                 // 117
        "right_bg"             => right_bgs,                // 0
        "n_deletions"          => n_deletionses,            // 0
        "n_insertions"         => n_insertionses,           // 0
        "n_mismatches"         => n_mismatcheses,           // 0
        "n_modules"            => n_moduleses,              // 3
        "module_bases"         => module_baseses,           // "30,1,0"
        "module_repetitions"   => module_repetitionses,     // "1,0,0"
        "module_sequences"     => module_sequenceses,       // "GGGCGGAATGGGGACTGCAGCTGCGGCAGC,G,"
        "module_nomenclatures" => module_nomenclatureses,   // "GGGCGGAATGGGGACTGCAGCTGCGGCAGC[1],G[1],"
        "module_classes"       => module_classeses,         // "Flanking,Missing,Missing"
    ]?;
    return Ok(result);
}

// #[cfg(test)]
// pub fn parse_tsv_file(p: &Path) -> Result<DataFrame, Box<dyn Error>> {
//     let file = File::open(p)?;
//     let opts = CsvReadOptions::default().with_parse_options(CsvParseOptions::default().with_separator(b'\t'));
//     let df = CsvReader::new(file).with_options(opts).finish()?;
//     return Ok(df);
// }


