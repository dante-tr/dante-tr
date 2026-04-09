use polars::prelude::*;

use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::io::Write;
use std::ops::Range;

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

pub(crate) fn get_n_co_occurrences(motif_df: &DataFrame, idx1: usize, seq1: String, idx2: usize, seq2: String) -> usize {
    let f3 = |s: &str| s.split(",").nth(idx1 + 1).unwrap().to_string();
    let classes = motif_df.column("module_classes").unwrap().str().unwrap().iter().map(|o| o.map(f3));
    let classes1: Column = StringChunked::from_iter_options("classes1".into(), classes).into_series().into();
    let nomenclatures = motif_df.column("module_nomenclatures").unwrap().str().unwrap().iter().map(|o| o.map(f3));
    let nomenclatures1: Column = StringChunked::from_iter_options("nomenclatures1".into(), nomenclatures).into_series().into();

    let f3 = |s: &str| s.split(",").nth(idx2 + 1).unwrap().to_string();
    let classes = motif_df.column("module_classes").unwrap().str().unwrap().iter().map(|o| o.map(f3));
    let classes2: Column = StringChunked::from_iter_options("classes2".into(), classes).into_series().into();
    let nomenclatures = motif_df.column("module_nomenclatures").unwrap().str().unwrap().iter().map(|o| o.map(f3));
    let nomenclatures2: Column = StringChunked::from_iter_options("nomenclatures2".into(), nomenclatures).into_series().into();

    let cooccurrences_df = DataFrame::new_infer_height(vec![nomenclatures1, nomenclatures2, classes1, classes2]).unwrap();

    let f = |o: Option<&str>| { let x = o.unwrap(); x == "Spanning" };
    let mask1: BooleanChunked = cooccurrences_df.column("classes1").unwrap().str().unwrap().iter().map(f).collect();
    let mask2: BooleanChunked = cooccurrences_df.column("classes2").unwrap().str().unwrap().iter().map(f).collect();
    let mask = mask1 & mask2;
    let cooccurrences_df = cooccurrences_df.filter(&mask).unwrap();

    let mut agg: DataFrame = cooccurrences_df
        .group_by(["nomenclatures1", "nomenclatures2"]).unwrap()
        .select(["nomenclatures1"])  // this is required, because otherwise count does not know how to call the new column
        .count().unwrap();
    agg.rename("nomenclatures1_count", "n_occ".into()).unwrap();
    let agg = agg.select(["n_occ", "nomenclatures1", "nomenclatures2"]).unwrap();
    let mask1: BooleanChunked = agg.column("nomenclatures1").unwrap().str().unwrap().iter()
        .map(|o: Option<&str>| { let x = o.unwrap(); x == seq1 }).collect();
    let mask2: BooleanChunked = agg.column("nomenclatures2").unwrap().str().unwrap().iter()
        .map(|o: Option<&str>| { let x = o.unwrap(); x == seq2 }).collect();

    let mask = mask1 & mask2;
    let result = agg.filter(&mask).unwrap();

    if result.height() == 0 {
        return 0;
    } else {
        let x = result.column("n_occ").unwrap().u32().unwrap().get(0).unwrap();
        return x.try_into().unwrap();
    }
}

pub(crate) fn get_nomenclature_counts(df: &DataFrame) -> DataFrame {
    let cols_in = ["lengths", "counts", "classes", "nomenclatures"];
    let cols_out = ["n_occ", "counts", "nomenclatures"];

    debug_assert!(df.get_column_names() == cols_in);

    // polars has a bit weird interface...
    let mut agg: DataFrame = df
        .group_by(cols_in).unwrap()
        .select(["nomenclatures"])  // this is required, because otherwise count does not know how to call the new column
        .count().unwrap();
    agg.rename("nomenclatures_count", "n_occ".into()).unwrap();
    let agg = agg.select(cols_out).unwrap();
    let sopt = SortMultipleOptions::new().with_order_descending(true);
    let agg = agg.sort(cols_out, sopt).unwrap();
    return agg;
}

#[allow(clippy::type_complexity)]
pub(crate) fn extract_from_df(df: &DataFrame, idx: usize) -> Result<(Vec<u64>, Vec<u64>, Vec<bool>, u64, u64), Box<dyn Error>> {
    // collect only relevant columns
    let module_df = get_module_df(df, idx)?;

    // filter only relevant rows
    let f = |o: Option<&str>| { let x = o.unwrap(); x == "Spanning" };
    let mask: BooleanChunked = module_df.column("classes")?.str()?.iter().map(f).collect();
    let spanning_df = module_df.filter(&mask)?;

    let f = |o: Option<&str>| { let x = o.unwrap(); x == "Flanking" };
    let mask: BooleanChunked = module_df.column("classes")?.str()?.iter().map(f).collect();
    let flanking_df = module_df.filter(&mask)?;

    let relevant_df = DataFrame::vstack(&spanning_df, &flanking_df)?;

    // extract to required datastructures
    let max_spanning_reps: u64 = spanning_df["counts"].u64()?.iter().max().unwrap().unwrap();
    let max_overall_reps: u64  = relevant_df["counts"].u64()?.iter().max().unwrap().unwrap();
    let counts: Vec<u64>       = relevant_df["counts"].u64()?.iter().map(|x| x.unwrap()).collect();
    let lengths: Vec<u64>      = relevant_df["lengths"].u64()?.iter().map(|x| x.unwrap()).collect();
    let is_spanning: Vec<bool> = relevant_df["classes"].str()?.iter().map(|x| x.unwrap() == "Spanning").collect();

    return Ok((counts, lengths, is_spanning, max_spanning_reps, max_overall_reps));
}

/// On success return DataFrame with columns (lengths: u64, counts: u64, classes: str,
/// nomenclatures: str)
pub(crate) fn get_module_df(df: &DataFrame, idx: usize) -> Result<DataFrame, Box<dyn Error>> {
    let f1 = |s: &str| s.len().try_into().unwrap();
    let lengths = df.column("read")?.str()?.iter().map(|o| o.map(f1));
    let lengths: Column = UInt64Chunked::from_iter_options("lengths".into(), lengths).into_series().into();

    let f2 = |s: &str| s.split(",").nth(idx).unwrap().parse().unwrap();
    let counts = df.column("module_repetitions")?.str()?.iter().map(|o| o.map(f2));
    let counts: Column = UInt64Chunked::from_iter_options("counts".into(), counts).into_series().into();

    let f3 = |s: &str| s.split(",").nth(idx).unwrap().to_string();
    let classes = df.column("module_classes")?.str()?.iter().map(|o| o.map(f3));
    let classes: Column = StringChunked::from_iter_options("classes".into(), classes).into_series().into();
    let nomenclatures = df.column("module_nomenclatures")?.str()?.iter().map(|o| o.map(f3));
    let nomenclatures: Column = StringChunked::from_iter_options("nomenclatures".into(), nomenclatures).into_series().into();

    let module_df = DataFrame::new_infer_height(vec![lengths, counts, classes, nomenclatures])?;
    return Ok(module_df);
}

pub(crate) fn get_single_string_prediction(df: &DataFrame, x: usize) -> String {
    let nomenclature = df
        .filter(&df.column("counts").unwrap().u64().unwrap().equal(x)).unwrap()
        .column("nomenclatures").unwrap()
        .str().unwrap()
        .get(0).unwrap()
        .to_string();
    return nomenclature;
}

pub(crate) fn select_spanning_reads(module_df: &DataFrame) -> DataFrame {
    // select spanning reads
    let f = |o: Option<&str>| { let x = o.unwrap(); x == "Spanning" };
    let mask: BooleanChunked = module_df.column("classes").unwrap().str().unwrap().iter().map(f).collect();
    let spanning_df = module_df.filter(&mask).unwrap();
    return spanning_df;
}

pub(crate) fn get_noms_and_occs(df: &DataFrame, a: usize) -> (Vec<String>, Vec<u32>) {
    // similar to ./../../dante_py/dante_remastr_simple.py:361:1  (use gF)
    let relevant_nomenclatures = df.filter(&df.column("counts").unwrap().u64().unwrap().equal(a)).unwrap();

    let noms = relevant_nomenclatures.column("nomenclatures").unwrap().str().unwrap().head(Some(2));
    let noms: Vec<String> = noms.iter().flatten().map(|x| x.to_string()).collect();

    let occs = relevant_nomenclatures.column("n_occ").unwrap().u32().unwrap().head(Some(2));
    let occs: Vec<u32> = occs.iter().flatten().collect();

    return (noms, occs);
}

pub(crate) fn get_nomenclatures(main_df: &DataFrame, idxs: Range<usize>) -> Vec<(u64, Vec<String>)> {
    // TODO:
    vec![
        (156, vec!["GGC[8]".to_string()]),
        (6, vec!["GGC[7]GGG[1]".to_string()]),
        (1, vec!["GGC[1]GGG[1]GGC[3]GGG[3]".to_string()]),
        (1, vec!["GGC[2]GGG[1]GGC[4]GGG[1]".to_string()]),
        (1, vec!["GGC[3]GGG[1]GGC[1]GGG[2]GGC[1]".to_string()]),
        (1, vec!["Finish me".to_string()]),
        // and 13 more
    ]
}




