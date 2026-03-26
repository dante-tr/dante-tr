#![allow(dead_code)]

use itertools::izip;
use std::error::Error;

use polars::prelude::*;
use statrs::{distribution::{Binomial, Discrete}, statistics::Statistics};
use ndarray::{self, s, Array};
use serde::{Serialize, Deserialize};

use crate::hmm::Module;

// #[test]
// fn test_genotyping_from_dataframe() {
//     use std::path::PathBuf;
//     use crate::annotation::parse_tsv_file;
// 
//     // /home/balaz/data/projects/STRs3/tools/remastr_dev/dante_lib/DM2.annotations.tsv
//     let tsv_file = PathBuf::from("/home/balaz/data/projects/STRs3/tools/remastr_dev/dante_lib/DM2.annotations.tsv");
//     let df: DataFrame = parse_tsv_file(&tsv_file).unwrap();
//     // let result = genotype(&df);
//     // let json = serde_json::to_string(&result).unwrap();
//     // println!("{}", json);
// }

#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct GenotypingResults {
    pub(crate) modules: Vec<ModuleResult>
}

#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct ModuleResult {
    pub(crate) predictions_enum: (Prediction, Prediction),
    pub(crate) predictions_seq: (String, String),
    confidences: [f64; 7],
    likelihoods: ndarray::Array2<f64>,
}

#[derive(Debug, Serialize, Deserialize, Clone, Copy)]
pub(crate) enum Prediction {
    Num(usize),     // change breaks python parsing
    Expansion,      // change breaks python parsing
    Background
}

pub(crate) fn genotype(df: &DataFrame, modules: &[Module]) -> GenotypingResults {
    // let n_modules: usize = df["n_modules"].get(0).unwrap().try_extract().unwrap();
    let n_modules: usize = modules.len();
    let mut gt_result = Vec::new();
    for (i, module) in modules.iter().enumerate().take(n_modules-1).skip(1) {
        let data = extract_from_df(df, i).unwrap();
        let (counts, lengths, is_spanning, max_spanning_reps, max_overall_reps) = data;

        let model = Model::new(&lengths, max_spanning_reps as usize, max_overall_reps as usize);
        let likelihoods = model.evaluate(&counts, &lengths, &is_spanning);
        let predictions_enum = model.predict_enum(likelihoods.clone());
        let confidences = model.get_conf(likelihoods.clone());

        let module_df = get_module_df(df, i).unwrap();
        let predictions_seq = get_predictions_seqs(&module_df, module, predictions_enum);

        let result = ModuleResult{ predictions_enum, predictions_seq, confidences, likelihoods };
        gt_result.push(result);
    }
    let gt_result = GenotypingResults{modules: gt_result};
    return gt_result;
}

// pub(crate) fn print_json_file(gt_results: &GenotypingResults, p: &Path) -> Result<(), Box<dyn Error>> {
//     let mut file = File::create(p)?;
//     let json = serde_json::to_string(&gt_results)?;
//     write!(file, "{}", json)?;
//     return Ok(());
// }

fn get_predictions_seqs(module_df: &DataFrame, module: &Module, prediction: (Prediction, Prediction)) -> (String, String) {
    // select spanning reads
    let f = |o: Option<&str>| { let x = o.unwrap(); x == "Spanning" };
    let mask: BooleanChunked = module_df.column("classes").unwrap().str().unwrap().iter().map(f).collect();
    let spanning_df = module_df.filter(&mask).unwrap();

    let nomenclatures_df = get_nomenclature_counts(&spanning_df);

    use Prediction as P;
    let result = match prediction {
        (P::Num(a), P::Num(b)) => {
            get_string_predictions(&nomenclatures_df, a, b)
        },
        (P::Num(a), y) => {
            let nomenclature = get_single_string_prediction(&nomenclatures_df, a);
            let nomenclature2 = format_nonnumeric_prediction(module, y);
            (nomenclature, nomenclature2)
        },
        (x, P::Num(b)) => {
            let nomenclature1 = format_nonnumeric_prediction(module, x);
            let nomenclature = get_single_string_prediction(&nomenclatures_df, b);
            (nomenclature1, nomenclature)
        },
        (x, y) => {
            let nomenclature1 = format_nonnumeric_prediction(module, x);
            let nomenclature2 = format_nonnumeric_prediction(module, y);
            (nomenclature1, nomenclature2)
        }
    };

    // println!("{:?}", nomenclatures_df);
    // println!("{:?}", prediction);
    // println!("{:?}", result);
    return result;
}

fn format_nonnumeric_prediction(module: &Module, prediction: Prediction) -> String {
    let module_str = match module {
        Module::Sequence(x) => str::from_utf8(x).unwrap(),
        Module::Repeat((x, _)) => str::from_utf8(x).unwrap(),
    };

    use Prediction as P;
    match prediction {
        P::Expansion => { return format!("{}[{}]", module_str, "E") },
        P::Background => { return format!("{}[{}]", module_str, "B") },
        P::Num(_) => { panic!("Oopsie.") /* This should be unreachable, but I do not have the guts to make it unreachable, because it relies on the caller. */ }
    }
}

fn get_single_string_prediction(df: &DataFrame, x: usize) -> String {
    let nomenclature = df
        .filter(&df.column("counts").unwrap().u64().unwrap().equal(x)).unwrap()
        .column("nomenclatures").unwrap()
        .str().unwrap()
        .get(0).unwrap()
        .to_string();
    return nomenclature;
}

fn get_string_predictions(df: &DataFrame, a: usize, b: usize) -> (String, String) {
    let cols_in = ["n_occ", "counts", "nomenclatures"];
    debug_assert!(df.get_column_names() == cols_in);

    if a == b {
        // similar to ./../../dante_py/dante_remastr_simple.py:361:1  (use gF)
        let relevant_nomenclatures = df.filter(&df.column("counts").unwrap().u64().unwrap().equal(a)).unwrap();

        let noms = relevant_nomenclatures.column("nomenclatures").unwrap().str().unwrap().head(Some(2));
        let noms: Vec<&str> = noms.iter().flatten().collect();

        let occs = relevant_nomenclatures.column("n_occ").unwrap().u32().unwrap().head(Some(2));
        let occs: Vec<u32> = occs.iter().flatten().collect();

        match noms.len() {
            2 => {
                const ASSIGNMENT_FACTOR: u32 = 5;
                if occs[0] >= occs[1] * ASSIGNMENT_FACTOR {
                    // [5, 1] -> 3 to a1, 2 to a2, 1 to err
                    // [4, 1] -> 4 to a1, 1 to a2
                    return (noms[0].to_string(), noms[0].to_string());
                } else {
                    return (noms[0].to_string(), noms[1].to_string());
                }
            },
            1 => {
                return (noms[0].to_string(), noms[0].to_string());
            },
            0 => { panic!("While theoretically possible, practically you should never get here."); }
            _ => { panic!("Unexpected number of nomenclatures."); }
        }
    } else {
        let allele1_nomenclature = get_single_string_prediction(df, a);
        let allele2_nomenclature = get_single_string_prediction(df, b);
        return (allele1_nomenclature, allele2_nomenclature);
    }
}

fn get_nomenclature_counts(df: &DataFrame) -> DataFrame {
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
fn extract_from_df(df: &DataFrame, idx: usize) -> Result<(Vec<u64>, Vec<u64>, Vec<bool>, u64, u64), Box<dyn Error>> {
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
fn get_module_df(df: &DataFrame, idx: usize) -> Result<DataFrame, Box<dyn Error>> {
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

#[test]
#[allow(non_snake_case)]
fn test_genotyping_ALS_motif() {
    use serde_json::{from_value, from_str, Value};
    let data = include_str!("./../genotyping_test_data.txt");
    let lines: Vec<_> = data.split("\n").collect();
    let tmp1 = lines[0].strip_prefix("JSON: ").unwrap();
    let tmp2: Value = from_str(tmp1).unwrap();
    let obj = tmp2.as_object().unwrap();

    let fl_counts:  Vec<u64> = from_value(obj["flanking_observed_counts"].clone()).unwrap();
    let sp_counts:  Vec<u64> = from_value(obj["spanning_observed_counts"].clone()).unwrap();
    let fl_lengths: Vec<u64> = from_value(obj["flanking_read_lengths"].clone()).unwrap();
    let sp_lengths: Vec<u64> = from_value(obj["spanning_read_lengths"].clone()).unwrap();
    let prediction = obj["prediction"].as_array().unwrap();

    let lengths: Vec<u64> = [sp_lengths, fl_lengths].concat();
    let is_spanning: Vec<bool> = [vec![true; sp_counts.len()], vec![false; fl_counts.len()]].concat();
    let max_spanning_reps: u64 = *sp_counts.iter().max().unwrap();
    let counts: Vec<u64> = [sp_counts, fl_counts].concat();
    let max_overall_reps: u64 = *counts.iter().max().unwrap();

    // println!("{:?}", counts);
    // println!("{:?}", lengths);
    // println!("{:?}", is_spanning);
    // println!("{:?}", max_spanning_reps);
    // println!("{:?}", max_overall_reps);
    // println!("{:?}", is_monoa);
    println!("{:?}", prediction);


    let model = Model::new(&lengths, max_spanning_reps as usize, max_overall_reps as usize);
    println!("{:?}", model);
    let likelihoods = model.evaluate(&counts, &lengths, &is_spanning);
    // let pred_sym = model.predict_sym(likelihoods.clone());
    let confidences = model.get_conf(likelihoods.clone());
    println!("{:?}", likelihoods);
    // println!("{:?}", pred_sym);
    println!("{:?}", confidences);
}

#[test]
fn test_genotyping_all_motifs() {
    use serde_json::{from_value, from_str, Value};
    let data = include_str!("./../genotyping_test_data.txt");
    let lines: Vec<_> = data.trim().split("\n").collect();
    for line in lines {
        let tmp1 = line.strip_prefix("JSON: ").unwrap();
        let tmp2: Value = from_str(tmp1).unwrap();
        let obj = tmp2.as_object().unwrap();

        let fl_counts:  Vec<u64> = from_value(obj["flanking_observed_counts"].clone()).unwrap();
        let sp_counts:  Vec<u64> = from_value(obj["spanning_observed_counts"].clone()).unwrap();
        let fl_lengths: Vec<u64> = from_value(obj["flanking_read_lengths"].clone()).unwrap();
        let sp_lengths: Vec<u64> = from_value(obj["spanning_read_lengths"].clone()).unwrap();
        let prediction = obj["prediction"].as_array().unwrap();

        let lengths = [sp_lengths, fl_lengths].concat();
        let is_spanning = [vec![true; sp_counts.len()], vec![false; fl_counts.len()]].concat();
        let max_spanning_reps = *sp_counts.iter().max().unwrap();
        let counts = [sp_counts, fl_counts].concat();
        let max_overall_reps = *counts.iter().max().unwrap();

        let model = Model::new(&lengths, max_spanning_reps as usize, max_overall_reps as usize);
        let likelihoods = model.evaluate(&counts, &lengths, &is_spanning);
        // let pred_sym = model.predict_sym(likelihoods.clone());
        let confidences = model.get_conf(likelihoods.clone());

        println!("{:?}", prediction);
        // println!("{:?}", pred_sym);
        println!("{:?}", confidences);
        println!();
    }
}

#[derive(Debug)]
struct Model {
     max_rep: usize,
    max_frep: usize,
     exp_idx: usize,
     bkg_idx: usize,

       rdist: Vec<f64>,
      mprobs: Vec<f64>,
      models: Vec<Vec<f64>>
}

impl Model {
    const L_OTHERS: f64 = 1.00;
    const L_EXP:    f64 = 1.01;
    const L_BKG:    f64 = 0.01;

    const P_DEL1:   f64 = 0.0001;
    const P_DEL2:   f64 = 0.0001;
    const P_INS:    f64 = 0.0001;

    fn new(lengths: &[u64], max_spanning_rep: usize, max_flanking_rep: usize) -> Self {
        Self { 
            max_rep: max_spanning_rep,
            max_frep: max_flanking_rep,
            exp_idx: max_spanning_rep + 1,
            bkg_idx: max_spanning_rep + 2,
            rdist: Self::construct_rdist(lengths),
            mprobs: Self::construct_mprobs(max_spanning_rep),
            models: Self::construct_models(max_spanning_rep, max_flanking_rep)
        }
    }

    fn construct_rdist(lengths: &[u64]) -> Vec<f64> {
        let n = *lengths.iter().max().unwrap();
        let mut result = vec![0.0; (n + 1) as usize];
        for item in lengths {
            let e = *item as usize;
            result[e] += 1.0;
        }
        let m = lengths.len() as u64;
        for item in &mut result { *item /= m as f64; }
        return result;
    }

    fn construct_mprobs(max_spanning_rep: usize) -> Vec<f64> {
        let mut mprobs = Vec::with_capacity(max_spanning_rep + 3);
        for _ in 0..=max_spanning_rep { mprobs.push(Self::L_OTHERS); }
        mprobs.push(Self::L_EXP);
        mprobs.push(Self::L_BKG);
        return mprobs;
    }

    fn construct_models(max_spanning_rep: usize, max_flanking_rep: usize) -> Vec<Vec<f64>> {
        let mut models = Vec::with_capacity(max_spanning_rep + 3);
        for i in 0..=max_spanning_rep { models.push(Self::model_full(max_flanking_rep, i)); }
        models.push(Self::model_expn(max_flanking_rep, max_spanning_rep));
        models.push(Self::model_bckg(max_flanking_rep));
        return models;
    }

    fn model_full(max_flanking_rep: usize, gt: usize) -> Vec<f64> { 
        let p_del = (Model::P_DEL1 + Model::P_DEL2 * (gt as f64)).clamp(0.0, 1.0);
        let deletes: Vec<f64> = {
            let gt: u64 = gt.try_into().unwrap();
            let dist = Binomial::new(p_del, gt).unwrap();
            (0..=gt).map(|i| dist.pmf(i)).collect()
        };

        let p_ins = Model::P_INS;
        let inserts: Vec<f64> = {
            let gt: u64 = gt.try_into().unwrap();
            let dist = Binomial::new(p_ins, gt).unwrap();
            (0..=gt).map(|i| dist.pmf(i)).collect()
        };

        let mut result = vec![0.0; max_flanking_rep + 1];
        let mut deletes = deletes;
        deletes.reverse();
        let r = convolve(&inserts, &deletes);
        let x = usize::min(result.len(), r.len());
        result[..x].copy_from_slice(&r[..x]);
        return result
    }

    fn model_expn(max_flanking_rep: usize, max_spanning_rep: usize) -> Vec<f64> {
        let mut result = vec![0.0; max_flanking_rep + 1];

        for i in (max_spanning_rep + 1)..=(max_flanking_rep + 1) {
            let tmp = Model::model_full(max_flanking_rep, i);
            for j in 0..tmp.len() { result[j] += tmp[j]; }
        }
        let n: f64 = (max_flanking_rep - max_spanning_rep + 1) as f64;
        for item in &mut result { *item /= n; }
        return result;
    }

    fn model_bckg(max_flanking_rep: usize) -> Vec<f64> {
        let x = 1.0 / ((max_flanking_rep + 1) as f64);
        let result = vec![x; max_flanking_rep + 1];
        return result;
    }
}

fn convolve(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0; a.len() + b.len() - 1];
    for j in 0..b.len() {
        for i in 0..a.len() {
            result[j + i] += b[j] * a[i];
        }
    }
    return result;
}

impl Model {
    fn evaluate(&self, observed: &[u64], rlengths: &[u64], spanning: &[bool]) -> ndarray::Array2<f64> {
        let n = self.max_rep + 3;   // 0, 1, ..., n, E, B
        let mut result = Array::from_elem((n, n), f64::NEG_INFINITY);

        for g1_idx in 0..n {
            for g2_idx in g1_idx..n {
                result[[g1_idx, g2_idx]] = self.loglikelihood_of_D_given_G(observed, rlengths, spanning, g1_idx, g2_idx);
            }
        }
        return result;
    }

    // def loglikelihood_of_D_given_G(self, obs_counts: list[int], read_lengths: list[int], is_spanning: list[bool], g1_idx: int, g2_idx: int) -> float:
    #[allow(non_snake_case)]
    fn loglikelihood_of_D_given_G(
        &self, observed: &[u64], rlengths: &[u64], spanning: &[bool], g1_idx: usize, g2_idx: usize
    ) -> f64 {
        let mut m_lh = 0.0;
        for (&oc, &rl, &sf) in izip!(observed, rlengths, spanning) {
            let bckgrnd_l = self.l_read_given_genotype(oc, rl, sf, self.bkg_idx);
            let allele1_l = self.l_read_given_genotype(oc, rl, sf, g1_idx);
            let allele2_l = self.l_read_given_genotype(oc, rl, sf, g2_idx);
            m_lh += (bckgrnd_l + allele1_l + allele2_l).ln();
        }
        return m_lh;
    }

    #[allow(clippy::needless_late_init)]
    fn l_read_given_genotype(&self, oc: u64, rl: u64, is_spanning: bool, gt_idx: usize) -> f64 {
        let lh_cover: f64 = 1.0 / (rl as f64);          // P(b_i | a_i, r_i) # incorrect, but does something
        let lh_r_len: f64 = self.rdist[rl as usize];    // P(r_i)            # correct, but does nothing
        let lh_model: f64;                              // P(a_i | g_i)
        let lh_mprob: f64 = self.mprobs[gt_idx];        // P(g_i)

        lh_model = if is_spanning {
            self.models[gt_idx][oc as usize]
        } else {
            let tmp = &self.models[gt_idx][(oc as usize)..];
            let numerator: f64 = tmp.iter().sum();
            let denominator: f64 = tmp.len() as f64;
            numerator / denominator
        };

        return lh_cover * lh_r_len * lh_model * lh_mprob;
    }

    fn predict(llmatrix: ndarray::Array2<f64>) -> (usize, usize) {
        let result = llmatrix
            .indexed_iter()
            .filter(|(_, x)| !x.is_nan())
            .max_by(|(_, x), (_, y)| x.partial_cmp(y).unwrap())
            .map(|((r, c), _)| (r, c))
            .unwrap();

        return result;
    }

    fn predict_enum(&self, llmatrix: ndarray::Array2<f64>) -> (Prediction, Prediction) {
        let (row, col) = Self::predict(llmatrix);
        let to_enum = |x: usize| {
            match x {
                x if x == self.bkg_idx => { Prediction::Background },
                x if x == self.exp_idx => { Prediction::Expansion },
                x => { Prediction::Num(x) } 
            }
        };
        return (to_enum(row), to_enum(col));
    }

    fn get_conf(&self, llmatrix: ndarray::Array2<f64>) -> [f64; 7] {
        // llmatrix = llmatrix - np.max(llmatrix)
        let llmatrix2 = (&llmatrix) - (&llmatrix).max();

        // softmax
        let llmatrix_exp = llmatrix2.mapv_into(|x| x.exp());
        let prob = &llmatrix_exp / llmatrix_exp.sum();

        let pred = Self::predict(llmatrix);
        let conf_pred = prob[pred];
        let conf_al_1 = prob.slice(s![pred.0, ..]).sum() + prob.slice(s![.., pred.0]).sum() - prob[[pred.0, pred.0]];
        let conf_al_2 = prob.slice(s![pred.1, ..]).sum() + prob.slice(s![.., pred.1]).sum() - prob[[pred.1, pred.1]];

        let bkg = self.bkg_idx;
        let conf_bckg = prob[[bkg, bkg]];
        let conf_bg_t = prob.slice(s![bkg, ..]).sum() + prob.slice(s![.., bkg]).sum() - prob[[bkg, bkg]];

        let exp = self.exp_idx;
        let conf_expn = prob[[exp, exp]];
        let conf_ex_t = prob.slice(s![exp, ..]).sum() + prob.slice(s![.., exp]).sum() - prob[[exp, exp]];

        return [conf_pred, conf_al_1, conf_al_2, conf_bckg, conf_bg_t, conf_expn, conf_ex_t];
    }
}
