use polars::prelude::*;
use serde::{Serialize, Deserialize};

use crate::genotyping::Prediction;

#[derive(Debug, Serialize, Deserialize, Default)]
pub(crate) struct PhasingResults {
    seq1: Vec<String>,
    seq2: Vec<String>,
    num1: Vec<Prediction>,
    num2: Vec<Prediction>
}

impl PhasingResults {
    fn from_gt_results(genotypes: &GenotypingResults) -> PhasingResults {
        let mut seq1 = Vec::new();
        let mut seq2 = Vec::new();
        let mut num1 = Vec::new();
        let mut num2 = Vec::new();
        for module in &genotypes.modules {
            let (s1, s2) = &module.predictions_seq;
            seq1.push(s1.clone());
            seq2.push(s2.clone());
            let (n1, n2) = module.predictions_enum;
            num1.push(n1);
            num2.push(n2);
        }
        return PhasingResults { seq1, seq2, num1, num2 };
    }

    fn swap_at(&mut self, idx: usize) {
        std::mem::swap(&mut self.num1[idx], &mut self.num2[idx]);
        std::mem::swap(&mut self.seq1[idx], &mut self.seq2[idx]);
    }

    fn is_homo_at(&self, idx: usize) -> bool {
        return self.seq1[idx] == self.seq2[idx];
    }

    fn is_crossing(&self, idx1: usize, idx2: usize, count: impl Fn(usize, String, usize, String) -> usize) -> bool {
        let n_seq1_idx1_seq1_idx2 = count(idx1, self.seq1[idx1].clone(), idx2, self.seq1[idx2].clone());
        let n_seq2_idx1_seq2_idx2 = count(idx1, self.seq2[idx1].clone(), idx2, self.seq2[idx2].clone());
        let n_seq1_idx1_seq2_idx2 = count(idx1, self.seq1[idx1].clone(), idx2, self.seq2[idx2].clone());
        let n_seq2_idx1_seq1_idx2 = count(idx1, self.seq2[idx1].clone(), idx2, self.seq1[idx2].clone());
        let crossing_score = n_seq1_idx1_seq2_idx2 + n_seq2_idx1_seq1_idx2;
        let straight_score = n_seq1_idx1_seq1_idx2 + n_seq2_idx1_seq2_idx2;

        if crossing_score > straight_score {
            return true;
        } else if crossing_score < straight_score {
            return false;
        } else {
            println!("Cannot decide. Putting noncrossing.");
            return false;
        }
    }
}


fn get_n_co_occurrences(motif_df: &DataFrame, idx1: usize, seq1: String, idx2: usize, seq2: String) -> usize {
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

use crate::genotyping::GenotypingResults;
pub(crate) fn phase(motif_df: &DataFrame, genotypes: &GenotypingResults) -> PhasingResults {

    let mut phasing_results = PhasingResults::from_gt_results(genotypes);
    let n = genotypes.modules.len();
    if n == 1 { // There is nothing to phase
        return phasing_results;
    }

    let get_co_occurrences = |a, b, c, d| get_n_co_occurrences(&motif_df.clone(), a, b, c, d);
    let mut i = 0;
    let mut j;
    loop {
        while i < n && phasing_results.is_homo_at(i) { i += 1; }
        j = i + 1;
        while j < n && phasing_results.is_homo_at(j) { j += 1; }
        if j >= n { break }

        if phasing_results.is_crossing(i, j, get_co_occurrences) {
            phasing_results.swap_at(j);
        }
        i = j;
    }

    return phasing_results;
}

pub(crate) fn phase2(motif_df: &DataFrame, genotypes: &GenotypingResults) -> GenotypingResults {
    let mut result = genotypes.clone();
    let n = genotypes.modules.len();
    if n == 1 { // There is nothing to phase
        return result;
    }

    let get_co_occurrences = |a, b, c, d| get_n_co_occurrences(&motif_df.clone(), a, b, c, d);
    let mut i = 0;
    let mut j;
    loop {
        while i < n && result.is_homo_at(i) { i += 1; }
        j = i + 1;
        while j < n && result.is_homo_at(j) { j += 1; }
        if j >= n { break }

        if result.is_crossing(i, j, get_co_occurrences) {
            result.swap_at(j);
        }
        i = j;
    }

    return result;

}
