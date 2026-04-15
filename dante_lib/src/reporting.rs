use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize)]
pub(crate) struct MotifData {
    motif_id: String,
    sequence: Vec<String>,
    nomenclatures: Vec<(u64, Vec<String>)>,
    modules: Vec<ModuleData>
}

#[derive(Serialize, Deserialize)]
struct ModuleData {
    module_id: usize,
    allele_1: AlleleData,
    allele_2: AlleleData,
    overall: OverallData,
    nomenclatures: Vec<(u64, Vec<String>)>,
    histogram_data: HistogramData,
    heatmap_data: HeatmapData
}

#[derive(Serialize, Deserialize)]
struct AlleleData {
    num_pred: String, // Can be usize or E or B
    num_conf: f64,
    num_reads_spanning: u64,
    seq_pred: String,
    seq_dist: Option<u64>,
    seq_reads_spanning: u64
}

#[derive(Serialize, Deserialize)]
struct OverallData {
    conf: f64,
    reads_spanning_num_nonspec: u64,
    reads_spanning_seq_nonspec: u64,
    reads_flanking: u64,
    reads_inrepeat: u64,
    reads_total: usize
}

#[derive(Serialize, Deserialize)]
struct HistogramData {
    spanning: Vec<u64>,
    flanking: Vec<u64>,
    inrepeat: Vec<u64>
}

#[derive(Serialize, Deserialize)]
struct HeatmapData {
    matrix: Vec<Vec<f64>>,
    xlim: Vec<usize>,
    ylim: Vec<usize>
}

// -----------------------------------------------------------------------------
use polars::prelude::DataFrame;
use crate::genotyping::{GenotypingResults, ModuleResult};
use crate::io::TRRecord;

impl MotifData {
    pub(crate) fn create(tr: &TRRecord, df: &DataFrame, gt: &GenotypingResults) -> MotifData {
        let n = tr.copy_unit.len();

        let motif_id = tr.name.clone();
        let sequence = tr.to_sequence();
        let nomenclatures = crate::df_ops::get_nomenclatures(df, 0..n);

        let mut modules = Vec::with_capacity(n);
        for i in 0..n { modules.push(
            ModuleData::create(i, &gt.modules[i], df, tr)
        ); }

        let result = MotifData { motif_id, sequence, nomenclatures, modules };
        return result;
    }
}

impl ModuleData {
    fn create(idx: usize, mr: &ModuleResult, df: &DataFrame, tr: &TRRecord) -> ModuleData {
        ModuleData {
            module_id: idx,
            allele_1: AlleleData::create(0, mr, df, idx, tr),
            allele_2: AlleleData::create(1, mr, df, idx, tr),
            overall: OverallData::create(df, idx, mr),
            nomenclatures: crate::df_ops::get_nomenclatures(df, idx..(idx+1)),
            histogram_data: HistogramData::create(df, idx),
            heatmap_data: HeatmapData::create(mr)
        }
    }
}

impl AlleleData {
    fn create(allele: usize, mr: &ModuleResult, df: &DataFrame, idx: usize, tr: &TRRecord) -> AlleleData {
        use crate::df_ops::get_num_reads_spanning;
        use crate::df_ops::get_seq_reads_spanning;
        match allele {
            0 => {
                let num_pred_raw = mr.predictions_enum.0;
                let num_pred = num_pred_raw.to_string();
                let seq_pred = mr.predictions_seq.0.clone();
                let num_conf = mr.confidences[1];
                let num_reads_spanning = get_num_reads_spanning(df, idx, num_pred_raw);
                let seq_reads_spanning = get_seq_reads_spanning(df, idx, &seq_pred);
                let seq_dist = compute_distance(tr, num_pred_raw, &seq_pred);               // TODO
                AlleleData {
                    num_pred, num_conf, num_reads_spanning,
                    seq_pred, seq_dist, seq_reads_spanning,
                }
            },
            1 => {
                let num_pred_raw = mr.predictions_enum.1;
                let num_pred = num_pred_raw.to_string();
                let seq_pred = mr.predictions_seq.1.clone();
                let num_conf = mr.confidences[2];
                let num_reads_spanning = get_num_reads_spanning(df, idx, num_pred_raw);
                let seq_reads_spanning = get_seq_reads_spanning(df, idx, &seq_pred);
                let seq_dist = compute_distance(tr, num_pred_raw, &seq_pred);               // TODO
                AlleleData {
                    num_pred, num_conf, num_reads_spanning,
                    seq_pred, seq_dist, seq_reads_spanning,
                }
            },
            _ => { unreachable!() }
        }
    }
}

impl OverallData {
    fn create(df: &DataFrame, idx: usize, mr: &ModuleResult) -> OverallData {
        let conf = mr.confidences[0];
        let reads_spanning_num_nonspec = crate::df_ops::get_reads_spanning_num_nonspec(df, idx, mr.predictions_enum);
        let reads_spanning_seq_nonspec = crate::df_ops::get_reads_spanning_seq_nonspec(df, idx, &mr.predictions_seq);
        let reads_flanking = crate::df_ops::get_reads_flanking(df, idx);
        let reads_inrepeat = crate::df_ops::get_reads_inrepeat(df, idx);
        let reads_total = crate::df_ops::get_reads_total(df);
        OverallData {
            conf, reads_spanning_num_nonspec, reads_spanning_seq_nonspec,
            reads_flanking, reads_inrepeat,
            reads_total, /* 177(Spanning) + 1(Spanning) + 68(Flanking) + 0(Inrepeat) + 10(Missing)*/
        }
    }
}

impl HistogramData {
    fn create(df: &DataFrame, idx: usize) -> HistogramData {
        use crate::annotation::AClass;
        let mut spanning = crate::df_ops::get_histogram(df, idx, AClass::Spanning);
        let mut flanking = crate::df_ops::get_histogram(df, idx, AClass::Flanking);
        let mut inrepeat = crate::df_ops::get_histogram(df, idx, AClass::InRepeat);
        let new_len = *[spanning.len(), flanking.len(), inrepeat.len(), 11].iter().max().unwrap();
        spanning.resize(new_len, 0);
        flanking.resize(new_len, 0);
        inrepeat.resize(new_len, 0);
        HistogramData { spanning, flanking, inrepeat }
    }
}

impl HeatmapData {
    fn create(mr: &ModuleResult) -> HeatmapData {
        let (matrix, xlim, ylim) = mr.to_matrix();
        HeatmapData { matrix, xlim, ylim }
    }
}

// -----------------------------------------------------------------------------
fn compute_distance(_tr: &TRRecord, _pred: crate::genotyping::Prediction, _seq_pred: &str) -> Option<u64> {
    return None;
}
