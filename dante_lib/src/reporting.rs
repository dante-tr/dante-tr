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
    seq_dist: u64,
    seq_reads_spanning: u64
}

#[derive(Serialize, Deserialize)]
struct OverallData {
    conf: f64,
    reads_spanning_num_nonspec: u64,
    reads_spanning_seq_nonspec: u64,
    reads_flanking: u64,
    reads_inrepeat: u64,
    reads_total: u64
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
    xlim: Vec<u64>,
    ylim: Vec<u64>
}

// -----------------------------------------------------------------------------
use polars::prelude::DataFrame;
use crate::genotyping::{GenotypingResults, ModuleResult};
use crate::io::TRRecord;

impl MotifData {
    pub(crate) fn create(tr: &TRRecord, df: &DataFrame, gt: &GenotypingResults) -> MotifData {
        use crate::df_ops::get_nomenclatures;
        let n = tr.copy_unit.len();

        let motif_id = tr.name.clone();
        let sequence = tr.to_sequence();
        let nomenclatures = get_nomenclatures(df, 0..n); // TODO

        let mut modules = Vec::with_capacity(n);
        for i in 0..n { modules.push(
            ModuleData::create(i, &gt.modules[i], df)
        ); }

        let result = MotifData { motif_id, sequence, nomenclatures, modules };
        return result;
    }
}

impl ModuleData {
    fn create(idx: usize, mr: &ModuleResult, df: &DataFrame) -> ModuleData {
        use crate::df_ops::get_nomenclatures;
        ModuleData {
            module_id: idx,
            allele_1: AlleleData::create(0, mr),
            allele_2: AlleleData::create(1, mr),
            overall: OverallData::create(),
            nomenclatures: get_nomenclatures(df, idx..(idx+1)), // TODO
            histogram_data: HistogramData::create(),
            heatmap_data: HeatmapData::create()
        }
    }
}

impl AlleleData {
    fn create(allele: usize, mr: &ModuleResult) -> AlleleData {
        match allele {
            0 => {
                let num_pred = mr.predictions_enum.0.to_string();
                let seq_pred = mr.predictions_seq.0.clone();
                AlleleData {
                    num_pred,
                    num_conf: 1.0,
                    num_reads_spanning: 177,
                    seq_pred,
                    seq_dist: 0,
                    seq_reads_spanning: 156,
                }
            },
            1 => {
                let num_pred = mr.predictions_enum.1.to_string();
                let seq_pred = mr.predictions_seq.1.clone();
                AlleleData {
                    num_pred,
                    num_conf: 1.0,
                    num_reads_spanning: 177,
                    seq_pred,
                    seq_dist: 0,
                    seq_reads_spanning: 156,
                }
            },
            _ => { unreachable!() }
        }
    }
}

impl OverallData {
    fn create() -> OverallData {
        OverallData {
            conf: 1.0,
            reads_spanning_num_nonspec: 1,
            reads_spanning_seq_nonspec: 22,
            reads_flanking: 68,
            reads_inrepeat: 0,
            reads_total: 256, /* 177(Spanning) + 1(Spanning) + 68(Flanking) + 0(Inrepeat) + 10(Missing)*/
        }
    }
}

impl HistogramData {
    fn create() -> HistogramData {
        HistogramData {
            spanning: vec![0, 0, 0, 0, 0 , 1 , 0, 0, 177, 0, 0],
            flanking: vec![0, 0, 6, 7, 11, 12, 9, 6, 16 , 1, 0],
            inrepeat: vec![0, 0, 0, 0, 0 , 0 , 0, 0, 0  , 0, 0],
        }
    }
}

impl HeatmapData {
    fn create() -> HeatmapData {
        const NINF: f64 = f64::NEG_INFINITY;
        HeatmapData {
            matrix: vec![
                vec![-2928.62, -2928.54, -2899.41, -2864.29, -2807.24, -2736.91, -2692.38, -2562.01, -1338.86, -2296.59, -2758.10],
                vec![NINF    , -2928.47, -2899.41, -2864.29, -2807.24, -2736.91, -2692.38, -2562.01, -1338.86, -2296.58, -2758.07],
                vec![NINF    , NINF    , -2895.08, -2860.16, -2803.11, -2732.78, -2688.25, -2557.87, -1334.73, -2290.14, -2733.10],
                vec![NINF    , NINF    , NINF    , -2854.82, -2798.28, -2727.95, -2683.42, -2553.05, -1329.90, -2282.60, -2702.96],
                vec![NINF    , NINF    , NINF    , NINF    , -2789.58, -2720.36, -2675.76, -2545.13, -1321.98, -2270.42, -2653.81],
                vec![NINF    , NINF    , NINF    , NINF    , NINF    , -2710.43, -2662.29, -2530.29, -1307.13, -2250.91, -2592.33],
                vec![NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , -2659.78, -2529.31, -1306.17, -2246.45, -2553.85],
                vec![NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , -2462.26, -1303.52, -2220.61, -2469.42],
                vec![NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , -1169.61, -1304.86, -1337.24],
                vec![NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , -2143.83, -2265.88],
                vec![NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , NINF    , -2658.36]
            ],
            xlim: vec![3, 11],
            ylim: vec![3, 11]
        }
    }
}


