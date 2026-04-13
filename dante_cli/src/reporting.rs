use std::{error::Error, fs::File, io::Read};
use std::io::Write;

use minijinja::{context, Environment};
use serde_json::Value;
use serde::{Serialize, Deserialize};

pub(crate) fn report(args: &crate::ArgsNew) -> Result<(), Box<dyn Error>> {
    let motif_names = [  // get these from args.motif_file
        "HMNR7", "SCA37", "NIID_ETM6", "FAME2", "SPD", "GAD", "SCA7",
        "DM2", "BPES", "FAME4", "HD", "CANVAS", "CCHS", "FAME7", "FAME3", 
        "SCA12", "SCA1", "CCD", "SCA17", "HGF", "CF", "OPDM1", "FAME1",
        "FTD_ALS", "FRDA", "HSAN8", "OPML1", "JBS", "DRPLA", "SCA2", "OPDM4",
        "SCA8", "HPE5", "SCA27B", "OPMD", "SCA3", "ALS", "BSS", "FAME6", "SCA",
        "SCA4", "HDL2", "RCPS", "FECD3", "SCA6", "OPDM2", "PSACH_MED", "DM1",
        "SCA36", "CJD", "EPM1", "TOF", "SCA10", "DEE1_MRXARX_PRTS", "DEE1_MRXARX",
        "DMD", "SBMA", "VACTERLX", "PHPX_XLMR", "FRAXA_FXTAS_FXPOI", "FRAXE", "FRAXF"
    ];

    let mut motifs = Vec::new();
    let mut buf = String::new();
    for motif in motif_names {
        let filename = args.output.join(motif.to_owned() + ".motif.json");
        let mut fp = File::open(filename).unwrap();
        let _ = fp.read_to_string(&mut buf).unwrap();
        let json: Value = serde_json::from_str(&buf).expect("JSON was not well-formatted");
        motifs.push(json);
        buf.clear();
    };

    let data = context! {
        dante_version => "0.14.0",  // TODO get_version from toml
        dante_params => context! {
            file_bam => args.bam_file,
            file_motif => args.motif_file,
            // is_male => true,        // TODO
            // max_noms => 5,          // TODO
        },
        motifs => motifs
    };

    let mut env = Environment::new();
    env.add_template("report_static", include_str!("../templates/report_static.html")).unwrap();
    let tmpl = env.get_template("report_static").unwrap();
    let result = tmpl.render(context! { data => data }).unwrap();
    let output_file = "report_static.html";  // TODO
    let mut file = File::create(output_file).unwrap();
    file.write_all(result.as_bytes()).unwrap();

    Ok(())
}

// pub(crate) fn report(output: &Path) -> Result<(), Box<dyn Error>> {
//     let annotations: Vec<String> = vec!["./output/motifs/ALS.annotations.tsv".to_string()];
//     let genotypes: Vec<String>   = vec!["./output/motifs/ALS.genotypes.json".to_string()];
//     // read tsv
//     let file = File::open(&annotations[0]).unwrap();
//     let opts = CsvReadOptions::default().with_parse_options(CsvParseOptions::default().with_separator(b'\t'));
//     let df = CsvReader::new(file).with_options(opts).finish().unwrap();
// 
//     // read json
//     let mut x = File::open(&genotypes[0]).unwrap();
//     let mut buf = String::new();
//     let y = x.read_to_string(&mut buf).unwrap();
//     let json: Value = serde_json::from_str(&buf).expect("JSON was not well-formatted");
// 
//     // collect data
//     let ctx = context!(name => "John");
// 
//     // create alignment reports
// 
//     // create main report
//     let mut env = Environment::new();
//     // let tmp = include_str!("./../../dante_py/templates/report_template2_static.html");
//     // let tmp = include_str!("./report_static.html");
//     env.add_template("hello", "Hello {{ name }}!").unwrap();
//     let tmpl = env.get_template("hello").unwrap();
//     let result = tmpl.render(ctx).unwrap();
// 
//     // ...profit
//     println!("{:?}", df);
//     println!("{:?}", json);
//     println!("{:?}", y);
//     println!("{:?}", annotations);
//     println!("{:?}", genotypes);
//     println!("{}", result);
//     return Ok(());
// }

#[test]
fn test_reporting() {
    use std::io::Write;
    // cargo test report
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/ALS.annotations.tsv
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/ALS.genotypes.json
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/DM2.annotations.tsv
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/DM2.genotypes.json
    //
    // ../../../../analyses/2025-12-15_gen_validation_v3/dante_out/report/report_static.html

    let data = context! {
        dante_version => "0.14.0",
        dante_params => context! {
            file_bam => "../../../../analyses/2026-03-04_validation_v4/data/in_HG002.GRCh38.selected_w_pairs.bam",
            file_motif => "../../../../analyses/2026-03-04_validation_v4/data/01_in_dante_nomenclatures_predominant.tsv",
            is_male => true,
            max_noms => 5,
        },
        motifs => generate_motifs()
    };

    let mut env = Environment::new();
    env.add_template("report_static", include_str!("../templates/report_static.html")).unwrap();
    let tmpl = env.get_template("report_static").unwrap();
    let result = tmpl.render(context! { data => data }).unwrap();
    let mut file = File::create("report_static.html").unwrap();
    file.write_all(result.as_bytes()).unwrap();
}

#[cfg(test)]
fn generate_motifs() -> Vec<minijinja::Value> {
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/ALS.annotations.tsv
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/ALS.genotypes.json
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/DM2.annotations.tsv
    // ../../../../analyses/2026-03-04_validation_v4/output/motifs/DM2.genotypes.json
    // assert_eq!(generate_ALS(), generate_ALS2());
    let motifs = vec![
        generate_ALS(),
        generate_ALS2(),
    ];
    return motifs;
}

#[allow(non_snake_case)]
#[cfg(test)]
fn generate_ALS() -> minijinja::Value {
    const NINF: f64 = f64::NEG_INFINITY;
    context! {
        // ALS chr15:g.22786680_22786703GGC[8] GGGCGGAATGGGGACTGCAGCTGCGGCAGC CGGGGAGGGGGCGCGTAGCCCGAGCCCCGC
        motif_id => "ALS",
        sequence => vec![
            "GGGCGGAATGGGGACTGCAGCTGCGGCAGC",
            "GGC[8]",
            "CGGGGAGGGGGCGCGTAGCCCGAGCCCCGC"
        ],
        nomenclatures => vec![
            (156, vec!["GGC[8]"]),
            (6, vec!["GGC[7]GGG[1]"]),
            (1, vec!["GGC[1]GGG[1]GGC[3]GGG[3]"]),
            (1, vec!["GGC[2]GGG[1]GGC[4]GGG[1]"]),
            (1, vec!["GGC[3]GGG[1]GGC[1]GGG[2]GGC[1]"]),
            // and 13 more
        ],
        modules => vec![
            context! {
                module_id => 0,
                allele_1 => context! {
                    num_pred => 8.to_string(),
                    num_conf => 1.0,
                    num_reads_spanning => 177,
                    seq_pred => "GGC[8]",
                    seq_dist => 0,
                    seq_reads_spanning => 156,
                },
                allele_2 => context! {
                    num_pred => 8.to_string(),
                    num_conf => 1.0,
                    num_reads_spanning => 177,
                    seq_pred => "GGC[8]",
                    seq_dist => 0,
                    seq_reads_spanning => 156,
                },
                overall => context! {
                    conf => 1.0,
                    reads_spanning_num_nonspec => 1,
                    reads_spanning_seq_nonspec => 22,
                    reads_flanking => 68,
                    reads_inrepeat => 0,
                    reads_total => 256, /* 177(Spanning) + 1(Spanning) + 68(Flanking) + 0(Inrepeat) + 10(Missing)*/
                },
                nomenclatures => vec![
                    (156, vec!["GGC[8]"]),
                    (6, vec!["GGC[7]GGG[1]"]),
                    (1, vec!["GGC[1]GGG[1]GGC[3]GGG[3]"]),
                    (1, vec!["GGC[2]GGG[1]GGC[4]GGG[1]"]),
                    (1, vec!["GGC[3]GGG[1]GGC[1]GGG[2]GGC[1]"]),
                    // and 13 more
                ],
                histogram_data => context! {
                    spanning => vec![0, 0, 0, 0, 0 , 1 , 0, 0, 177, 0, 0],
                    flanking => vec![0, 0, 6, 7, 11, 12, 9, 6, 16 , 1, 0],
                    inrepeat => vec![0, 0, 0, 0, 0 , 0 , 0, 0, 0  , 0, 0],
                },
                heatmap_data => context! {
                    matrix => vec![
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
                    xlim => vec![3, 11],
                    ylim => vec![3, 11]
                }
            }
        ]
    }
}

#[allow(non_snake_case)]
#[cfg(test)]
fn generate_ALS2() -> minijinja::Value {
    const NINF: f64 = f64::NEG_INFINITY;
    let result = MotifData {
        motif_id: "ALS2".to_string(),
        sequence: vec![
            "GGGCGGAATGGGGACTGCAGCTGCGGCAGC".to_string(),
            "GGC[8]".to_string(),
            "CGGGGAGGGGGCGCGTAGCCCGAGCCCCGC".to_string()
        ],
        nomenclatures: vec![
            (156, vec!["GGC[8]".to_string()]),
            (6, vec!["GGC[7]GGG[1]".to_string()]),
            (1, vec!["GGC[1]GGG[1]GGC[3]GGG[3]".to_string()]),
            (1, vec!["GGC[2]GGG[1]GGC[4]GGG[1]".to_string()]),
            (1, vec!["GGC[3]GGG[1]GGC[1]GGG[2]GGC[1]".to_string()]),
            // and 13 more
        ],
        modules: vec![
            ModuleData {
                module_id: 0,
                allele_1: AlleleData {
                    num_pred: 8.to_string(),
                    num_conf: 1.0,
                    num_reads_spanning: 177,
                    seq_pred: "GGC[8]".to_string(),
                    seq_dist: 0,
                    seq_reads_spanning: 156,
                },
                allele_2: AlleleData {
                    num_pred: 8.to_string(),
                    num_conf: 1.0,
                    num_reads_spanning: 177,
                    seq_pred: "GGC[8]".to_string(),
                    seq_dist: 0,
                    seq_reads_spanning: 156,
                },
                overall: OverallData {
                    conf: 1.0,
                    reads_spanning_num_nonspec: 1,
                    reads_spanning_seq_nonspec: 22,
                    reads_flanking: 68,
                    reads_inrepeat: 0,
                    reads_total: 256, /* 177(Spanning) + 1(Spanning) + 68(Flanking) + 0(Inrepeat) + 10(Missing)*/
                },
                nomenclatures: vec![
                    (156, vec!["GGC[8]".to_string()]),
                    (6, vec!["GGC[7]GGG[1]".to_string()]),
                    (1, vec!["GGC[1]GGG[1]GGC[3]GGG[3]".to_string()]),
                    (1, vec!["GGC[2]GGG[1]GGC[4]GGG[1]".to_string()]),
                    (1, vec!["GGC[3]GGG[1]GGC[1]GGG[2]GGC[1]".to_string()]),
                    // and 13 more
                ],
                histogram_data: HistogramData {
                    spanning: vec![0, 0, 0, 0, 0 , 1 , 0, 0, 177, 0, 0],
                    flanking: vec![0, 0, 6, 7, 11, 12, 9, 6, 16 , 1, 0],
                    inrepeat: vec![0, 0, 0, 0, 0 , 0 , 0, 0, 0  , 0, 0],
                },
                heatmap_data: HeatmapData {
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
        ]
    };

    let x: serde_json::Value = serde_json::to_value(&result).unwrap();
    return minijinja::Value::from_serialize(x);
}

#[derive(Serialize, Deserialize)]
struct MotifData {
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
    num_pred: String,
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
