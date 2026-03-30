use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;
use std::error::Error;

use crate::repeats::HGVSNomenclature;
use crate::hmm::Module;

pub fn get_modules(
    tr_record: &TRRecord
    // left_flank: &[u8], repeat: &TandemRepeat, right_flank: &[u8]
) -> Vec<Module> {
    let mut modules = Vec::new();
    modules.push((*tr_record.flank_l).into());
    for i in 0..tr_record.copy_unit.len() {
        modules.push((&tr_record.copy_unit[i][..], tr_record.copy_number[i]).into())
    }
    modules.push((*tr_record.flank_r).into());
    return modules;
}

pub(crate) fn print_to_file(json_str: &String, p: &Path) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(p)?;
    write!(file, "{}", json_str)?;
    return Ok(());
}

fn read_motifs(filename: &Path) -> Vec<(Vec<u8>, HGVSNomenclature, Vec<u8>)> {
    let file = File::open(filename).expect("Cannot find nomenclature file.");
    let reader = BufReader::new(file);

    // let crash = |_| panic!("line {}: Nomenclature {} malformatted. Accepted format is <chr>:g.<start>_<end><sequence>[repetitions].", i+1, split[1])
    // assert!(split.len() == 4,
    // "Malformatted line, expected format is <name>\\t<left_flank>\\t<hgvs_nomenclature>\\t<right_flank>\\n.");
    // Accepted format is <chr>:g.<start>_<end><sequence>[repetitions].

    let mut result = Vec::new();
    for line in reader.lines().skip(1) {
        let line = line.expect("Cannot read line from nomenclature file.").trim().to_owned();
        let split: Vec<_> = line.split('\t').collect();

        let name = split[0].to_owned();
        let left_flank = split[2].as_bytes().to_owned();
        let mut repeat: HGVSNomenclature = split[1].parse().expect("Malformatted nomenclature found.");
        repeat.name = Some(name);
        let right_flank = split[3].as_bytes().to_owned();

        result.push((left_flank, repeat, right_flank));
    }
    return result;
}

#[derive(Default, Debug, PartialEq, Clone)]
pub(crate) struct TRRecord {
    pub(crate) name: String,
    pub(crate) reference: String,
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) copy_unit: Vec<Vec<u8>>,
    pub(crate) copy_number: Vec<usize>,
    pub(crate) flank_l: Vec<u8>,
    pub(crate) flank_r: Vec<u8>
}

impl TRRecord {
    pub(crate) fn to_hgvs_nomenclature(&self) -> String {
        // TODO: output as TandemRepeat
        let x = HGVSNomenclature{
            name: Some(self.name.clone()),
            reference: self.reference.clone(),
            start: self.start,
            end: self.end,
            copy_unit: self.copy_unit.clone(),
            copy_number: self.copy_number.clone(),
        };
        return format!("{}", x);
    }

    pub(crate) fn region(&self) -> String {
        let region_str = format!("{}:{}-{}", self.reference, self.start + 1, self.end);
        return region_str;
    }
}

pub(crate) fn read_motifs2(filename: &Path) -> Vec<TRRecord> {
    let v = read_motifs(filename);
    let mut result = Vec::new();
    for (flank_l, tr, flank_r) in v {
        let x = TRRecord {
            name: tr.name.unwrap(),
            reference: tr.reference,
            start: tr.start,
            end: tr.end,
            copy_unit: tr.copy_unit,
            copy_number: tr.copy_number,
            flank_l,
            flank_r
        };
        result.push(x);
    }
    return result;
}
