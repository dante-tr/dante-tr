use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;
use std::error::Error;

use crate::repeats::TandemRepeat;
use crate::hmm::Module;

pub fn get_modules(
    left_flank: &[u8], repeat: &TandemRepeat, right_flank: &[u8]
) -> Vec<Module> {
    let mut modules = Vec::new();
    modules.push(left_flank.into());
    modules_add_motif(&mut modules, repeat);
    modules.push(right_flank.into());
    return modules;
}

fn modules_add_motif(modules: &mut Vec<Module>, motif: &TandemRepeat) {
    for i in 0..motif.copy_unit.len() {
        modules.push((&motif.copy_unit[i][..], motif.copy_number[i]).into())
    }
}

pub(crate) fn print_to_file(json_str: &String, p: &Path) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(p)?;
    write!(file, "{}", json_str)?;
    return Ok(());
}

pub(crate) fn read_motifs(filename: &Path) -> Vec<(Vec<u8>, TandemRepeat, Vec<u8>)> {
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
        let mut repeat: TandemRepeat = split[1].parse().expect("Malformatted nomenclature found.");
        repeat.name = Some(name);
        let right_flank = split[3].as_bytes().to_owned();

        result.push((left_flank, repeat, right_flank));
    }
    return result;
} 
