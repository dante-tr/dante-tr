use noodles::fasta;
use std::collections::HashMap;
use std::str;

use crate::TandemRepeat;
use crate::Module;

pub fn read_reference(filename: &str) -> HashMap<String, Vec<u8>> {
    let mut reader = fasta::reader::Builder.build_from_path(filename).unwrap();

    let mut result = HashMap::new();
    for record in reader.records() {
        let record = record.unwrap();

        result.insert(record.name().to_string(), (record.sequence()[..]).to_vec());
        // Is there a better way to get Vec<u8> than this? --------^
        // Do I need Vec<u8>? Cannot I leave it as Sequence?
    }
    return result;
}

pub fn get_modules(
    repeat: &TandemRepeat, refs: &HashMap<String, Vec<u8>>, flank_size: usize
) -> Vec<Module> {
    let refseq = refs.get(&repeat.reference).unwrap(); // safe due to nomenclature check
    assert!(repeat.start >= flank_size,
        "Cannot create left flank of size {flank_size} for repeat {repeat}.");
    let left_flank = &refseq[(repeat.start-flank_size)..repeat.start];
    assert!(repeat.end + flank_size <= refseq.len(),
        "Cannot create right flank of size {flank_size} for repeat {repeat}.");
    let right_flank = &refseq[repeat.end..(repeat.end+flank_size)];

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

