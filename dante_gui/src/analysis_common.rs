use std::path::Path;
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::collections::HashSet;

pub(super) fn parse_motifs(path: &Path) -> Vec<(bool, String, Vec<String>, String)> {
    let file = File::open(path).expect("Cannot find motif file.");
    let reader = BufReader::new(file);

    let mut result = Vec::new();
    for line in reader.lines().skip(1) {
        let line = line.expect("Cannot read line from motif file.").trim().to_string();
        let split: Vec<_> = line.split('\t').collect();

        let id = split[0].to_string();
        let groups = split[4].split(',').map(|x| x.to_string()).collect();
        let description = split[5].to_string();
        result.push((false, id, groups, description));
    }

    return result;
}

pub(crate) fn get_groups(motifs: &[(bool, String, Vec<String>, String)]) -> Vec<(bool, String)> {
    let groups: HashSet<(bool, String)> = motifs.iter()
        .flat_map(|x| &x.2)
        .map(|x| (false, x.to_string()))
        .collect();
    let mut groups: Vec<_> = groups.into_iter().collect();
    groups.sort();
    return groups;
}
