use std::collections::HashMap;
use std::collections::HashSet;

use crate::repeats::TandemRepeat;

fn resolve_inconsistency() {
    // r1 = reference used to create bam
    // r2 = input reference
    // r3 = reference, from which the hgvs was extracted

    // r3 must be subset of r2, otherwise we cannot extract flank
    // if r2 is not subset of r1, some motifes may have 0 reads mapped
}

fn hgvs_is_subset_of_reference(repeats: &[TandemRepeat], references: &HashMap<String, Vec<u8>>) -> bool {
    let tr_ids: HashSet<String> = repeats.iter().map(|x| x.reference.to_owned()).collect();
    let ref_ids: HashSet<String> = references.keys().map(|x| x.to_owned()).collect();
    tr_ids.is_subset(&ref_ids)
}

fn hgvs_wrt_ref_is_valid(repeats: &[TandemRepeat], references: &HashMap<String, Vec<u8>>) -> bool {
    for tr in repeats {
        let seq = match references.get(&tr.reference) {
            None => {
                println!("{} not found in reference.", tr.reference); 
                return false;
            }
            Some(s) => { s }
        };
        if tr.end > seq.len() { 
            println!("{}'s end is longer than reference sequence", tr);
            return false;
        }
    }
    return true;
}

fn hgvs_wrt_bam_is_valid(bam_refs: &HashMap<String, usize>, repeats: &[TandemRepeat]) -> bool {
    for tr in repeats {
        let len = match bam_refs.get(&tr.reference) {
            Some(n) => { *n },
            None => {
                println!("{} not found in bam.", &tr.reference);
                return false;
            }
        };
        if len < tr.end {
            println!("reference sequence is shorter than {}'s end", tr);
            return false;
        }
    }
    return true;
}

fn ref_wrt_bam_is_valid(bam_refs: &HashMap<String, usize>, references: &HashMap<String, Vec<u8>>) -> bool {
    for (id, seq) in references {
        let len = match bam_refs.get(id) {
            Some(n) => { *n },
            None => {
                println!("{} not found in bam.", id);
                return false;
            }
        };
        if len != seq.len() {
            println!("{} lengths in bam and fasta differ.", id);
            return false;
        }
    }
    return true;
}

// fn bam_wrt_ref_is_valid(bam_refs: &HashMap<String, usize>, references: &HashMap<String, Vec<u8>>) -> bool {
//     for (id, &len) in bam_refs {
//         let seq = match references.get(id) {
//             None => {
//                 println!("{} not found in reference.", id);
//                 return false;
//             },
//             Some(s) => { s }
//         };
//         if len != seq.len() {
//             println!("{} lengths differ in bam and fasta.", id);
//             return false;
//         }
//     } 
//     return true;
// }

fn fix_reference(
    references: HashMap<String, Vec<u8>>, bam_refs: &HashMap<String, usize>, mapping: &str
) -> HashMap<String, Vec<u8>> {
    let mut result = HashMap::new();

    let m = HashMap::from([("NC_000023.11", "chrX")]);
    let mut br = bam_refs.clone();
    for (id, seq) in references {
        let new_id = match m.get(&id[..]) {
            Some(new_id) => { new_id },
            None => { panic!(); } 
        };
        result.insert(new_id.to_string(), seq);
        br.remove(&new_id[..]);
    }
    return result;
}

fn fix_repeats(repeats: Vec<TandemRepeat>, mapping: &str) -> Vec<TandemRepeat> {
    let mut result = Vec::with_capacity(repeats.len());

    println!("Loading from {}", mapping);
    let m = HashMap::from([("NC_000023.11", "chrX")]);

    for tr in repeats {
        let new_id = match m.get(&tr.reference[..]) {
            Some(x) => { x.to_string() },
            None => { panic!(); }
        };
        result.push(TandemRepeat{reference: new_id, ..tr})
    }
    return result;
}

#[cfg(test)]
mod tests {
    use clap::Parser;
    use super::*;

    use crate::Args;
    use crate::read_nomenclature;
    use crate::read_reference;
    use crate::read_bam_refs;


    #[test]
    fn test_main() {
        let args = Args::try_parse_from([
            "remastr",
            "-f", "data/chromosomeX.fna",
            "-b", "data/real/twist_S22-157-01_S1.bam",
            "-n", "data/HGVS.txt",
            "--map1", "data/chromosomeX_to_bam_map.txt"
        ]).unwrap();

        let bam_refs = read_bam_refs(&args.bam_file);
        let mut references = read_reference(&args.ref_file);
        let mut repeats = read_nomenclature(&args.hgvs_file);
        // correct reference
        // correct repeats
        //

    }

    #[test]
    fn input_can_be_remapped() {
        let args = Args::try_parse_from([
            "remastr",
            "-f", "data/chromosomeX.fna",
            "-b", "data/real/twist_S22-157-01_S1.bam",
            "-n", "data/HGVS.txt",
            "--map1", "data/chromosomeX_to_bam_map.txt"
        ]).unwrap();

        let bam_refs = read_bam_refs(&args.bam_file);
        let mut references = read_reference(&args.ref_file);

        if ! ref_wrt_bam_is_valid(&bam_refs, &references) {
            println!("IDs in BAM and reference differ. Attempting correction of reference... ");
            if let Some(ref2bam) = args.ref2bam {
                println!("Correcting with map provided in {}.", ref2bam);
                references = fix_reference(references, &bam_refs, &ref2bam);
            } else {
                println!("Correcting with best effort heuristic.");
            }

            match ref_wrt_bam_is_valid(&bam_refs, &references) {
                true => { println!("Success!"); }
                false => { panic!("Unable to correct reference!"); }
            }
        }
    }

    #[test]
    fn nomenclature_can_be_remapped() {
        let args = Args::try_parse_from([
            "remastr",
            "-f", "data/chromosomeX.fna",
            "-b", "data/real/twist_S22-157-01_S1.bam",
            "-n", "data/HGVS.txt",
            "--map1", "data/chromosomeX_to_bam_map.txt"
        ]).unwrap();

        let bam_refs = read_bam_refs(&args.bam_file);
        let mut repeats = read_nomenclature(&args.hgvs_file);

        if ! hgvs_wrt_bam_is_valid(&bam_refs, &repeats) {
            println!("IDs in BAM and nomenclature differ. Attempting correction of nomenclature... ");
            if let Some(hgvs2bam) = args.hgvs2bam {
                println!("Correcting with map provided in {}.", hgvs2bam);
            } else if let Some(ref2bam) = args.ref2bam {
                println!("Correcting with map provided in {}.", ref2bam);
                repeats = fix_repeats(repeats, &ref2bam);
            } else {
                println!("Correcting with best effort heuristic.");
            }

            match hgvs_wrt_bam_is_valid(&bam_refs, &repeats) {
                true => { println!("Success!"); }
                false => { panic!("Unable to correct reference!"); }
            }
        }
    }
}
