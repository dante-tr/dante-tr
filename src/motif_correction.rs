use std::collections::HashMap;
use std::str;
use ndarray::ArrayView1;
use ndarray::{Array, Array2};
use ndarray::s;
use crate::repeats::TandemRepeat;
use crate::modules_add_motif;
use crate::HMM;

const FLANK:usize = 20;

pub fn correct_repeats(refs: &HashMap<String, Vec<u8>>, repeats: &Vec<TandemRepeat>) -> Vec<TandemRepeat> {
    let mut valid_repeats = Vec::new();
    for motif in repeats.iter() {
        if is_present(&motif, &refs) {
            valid_repeats.push(motif.clone());
        } else {
            // let from = motif.start - FLANK;
            // let to = motif.end + FLANK;
            // let seq = ref_region(refs, &motif.reference, from, to)
            //     .expect("Unable to get reference region.");

            // let corrected_motif = correct_motif(&seq, &motif, FLANK);
            // println!(
            //     "{} -> {}\n{}\n{}\n{}\n",
            //     motif, corrected_motif,
            //     str::from_utf8(&seq).unwrap(),
            //     str::from_utf8(&motif.view(from, to)).unwrap(),
            //     str::from_utf8(&corrected_motif.view(from, to)).unwrap(),
            // );

            // valid_repeats.push(corrected_motif);
        }
    }
    return valid_repeats;
}

fn ref_region<'a>(
    refseq: &'a HashMap<String, Vec<u8>>, id: &str, start: usize, end: usize
) -> Option<&'a[u8]> {
    let seq = match refseq.get(id) {
        None => { return None; },
        Some(x) => { x },
    };
    return Some(&seq[start..end]);
}

fn is_present(tr: &TandemRepeat, seq: &HashMap<String, Vec<u8>>) -> bool {
    let ref_repeat = match ref_region(seq, &tr.reference, tr.start, tr.end) {
        None => { return false; },
        Some(x) => { x },
    };
    let hgvs_repeat = &tr.sequence();
    if ref_repeat != hgvs_repeat {
        return false;
    }
    return true;
}

fn correct_motif(seq: &[u8], repeat: &TandemRepeat, flank: usize) -> TandemRepeat {
    let qual = b"~".repeat(seq.len());

    let mut modules = Vec::new();
    modules_add_motif(&mut modules, &repeat); // TODO: make this function?

    let model = HMM::from(&modules).log();
    let (_, annotation) = model.log_predict(&seq, &qual);

    let suggested_repeat = {
        let mut new_repeat = repeat.clone();
        let start = match annotation.iter().position(|&x| x != 0) {
            None => {
                eprintln!("Unable to match with reference.");
                return new_repeat;
            },
            Some(x) => { x }
        };
        let start = repeat.start - flank + start;
        new_repeat.start = start;
        new_repeat.end = start + repeat.sequence().len();
        new_repeat
    };

    let mut orig_motif = b"-".repeat(flank);
    orig_motif.extend_from_slice(&repeat.sequence());
    orig_motif.extend_from_slice(&b"-".repeat(flank));

    return suggested_repeat;
}

const INDEL:u8 = 1;

fn fill_dp_table(target: &[u8], query: &[u8]) -> Array2<u8> {
    let n = query.len() + 1;
    let m = target.len() + 1;
    let mut dp = Array::zeros((n, m));

    for i in 0..n { dp[[i, 0]] = i as u8; }
    for j in 0..m { dp[[0, j]] = 0; }

    for i in 1..n {
        for j in 1..m {
            let edit = (query[i-1] != target[j-1]) as u8;
            dp[[i, j]] = *[
                dp[[i-1, j-1]] + edit,
                dp[[i-1, j]] + INDEL,
                dp[[i, j-1]] + INDEL
            ].iter().min().unwrap();
        }
    }
    return dp;
}

use counter::Counter;
fn sgalign(target: &[u8], query: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let dp = fill_dp_table(target, query);
    let cigar = get_cigar(dp, target, query);
    let count: Counter<_> = cigar.iter().collect();
    println!("{:?}", count);
    let (aligned_target, aligned_query) = apply_cigar(&cigar, target, query);
    return (aligned_target, aligned_query);
}

fn apply_cigar(cigar: &[Cigar], target: &[u8], query: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut new_target = Vec::new();
    let mut new_query = Vec::new();

    let mut i = 0;
    let mut j = 0;
    for c in cigar {
        match c {
            Cigar::M | Cigar::X => {
                new_target.push(target[i]); i += 1;
                new_query.push(query[j]); j += 1;
            },
            Cigar::D => {
                new_target.push(b'_');
                new_query.push(query[j]); j += 1;
            },
            Cigar::I => {
                new_target.push(target[i]); i += 1;
                new_query.push(b'_');
            },
            Cigar::N => {
                new_target.push(target[i]); i += 1;
                new_query.push(b'-');
            }
        }
    }
    return (new_target, new_query);
}

#[derive(Debug, PartialEq, Eq, Hash)]
enum Cigar {
    M,  // match
    X,  // substitution
    I,  // insertion to target
    D,  // deletion from target
    N   // clipping in target
}

fn argmin(a: ArrayView1<u8>) -> usize {
    let mut min_pos = 0;
    let mut min_val = u8::MAX;
    for (p, &v) in a.iter().enumerate() {
        if v < min_val { min_pos = p; min_val = v; }
    }
    return min_pos;
}

fn get_cigar(a: Array2<u8>, target: &[u8], query: &[u8]) -> Vec<Cigar> {
    let n = a.shape()[0];
    let m = a.shape()[1];
    let end = argmin(a.slice(s![n-1, ..]));

    let mut cigar = Vec::new();
    for _ in 0..m-1-end { cigar.push(Cigar::N); }

    let mut i = n-1;
    let mut j = end;
    while i != 0 {
        // if j == 0: push D and i--; but first write failing test;
        let edit = (query[i-1] != target[j-1]) as u8;

        if a[[i, j]] == a[[i-1, j]] + INDEL {
            cigar.push(Cigar::D);
            i -= 1;
        } else if a[[i, j]] == a[[i-1, j-1]] + edit {
            if edit == 0 { cigar.push(Cigar::M); }
            else { cigar.push(Cigar::X); }
            i -= 1; j -= 1;
        } else if a[[i, j]] == a[[i, j-1]] + INDEL {
            cigar.push(Cigar::I);
            j -= 1;
        }
    }
    for _ in 0..j { cigar.push(Cigar::N); }
    cigar.reverse();
    return cigar;
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn can_move_motif() {
        let motif: TandemRepeat = "SEQ1:g.6_15CG[5]".parse().unwrap();
        let flank = 5;
        let from = motif.start - flank;
        let to = motif.end + flank;
        let seq = &b"AAAAAAACGCGCGCGCGAAA"[from..to];

        let expected_motif: TandemRepeat = "SEQ1:g.8_17CG[5]".parse().unwrap();
        let corrected_motif = correct_motif(&seq, &motif, flank);

        println!(
            "{} -> {}\n{}\n{}\n{}",
            motif, corrected_motif,
            str::from_utf8(&seq).unwrap(),
            str::from_utf8(&motif.view(from, to)).unwrap(),
            str::from_utf8(&corrected_motif.view(from, to)).unwrap(),
        );

        assert_eq!(expected_motif, corrected_motif);
    }

    #[test]
    fn test_semiglobal_align() {
        let mut motif: TandemRepeat = "S1:g.1_18AACCCT[3]".parse().unwrap();
        let reference = b"TGTAACCCGAAACCTCAAAGCCTAACCCTAACCCTAACCCCTACAGTTGAGGTCCCCC".to_vec();

        let seq_local = reference;

        for i in 3..10 {
            motif.copy_number[0] = i;
            let seq_global = motif.sequence();
            let (s1, s2) = sgalign(&seq_local, &seq_global);
            println!("{}", str::from_utf8(&s1).unwrap());
            println!("{}", str::from_utf8(&s2).unwrap());
            println!();
        }
    }

    #[test]
    fn test_semiglobal_align_simple() {
        let reference = b"ACCCA".to_vec();
        let query = b"CCC".to_vec();
        let (s1, s2) = sgalign(&reference, &query);
        println!("{}", str::from_utf8(&s1).unwrap());
        println!("{}", str::from_utf8(&s2).unwrap());
    }

    #[test]
    fn test_dp_to_cigar() {
        let target = b"ACCCA".to_vec();
        let query = b"CCC".to_vec();
        let dp = array![
            [0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 1],
            [2, 2, 1, 0, 0, 1],
            [3, 3, 2, 1, 0, 1]
        ];

        use Cigar as C;
        let exp_cigar = vec![C::N, C::M, C::M, C::M, C::N];
        let cigar = get_cigar(dp, &target, &query);
        assert_eq!(exp_cigar, cigar);
    }

    #[test]
    fn test_cigar_with_deletion() {
        let target = b"GGGGAACCCCTGGGG".to_vec();
        let query = b"AACCCT".to_vec();
        let dp = fill_dp_table(&target, &query);
        println!("{:?}", dp);
        let cigar = get_cigar(dp, &target, &query);
        println!("{:?}", cigar);

        // use Cigar as C;
        // let exp_cigar = vec![NNNNMMMMMDMNNNN];
    }
}

