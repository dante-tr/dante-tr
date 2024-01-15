use std::collections::HashMap;
use std::str;

#[cfg(test)]
use ndarray::{Array, Array2};
#[cfg(test)]
use ndarray::s;
#[cfg(test)]
use counter::Counter;
#[cfg(test)]
use ndarray::ArrayView1;
#[cfg(test)]
use ndarray::ArrayView2;

use crate::repeats::TandemRepeat;
// use crate::modules_add_motif;
// use crate::HMM;

// const FLANK:usize = 20;
#[cfg(test)]
type DPEntry = u16;
#[cfg(test)]
const INDEL:DPEntry = 1;


#[test]
fn can_move_motif() {
    let motif: TandemRepeat = "SEQ1:g.6_15CG[5]".parse().unwrap();
    let flank = 5;
    let from = motif.start - flank;
    let to = motif.end + flank;
    let seq = &b"AAAAAAACGCGCGCGCGAAA"[from..to];

    let expected_motif: TandemRepeat = "SEQ1:g.8_17CG[5]".parse().unwrap();
    let corrected_motif = correct_motif(seq, &motif, flank);

    println!(
        "{} -> {}\n{}\n{}\n{}",
        motif, corrected_motif,
        str::from_utf8(seq).unwrap(),
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

    let mut start_counter: Counter<usize, usize> = Counter::new();
    let mut end_counter: Counter<usize, usize> = Counter::new();
    for i in 3..10 {
        motif.copy_number[0] = i;
        let seq_global = motif.sequence();
        let dp = fill_dp_table(&seq_local, &seq_global);
        let cigar = get_cigar(dp.view(), &seq_local, &seq_global);
        let start_offset = get_start_offset(&cigar);
        start_counter[&start_offset] += 1;
        let end_offset = get_end_offset(&cigar);
        end_counter[&end_offset] += 1;
        let count: Counter<_> = cigar.iter().collect();
        let (aligned_target, aligned_query) = apply_cigar(&cigar, &seq_local, &seq_global);

        println!("{:?}\t{}\t{}", count, start_offset, end_offset);
        println!("{}", str::from_utf8(&aligned_target).unwrap());
        println!("{}", str::from_utf8(&aligned_query).unwrap());
        println!();
    }
    let n = seq_local.len();
    let recommended_start = start_counter.most_common()[0].0;
    let recommended_end = n - 1 - end_counter.most_common()[0].0;
    println!("{:?}", recommended_start);
    println!("{:?}", recommended_end);

    motif.start = recommended_start;
    motif.end = recommended_end;
    motif.copy_number[0] = (motif.end - recommended_start) / motif.copy_unit[0].len();
    println!("{:?}", motif);

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
    use ndarray::array;
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
    let cigar = get_cigar(dp.view(), &target, &query);
    assert_eq!(exp_cigar, cigar);
}

#[test]
fn test_cigar_with_deletion() {
    let target = b"GGGGAACCCCTGGGG".to_vec();
    let query = b"AACCCT".to_vec();
    let dp = fill_dp_table(&target, &query);
    println!("{:?}", dp);
    let cigar = get_cigar(dp.view(), &target, &query);
    println!("{:?}", cigar);

    // use Cigar as C;
    // let exp_cigar = vec![NNNNMMMMMDMNNNN];
}

#[test]
fn look_at_motif() {
    let target = b"NNNNNNNNNNNNNNNNNNNNCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTCTGAAAGTGGACCTATCA";
    let old_motif: TandemRepeat = "NC_000023.11:g.10001_10036AACCCT[6]".parse().unwrap();
    let new_motif: TandemRepeat = "NC_000023.11:g.10003_10040AACCCT[6]".parse().unwrap();

    let (s1, s2) = sgalign(&target[..], &old_motif.sequence());
    println!("{}", str::from_utf8(&s1).unwrap());
    println!("{}", str::from_utf8(&s2).unwrap());

    let (s1, s2) = sgalign(&target[..], &new_motif.sequence());
    println!("{}", str::from_utf8(&s1).unwrap());
    println!("{}", str::from_utf8(&s2).unwrap());
}

pub fn correct_repeats(refs: &HashMap<String, Vec<u8>>, repeats: &[TandemRepeat]) -> Vec<TandemRepeat> {
    println!("Correcting motifs.");
    let mut valid_repeats = Vec::new();
    for motif in repeats.iter() {
        if is_present(motif, refs) {
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

#[cfg(test)]
fn sgalign(target: &[u8], query: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let dp = fill_dp_table(target, query);
    let cigar = get_cigar(dp.view(), target, query);
    let count: Counter<_> = cigar.iter().collect();
    println!("{:?}", count);
    let (aligned_target, aligned_query) = apply_cigar(&cigar, target, query);
    return (aligned_target, aligned_query);
}

#[cfg(test)]
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

#[cfg(test)]
#[derive(Debug, PartialEq, Eq, Hash)]
enum Cigar {
    M,  // match
    X,  // substitution
    I,  // insertion to target
    D,  // deletion from target
    N   // clipping in target
}

#[cfg(test)]
fn correct_motif(seq: &[u8], repeat: &TandemRepeat, flank: usize) -> TandemRepeat {
    let mut motif = repeat.clone();

    let unit_len = motif.copy_unit[0].len();
    let max_k = seq.len() / unit_len;
    motif.copy_number[0] = max_k;
    let query = motif.sequence();
    let dp = fill_dp_table(seq, &query);

    let mut start_counter: Counter<usize, usize> = Counter::new();
    let mut end_counter: Counter<usize, usize> = Counter::new();

    for i in 1..=max_k {
        let dp_view = dp.slice(s![0..(unit_len * i) + 1, ..]);
        let cigar = get_cigar(dp_view, seq, &query);

        let start_offset = get_start_offset(&cigar);
        start_counter[&start_offset] += 1;
        let end_offset = get_end_offset(&cigar);
        end_counter[&end_offset] += 1;
    }
    let recommended_start = motif.start - flank + start_counter.most_common()[0].0;
    let recommended_end = motif.end + flank - end_counter.most_common()[0].0;

    motif.start = recommended_start;
    motif.end = recommended_end;
    motif.copy_number[0] = (motif.end - motif.start) / unit_len;
    return motif;
}

#[cfg(test)]
fn fill_dp_table(target: &[u8], query: &[u8]) -> Array2<DPEntry> {
    let n = query.len() + 1;
    let m = target.len() + 1;
    let mut dp = Array::zeros((n, m));

    for i in 0..n { dp[[i, 0]] = i as DPEntry; }
    for j in 0..m { dp[[0, j]] = 0; }

    for i in 1..n {
        for j in 1..m {
            let edit = (query[i-1] != target[j-1]) as DPEntry;
            dp[[i, j]] = *[
                dp[[i-1, j-1]] + edit,
                dp[[i-1, j]] + INDEL,
                dp[[i, j-1]] + INDEL
            ].iter().min().unwrap();
        }
    }
    return dp;
}

#[cfg(test)]
fn get_cigar(a: ArrayView2<DPEntry>, target: &[u8], query: &[u8]) -> Vec<Cigar> {
    let n = a.shape()[0];
    let m = a.shape()[1];
    let end = argmin(a.slice(s![n-1, ..]));

    let mut cigar = Vec::new();
    for _ in 0..m-1-end { cigar.push(Cigar::N); }

    let mut i = n-1;
    let mut j = end;
    while i != 0 {
        // if j == 0: push D and i--; but first write failing test;
        let edit = (query[i-1] != target[j-1]) as DPEntry;

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
fn argmin(a: ArrayView1<DPEntry>) -> usize {
    let mut min_pos = 0;
    let mut min_val = DPEntry::MAX;
    for (p, &v) in a.iter().enumerate() {
        if v < min_val { min_pos = p; min_val = v; }
    }
    return min_pos;
}

#[cfg(test)]
fn get_start_offset(cigar: &[Cigar]) -> usize {
    for (offset, c) in cigar.iter().enumerate() {
        if ! matches!(c, Cigar::N) { return offset; }
    }
    return usize::MAX;
}

#[cfg(test)]
fn get_end_offset(cigar: &[Cigar]) -> usize {
    for (offset, c) in cigar.iter().rev().enumerate() {
        if ! matches!(c, Cigar::N) { return offset; }
    }
    return usize::MAX;
}

