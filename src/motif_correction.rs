use std::collections::HashMap;
use std::str;
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
            let from = motif.start - FLANK;
            let to = motif.end + FLANK;
            let seq = ref_region(refs, &motif.reference, from, to)
                .expect("Unable to get reference region.");

            let corrected_motif = correct_motif(&seq, &motif, FLANK);
            println!(
                "{} -> {}\n{}\n{}\n{}\n",
                motif, corrected_motif,
                str::from_utf8(&seq).unwrap(),
                str::from_utf8(&motif.view(from, to)).unwrap(),
                str::from_utf8(&corrected_motif.view(from, to)).unwrap(),
            );

            valid_repeats.push(corrected_motif);
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

fn fn3(model: &HMM, annotation: &[usize]) -> (usize, usize) {
    let start: usize = 0;
    let end: usize = 6; //TODO: model.get_end();

    let mut m_start = usize::MIN;
    let mut m_end = usize::MAX;
    for (i, &state) in annotation.iter().enumerate() {
        if state == start { m_start = i; }
        if state == end && m_end == usize::MAX { m_end = i; }
    }
    return (m_start + 1, m_end);
}

#[cfg(test)]
mod tests {
    use super::*;

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
}

