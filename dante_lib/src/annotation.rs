use std::fmt;
use std::iter::zip;
use std::ops::Range;
use std::str;

use itertools::izip;
use polars::prelude::DataFrame;

use nom::AsBytes;
use noodles::bam;

use crate::df_ops;
use crate::io::TRRecord;
use crate::hmm::Hmm;

pub fn annotate<T>(in_reads: T, model: Hmm, repeat: &TRRecord) -> DataFrame
where
    T: Iterator<Item = bam::Record>,
{
    let mut seqs = Vec::new();
    let mut quals = Vec::new();
    let mut read_ids = Vec::new();
    let mut mate_orders = Vec::new();
    for read in in_reads {
        let seq: Vec<_> = read.sequence().iter().collect();
        let qual: Vec<_> = read.quality_scores().as_ref().iter().map(|x| x + 33).collect();
        let read_id = str::from_utf8(read.name().unwrap().as_bytes()).unwrap();
        let mate_order = mate_order(&read);
        if qual.is_empty() {
            println!("Read {read_id} does not have sequence and quality.");
            // should I do something else?
            continue;
        }
        seqs.push(seq);
        quals.push(qual);
        read_ids.push(read_id.to_string());
        mate_orders.push(mate_order);
    }

    let n = seqs.len();

    let mut read_sns = Vec::new();
    for i in 0..n {
        let read_sn = i;
        read_sns.push(read_sn as u64);
    }

    let mut names = Vec::new();
    for _ in 0..n {
        let name: String = repeat.name.clone();
        // let name: String = match &repeat.name {
        //     Some(x) => x.to_string(),
        //     None    => "None".to_string(),
        // };
        names.push(name);
    }

    let mut motifs = Vec::new();
    for _ in 0..n {
        let motif = repeat.to_hgvs_nomenclature();
        motifs.push(motif);
    }

    let mut qualities = Vec::new();
    for _ in 0..n {
        let qual_str: Vec<_> = "No quality".bytes().collect();
        let quality = str::from_utf8(&qual_str).unwrap();
        qualities.push(quality.to_string());
    }

    let mut reads = Vec::new();
    let mut references = Vec::new();
    let mut moduleses = Vec::new();
    let mut log_likelihoods = Vec::new();
    let mut left_bgs = Vec::new();
    let mut right_bgs = Vec::new();
    let mut n_moduleses = Vec::new();
    let mut module_baseses = Vec::new();
    let mut module_repetitionses = Vec::new();
    let mut module_sequenceses = Vec::new();
    let mut module_nomenclatureses = Vec::new();
    let mut mismatches_strs = Vec::new();
    let mut module_classeses = Vec::new();
    for (seq, qual) in izip!(seqs, quals) {

        let (log_likelihood, annotation) = model.log_predict(&seq, &qual);
        let (new_annot, reconstructed_read) = model.realign(&annotation, &seq);
        let reconstructed_reference = model.reconstruct_sequence(&new_annot);

        let (partition, mod_ids) = model.partition_to_units(&annotation);
        let mods = model.reconstruct_mod_ids(&new_annot);

        reads.push(str::from_utf8(&reconstructed_read).unwrap().to_string());
        references.push(str::from_utf8(&reconstructed_reference).unwrap().to_string());
        moduleses.push(str::from_utf8(&mods).unwrap().to_string());
        mismatches_strs.push(generate_mismatches(&reconstructed_read, &reconstructed_reference));
        log_likelihoods.push(log_likelihood);

        let left_bg = get_left_bg(&mods);
        let right_bg = get_right_bg(&mods);
        let n_modules = repeat.copy_number.len() + 2;

        let mut module_sequences: Vec<String> = Vec::with_capacity(n_modules);
        for i in 0..n_modules {
            let ms = get_module_sequences(&seq, &partition, &mod_ids, i);
            module_sequences.push(ms);
        }

        let mut module_bases: Vec<usize> = Vec::with_capacity(n_modules);
        for i in 0..n_modules {
            let mb = get_module_bases(&mods, i);
            module_bases.push(mb);
        }

        let mut module_nomenclatures: Vec<String> = Vec::with_capacity(n_modules);
        for i in 0..n_modules {
            let mn = get_module_nomenclature(&seq, &partition, &mod_ids, i);
            module_nomenclatures.push(mn);
        }

        let mut module_repetitions: Vec<usize> = Vec::with_capacity(n_modules);
        for i in 0..n_modules {
            let mr = get_module_repetitions(&mod_ids, i);
            module_repetitions.push(mr);
        }

        let module_classes = get_module_classes(left_bg, &module_bases, right_bg);
        let module_bases = module_bases.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        let module_repetitions = module_repetitions.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        let module_sequences = module_sequences.join(",");
        let module_nomenclatures = module_nomenclatures.join(",");
        let module_classes = module_classes.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");

        n_moduleses.push(n_modules as u64);
        left_bgs.push(left_bg as u64);
        right_bgs.push(right_bg as u64);
        module_baseses.push(module_bases);
        module_repetitionses.push(module_repetitions);
        module_sequenceses.push(module_sequences);
        module_nomenclatureses.push(module_nomenclatures);
        module_classeses.push(module_classes);
    }

    let mut n_deletionses = Vec::new();
    let mut n_insertionses = Vec::new();
    let mut n_mismatcheses = Vec::new();
    for mismatches_str in &mismatches_strs {
        let n_deletions = mismatches_str.bytes().filter(|x| *x == b'D').count();
        n_deletionses.push(n_deletions as u64);
        let n_insertions = mismatches_str.bytes().filter(|x| *x == b'I').count();
        n_insertionses.push(n_insertions as u64);
        let n_mismatches = mismatches_str.bytes().filter(|x| *x == b'M').count();
        n_mismatcheses.push(n_mismatches as u64);
    }

    // make this a builder pattern?
    let result = df_ops::construct_df(
        names,                    // "ALS"
        motifs,                   // "chr15:g.22786680_22786703GGC[8]"
                                  // motif modules
        read_sns,                 // 0
        read_ids,                 // "HISEQ1:29:HA2WPADXX:2:2202:2985:13224"
        mate_orders,              // "1"
                                  // TODO: add seq?
        qualities,                // "No quality" TODO: add qual?
        reads,                    // "CCTCTTCCTGCTCCTCCCCCACCCGTCCCCCTCCCCTCCCCCGCCCGCGCCTCCCGGTCACCCCCCATCCCGCCCCGCGGGGCGCGGCGCGCAGGCGCAGGCTCGGAGGGCGGGCGCGGGCGGAATGGGGACTGCAGCTGCGGCAGCG"
        references,               // "---------------------------------------------------------------------------------------------------------------------GGGCGGAATGGGGACTGCAGCTGCGGCAGCG"
        moduleses,                // "---------------------------------------------------------------------------------------------------------------------0000000000000000000000000000001"
        mismatches_strs,          // "____________________________________________________________________________________________________________________________________________________"
        log_likelihoods,          // -172.339767
        left_bgs,                 // 117
        right_bgs,                // 0
        n_deletionses,            // 0
        n_insertionses,           // 0
        n_mismatcheses,           // 0
        n_moduleses,              // 3
        module_baseses,           // "30,1,0"
        module_repetitionses,     // "1,0,0"
        module_sequenceses,       // "GGGCGGAATGGGGACTGCAGCTGCGGCAGC,G,"
        module_nomenclatureses,   // "GGGCGGAATGGGGACTGCAGCTGCGGCAGC[1],G[1],"
        module_classeses,         // "Flanking,Missing,Missing"
    ).expect("Cannot create dataframe");
    return result;
}

fn get_left_bg(mods: &[u8]) -> usize {
    let mlen = mods.len();
    let mut left_bg = 0;
    while left_bg < mlen && mods[left_bg] == b'-' { left_bg += 1; }
    return left_bg;
}

fn get_right_bg(mods: &[u8]) -> usize {
    let mlen = mods.len();
    let mut right_bg = 0;
    while right_bg < mlen && mods[(mlen - 1) - right_bg] == b'-' { right_bg += 1; }
    return right_bg;
}

#[derive(Clone, PartialEq)]
enum AClass {
    Spanning,
    Flanking,
    InRepeat, // In-repeat
    Missing,
    Filtered(String)
}

impl fmt::Display for AClass {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AClass::Missing => write!(f, "Missing"),
            AClass::Spanning => write!(f, "Spanning"),
            AClass::Flanking => write!(f, "Flanking"),
            AClass::InRepeat => write!(f, "In-repeat"),
            AClass::Filtered(x) => write!(f, "Filtered({x})")
        }
    }
}

fn get_module_classes(left_bg: usize, module_bases: &[usize], right_bg: usize) -> Vec<AClass> {
    let mut base_count = Vec::with_capacity(module_bases.len() + 2);
    base_count.push(left_bg);
    for &x in module_bases { base_count.push(x); }
    base_count.push(right_bg);

    // Filtering options:
    // 1)  addresses this case:
    //     modules               0111111111111222222222222222222222222222222...
    //     module_bases          1,12,33,30
    //     module_repetitions    1,4,11,1
    const MIN_MOD_LEN: usize = 3;
    // for x in &mut base_count { *x = (*x).saturating_sub(MIN_MOD_LEN); }
    // 2)  previous case could be better addressed by ignoring first and last ~3bp of read/annotation
    let mut to_remove = MIN_MOD_LEN;
    let mut i = 0;
    while to_remove > 0 {
        let tmp = [to_remove, base_count[i]];
        let m = *tmp.iter().min().unwrap();
        to_remove -= m; base_count[i] -= m;
        i += 1;
    }

    let mut to_remove = MIN_MOD_LEN;
    let mut i = base_count.len() - 1;
    while to_remove > 0 {
        let tmp = [to_remove, base_count[i]];
        let m = *tmp.iter().min().unwrap();
        to_remove -= m; base_count[i] -= m;
        i -= 1;
    }
    // 3)  if read has too many mismatches+indels, filter it out.
    //     This was implemented, but unused in the python code.
    // 4)  Filter out skipped modules. This should not happen now,
    //     but filtering option 1) can trigger it
    let first_nonzero = base_count.iter().position(|&x| x != 0);
    let last_nonzero  = base_count.iter().rposition(|&x| x != 0);
    let valid = match (first_nonzero, last_nonzero) {
        (Some(s), Some(e)) => base_count[s..=e].iter().all(|&x| x != 0),
        (None, None) => true,
        _ => unreachable!()
    };
    if ! valid { return vec![AClass::Filtered("Incorrect Annotation".to_string()); module_bases.len()] }
    // End of filtering options

    let mut result = Vec::with_capacity(module_bases.len());
    for i in 1..(base_count.len()-1) {
        if base_count[i] == 0 {
            result.push(AClass::Missing);
        } else if base_count[i-1] != 0 && base_count[i+1] != 0 {
            result.push(AClass::Spanning);
        } else {
            // bc[i] != 0 and (bc[i-1] == 0 or bc[i+1] == 0) 
            result.push(AClass::Flanking);
        }
    }
    for i in 1..(result.len()-1) {
        if result[i-1] == AClass::Missing && result[i] == AClass::Flanking && result[i+1] == AClass::Missing {
            result[i] = AClass::InRepeat;
        }
    }
    return result;
}

fn get_module_sequences(seq: &[u8], partition: &[Range<usize>], mod_ids: &[usize], idx: usize) -> String {
    let mut ms = Vec::new();
    for i in 0..mod_ids.len() {
        if mod_ids[i] == idx {
            let x = &seq[partition[i].clone()];
            ms.extend(x);
        }
    }
    let ms = str::from_utf8(&ms).unwrap().to_string();
    return ms;
}

fn get_module_nomenclature(seq: &[u8], partition: &[Range<usize>], mod_ids: &[usize], idx: usize) -> String {
    let mut module_nomenclature = Vec::new();
    let mut append_unit = |s, o|{
        if o != 0 {
            module_nomenclature.push(format!("{}[{}]", str::from_utf8(s).unwrap(), o));
        }
    };

    let mut prev: &[u8] = b"";
    let mut occ = 0;
    for i in 0..mod_ids.len() {
        if mod_ids[i] == idx {
            let curr = &seq[partition[i].clone()];
            if curr == prev {
                occ += 1;
            } else {
                append_unit(prev, occ);

                prev = curr;
                occ = 1;
            }
        }
    }
    append_unit(prev, occ);
    return module_nomenclature.join("");
}

fn get_module_repetitions(mod_ids: &[usize], idx: usize) -> usize {
    let mut result = 0;
    for &m_id in mod_ids {
        if m_id == idx {
            result += 1;
        }
    }
    return result;
}

fn get_module_bases(mods: &[u8], idx: usize) -> usize {
    const ASCII_ZERO: usize = 48;
    let idx: u8 = (idx + ASCII_ZERO).try_into().unwrap();
    let count = mods.iter().filter(|&&x| x == idx).count();
    return count;
}

fn generate_mismatches(read: &[u8], reference: &[u8]) -> String {
    let mut result = String::with_capacity(read.len());
    for (x, y) in zip(read, reference) {
        match (x, y) {
            (_,    b'-') => { result.push('_'); }
            (b'_', _   ) => { result.push('D'); }
            (_,    b'_') => { result.push('I'); }
            (x,    y   ) => {
                if x == y { result.push('_'); } else { result.push('M'); }
            }
        }
    }
    return result;
}

fn mate_order(read: &bam::Record) -> String {
    if read.flags().is_first_segment() { "1".to_string() }
    else if read.flags().is_last_segment() { "2".to_string() } 
    else {
        // println!("Read {} does not have pair information.", read.read_name().unwrap());
        "0".to_string()
    }
}

// #[test]
// fn tmp_fn_name() {
//     use crate::io::get_modules;
//     // CTTCCTCCTCCTCATCGGTGGCGGCGGCGGCGGCGTCAGGCCAGTGCCGCGGCTTTCTCTCCGCGCCTGTGTTCGCCGGGACGCATTCGGGGCGGGCGGCGGCGGCGGCAGCGGCGGCCGCGGCGGCGGGCGGGGCCGCCCCCCGCCT
//     // CCCFFFFFHHHHGIIJJIHIJIJJJHGHDDDBDDBBDDDDDDDDCCCC5@BDDBBCCCCACDBBDBDBBCCCCCBDDBBD@BDDBDDDDDDDDDD@5.5&)0)&5))0)&)2&55)0&&&&&)&&)&))&&&&&)&&&&)&&50)&&&
//     // -----------------------------------------------------------------00000000000000000000000000000011111111111111111111111111111111111111111111122222222
//     // GCG[4]GCA[1]GCG[2]GCC[1]GCG[3]G[1]GCG[1]GGGCCGCC[1]
//     //
//     // I think annotation is wrong
//     // Independent of annotation, seq nomenclature is wrong as well
//     //
//     // SPD     chr2:g.176093059_176093103GCG[15]       CCTGTGTTCGCCGGGACGCATTCGGGGCGG  TCCGGCTTTGCGTACCCCGGGACCTCTGAG
//     // result 15
// 
//     // HISEQ1:26:HA2RRADXX:1:1203:16720:7919
//     let seq:  Vec<u8> = b"CTTCCTCCTCCTCATCGGTGGCGGCGGCGGCGGCGTCAGGCCAGTGCCGCGGCTTTCTCTCCGCGCCTGTGTTCGCCGGGACGCATTCGGGGCGGGCGGCGGCGGCGGCAGCGGCGGCCGCGGCGGCGGGCGGGGCCGCCCCCCGCCT".to_vec();
//     let qual: Vec<u8> = b"CCCFFFFFHHHHGIIJJIHIJIJJJHGHDDDBDDBBDDDDDDDDCCCC5@BDDBBCCCCACDBBDBDBBCCCCCBDDBBD@BDDBDDDDDDDDDD@5.5&)0)&5))0)&)2&55)0&&&&&)&&)&))&&&&&)&&&&)&&50)&&&".to_vec();
// 
//     // SPD
//     let left_flank:  Vec<u8> = b"CCTGTGTTCGCCGGGACGCATTCGGGGCGG".to_vec();
//     let right_flank: Vec<u8> = b"TCCGGCTTTGCGTACCCCGGGACCTCTGAG".to_vec();
//     let repeat: TandemRepeat = "chr2:g.176093059_176093103GCG[15]".parse().expect("Malformatted nomenclature found.");
// 
//     let modules = get_modules(&left_flank, &repeat, &right_flank);
//     let model = Hmm::from(&modules).log();
//     let (_likelihood, annotation) = model.log_predict(&seq, &qual);
// 
//     for x in &annotation { print!("{}", x / 10); }
//     println!();
//     for x in &annotation { print!("{}", x % 10); }
//     println!();
//     println!("{}", str::from_utf8(&seq).unwrap());
//     println!("{}", str::from_utf8(&qual).unwrap());
// 
//     let (partition, mod_ids) = model.partition_to_units(&annotation);
// 
//     let exp_split: Vec<&[u8]> = vec![
//         b"CTTCCTCCTCCTCATCGGTGGCGGCGGCGGCGGCGTCAGGCCAGTGCCGCGGCTTTCTCTCCGCG",
//         b"CCTGTGTTCGCCGGGACGCATTCGGGGCGG",
//         b"GCG", b"GCG", b"GCG", b"GCG", b"GCA", b"GCG", b"GCG", b"GCC", b"GCG", b"GCG", b"GCG", b"GGC", b"GGG", b"GCC", b"GCC",
//         b"CCCCGCCT"
//     ];
// 
//     let exp_mod_ids: Vec<_> = vec![
//         usize::MAX, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2
//     ];
// 
//     println!();
//     let x = get_module_nomenclature(&seq, &partition, &mod_ids, 0);
//     println!("{}", x);
//     let x = get_module_nomenclature(&seq, &partition, &mod_ids, 1);
//     println!("{}", x);
//     let x = get_module_nomenclature(&seq, &partition, &mod_ids, 2);
//     println!("{}", x);
//     let x = get_module_nomenclature(&seq, &partition, &mod_ids, 3);
//     println!("{}", x);
//     println!();
// 
//     for (i, p) in partition.into_iter().enumerate() {
//         println!("{}", str::from_utf8(&seq[p]).unwrap());
//         println!("{}", str::from_utf8(exp_split[i]).unwrap());
//         println!("{} {}", mod_ids[i], exp_mod_ids[i])
//     }
// }
