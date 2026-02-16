use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io::Write;
use std::iter::zip;
use std::ops::Range;
use std::path::Path;
use std::str;

use polars::prelude::*;

use nom::AsBytes;
use noodles::bam;

use crate::repeats::TandemRepeat;
use crate::hmm::Hmm;

pub fn print_tsv_file(df: &mut DataFrame, p: &Path) -> Result<(), Box<dyn Error>> {
    let file = File::create(p)?;
    CsvWriter::new(file).with_separator(b'\t').finish(df)?;
    return Ok(());
}

pub fn print_dbg_file(df: &DataFrame, p: &Path) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(p)?;

    // use polars::frame::row::Row;
    // let mut row = Row::default();
    for i in 0..df.height() {
        let row = df.get_row(i)?.0;
        let col_names = df.columns();
        for (name, value) in std::iter::zip(col_names, row) {
            writeln!(file, "{}\t{}", name.name(), value)?;
        }
        writeln!(file)?;
    }
    return Ok(());
}

pub fn annotate_reads<T>(in_reads: T, model: Hmm, repeat: &TandemRepeat) -> DataFrame
where
    T: Iterator<Item = bam::Record>,
{
    // TODO: refactor this
    let mut names = Vec::new();
    let mut motifs = Vec::new();
    let mut read_sns = Vec::new();
    let mut read_ids = Vec::new();
    let mut mate_orders = Vec::new();
    let mut reads = Vec::new();
    let mut references = Vec::new();
    let mut moduleses = Vec::new();
    let mut qualities = Vec::new();
    let mut log_likelihoods = Vec::new();
    let mut left_bgs = Vec::new();
    let mut right_bgs = Vec::new();
    let mut n_moduleses = Vec::new();
    let mut module_baseses = Vec::new();
    let mut module_repetitionses = Vec::new();
    let mut module_sequenceses = Vec::new();
    let mut module_nomenclatureses = Vec::new();
    let mut n_deletionses = Vec::new();
    let mut n_insertionses = Vec::new();
    let mut n_mismatcheses = Vec::new();
    let mut mismatches_strs = Vec::new();
    let mut module_classeses = Vec::new();
    for (i, read) in in_reads.enumerate() {
        let seq: Vec<_> = read.sequence().iter().collect();
        let qual: Vec<_> = read.quality_scores().as_ref().iter().map(|x| x + 33).collect();
        let qual_mod = qual.clone();
        let qual_str: Vec<_> = "No quality".bytes().collect();

        let (likelihood, annotation) = model.log_predict(&seq, &qual_mod);
        let (partition, mod_ids) = model.partition_to_units(&annotation);
        let (new_annot, reconstructed_read) = model.realign(&annotation, &seq);
        let reconstructed_reference = model.reconstruct_sequence(&new_annot);
        let mods = model.reconstruct_mod_ids(&new_annot);

        let name: String = match &repeat.name {
            Some(x) => x.to_string(),
            None    => "None".to_string(),
        };
        names.push(name);
        let motif = repeat.to_string();
        motifs.push(motif);
        let read_sn = i;
        read_sns.push(read_sn as u64);
        let read_id = str::from_utf8(read.name().unwrap().as_bytes()).unwrap();
        read_ids.push(read_id.to_string());
        let mate_order = mate_order(&read);
        mate_orders.push(mate_order);
        let read = str::from_utf8(&reconstructed_read).unwrap();
        reads.push(read.to_string());
        let reference = str::from_utf8(&reconstructed_reference).unwrap();
        references.push(reference.to_string());
        let modules = str::from_utf8(&mods).unwrap();
        moduleses.push(modules.to_string());
        let quality = str::from_utf8(&qual_str).unwrap();
        qualities.push(quality.to_string());
        let log_likelihood = likelihood;
        log_likelihoods.push(log_likelihood);
        let mlen = mods.len();
        let mut left_bg = 0;
        while left_bg < mlen && mods[left_bg] == b'-' { left_bg += 1; }
        left_bgs.push(left_bg as u64);
        let mut right_bg = 0;
        while right_bg < mlen && mods[(mlen - 1) - right_bg] == b'-' { right_bg += 1; }
        right_bgs.push(right_bg as u64);
        let mismatches_str = generate_mismatches(&reconstructed_read, &reconstructed_reference);
        mismatches_strs.push(mismatches_str.clone());
        let n_deletions = mismatches_str.bytes().filter(|x| *x == b'D').count();
        n_deletionses.push(n_deletions as u64);
        let n_insertions = mismatches_str.bytes().filter(|x| *x == b'I').count();
        n_insertionses.push(n_insertions as u64);
        let n_mismatches = mismatches_str.bytes().filter(|x| *x == b'M').count();
        n_mismatcheses.push(n_mismatches as u64);
        let n_modules = repeat.copy_number.len() + 2;
        n_moduleses.push(n_modules as u64);
        let mut module_sequences: Vec<String> = Vec::with_capacity(n_modules);
        let mut module_nomenclatures: Vec<String> = Vec::with_capacity(n_modules);
        let mut module_bases: Vec<usize> = Vec::with_capacity(n_modules);
        let mut module_repetitions: Vec<usize> = Vec::with_capacity(n_modules);
        for i in 0..n_modules {
            let ms = get_module_sequences(&seq, &partition, &mod_ids, i);
            let mn = get_module_nomenclature(&seq, &partition, &mod_ids, i);
            let mb = get_module_bases(&mods, i);
            let mr = get_module_repetitions(mb, &repeat.copy_unit, i);
            module_sequences.push(ms);
            module_nomenclatures.push(mn);
            module_bases.push(mb);
            module_repetitions.push(mr);
        }
        let module_classes = get_module_classes(left_bg, &module_bases, right_bg);
        let module_bases = module_bases.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        module_baseses.push(module_bases);
        let module_repetitions = module_repetitions.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        module_repetitionses.push(module_repetitions);
        let module_sequences = module_sequences.join(",");
        module_sequenceses.push(module_sequences);
        let module_nomenclatures = module_nomenclatures.join(",");
        module_nomenclatureses.push(module_nomenclatures);
        let module_classes = module_classes.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        module_classeses.push(module_classes);
    }

    let result = df![
        "name"                 => names,
        "motif"                => motifs,
        "read_sn"              => read_sns,
        "read_id"              => read_ids,
        "mate_order"           => mate_orders,
        "quality"              => qualities,
        "log_likelihood"       => log_likelihoods,
        "read"                 => reads,
        "reference"            => references,
        "n_modules"            => n_moduleses,
        "left_bg"              => left_bgs,
        "module_bases"         => module_baseses,
        "right_bg"             => right_bgs,
        "module_repetitions"   => module_repetitionses,
        "module_sequences"     => module_sequenceses,
        "module_nomenclatures" => module_nomenclatureses,
        "modules"              => moduleses,
        "n_deletions"          => n_deletionses,
        "n_insertions"         => n_insertionses,
        "n_mismatches"         => n_mismatcheses,
        "mismatches_str"       => mismatches_strs,
        "module_classes"       => module_classeses,
    ].expect("Cannot create dataframe");
    return result;
}

#[derive(Clone)]
enum AClass {
    Spanning,
    Flanking,
    Missing,
    Filtered(String)
}

impl fmt::Display for AClass {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AClass::Missing => write!(f, "Missing"),
            AClass::Spanning => write!(f, "Spanning"),
            AClass::Flanking => write!(f, "Flanking"),
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
    return result;
}

fn get_module_nomenclature(seq: &[u8], partition: &[Range<usize>], mod_ids: &[usize], idx: usize) -> String {
    let mut mn = Vec::new();
    let mut append_unit = |s, o|{
        if o != 0 {
            mn.push(format!("{}[{}]", str::from_utf8(s).unwrap(), o));
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
    return mn.join("");
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

fn get_module_repetitions(mb: usize, copy_units: &[Vec<u8>], idx: usize) -> usize {
    // TODO: this should really reflect how many times the HMM passed through module.
    // Now it should be easy to implement, but for parity with python it is implemented like this.
    // But definitely, change it in the future.
    if mb == 0 { return 0; }
    if idx == 0 { return 1; }
    if idx == copy_units.len() + 1 { return 1; }
    if idx > copy_units.len() + 1 { panic!("This should never happen."); }
    let copy_len: u8 = copy_units[idx - 1].len().try_into().unwrap();
    let mb: f64 = mb as f64;  // f64 cannot represent all usize values. Potentially dangerous
    let copy_len: f64 = copy_len.into();
    let res = mb / copy_len;
    // fNN as iNN is defined to be a truncating cast, saturating out-of-range values and mapping NaN to 0.
    // If that's what you want, just writing as is currently the best way to get that behavior.
    // https://users.rust-lang.org/t/floor-and-cast-f64-to-usize-in-one-operation/88768/4
    return res.round_ties_even() as usize;
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

#[test]
fn tmp_fn_name() {
    use crate::get_modules;
    // CTTCCTCCTCCTCATCGGTGGCGGCGGCGGCGGCGTCAGGCCAGTGCCGCGGCTTTCTCTCCGCGCCTGTGTTCGCCGGGACGCATTCGGGGCGGGCGGCGGCGGCGGCAGCGGCGGCCGCGGCGGCGGGCGGGGCCGCCCCCCGCCT
    // CCCFFFFFHHHHGIIJJIHIJIJJJHGHDDDBDDBBDDDDDDDDCCCC5@BDDBBCCCCACDBBDBDBBCCCCCBDDBBD@BDDBDDDDDDDDDD@5.5&)0)&5))0)&)2&55)0&&&&&)&&)&))&&&&&)&&&&)&&50)&&&
    // -----------------------------------------------------------------00000000000000000000000000000011111111111111111111111111111111111111111111122222222
    // GCG[4]GCA[1]GCG[2]GCC[1]GCG[3]G[1]GCG[1]GGGCCGCC[1]
    //
    // I think annotation is wrong
    // Independent of annotation, seq nomenclature is wrong as well
    //
    // SPD     chr2:g.176093059_176093103GCG[15]       CCTGTGTTCGCCGGGACGCATTCGGGGCGG  TCCGGCTTTGCGTACCCCGGGACCTCTGAG
    // result 15

    // HISEQ1:26:HA2RRADXX:1:1203:16720:7919
    let seq:  Vec<u8> = b"CTTCCTCCTCCTCATCGGTGGCGGCGGCGGCGGCGTCAGGCCAGTGCCGCGGCTTTCTCTCCGCGCCTGTGTTCGCCGGGACGCATTCGGGGCGGGCGGCGGCGGCGGCAGCGGCGGCCGCGGCGGCGGGCGGGGCCGCCCCCCGCCT".to_vec();
    let qual: Vec<u8> = b"CCCFFFFFHHHHGIIJJIHIJIJJJHGHDDDBDDBBDDDDDDDDCCCC5@BDDBBCCCCACDBBDBDBBCCCCCBDDBBD@BDDBDDDDDDDDDD@5.5&)0)&5))0)&)2&55)0&&&&&)&&)&))&&&&&)&&&&)&&50)&&&".to_vec();

    // SPD
    let left_flank:  Vec<u8> = b"CCTGTGTTCGCCGGGACGCATTCGGGGCGG".to_vec();
    let right_flank: Vec<u8> = b"TCCGGCTTTGCGTACCCCGGGACCTCTGAG".to_vec();
    let repeat: TandemRepeat = "chr2:g.176093059_176093103GCG[15]".parse().expect("Malformatted nomenclature found.");

    let modules = get_modules(&left_flank, &repeat, &right_flank);
    let model = Hmm::from(&modules).log();
    let (_likelihood, annotation) = model.log_predict(&seq, &qual);

    for x in &annotation { print!("{}", x / 10); }
    println!();
    for x in &annotation { print!("{}", x % 10); }
    println!();
    println!("{}", str::from_utf8(&seq).unwrap());
    println!("{}", str::from_utf8(&qual).unwrap());

    let (partition, mod_ids) = model.partition_to_units(&annotation);

    let exp_split: Vec<&[u8]> = vec![
        b"CTTCCTCCTCCTCATCGGTGGCGGCGGCGGCGGCGTCAGGCCAGTGCCGCGGCTTTCTCTCCGCG",
        b"CCTGTGTTCGCCGGGACGCATTCGGGGCGG",
        b"GCG", b"GCG", b"GCG", b"GCG", b"GCA", b"GCG", b"GCG", b"GCC", b"GCG", b"GCG", b"GCG", b"GGC", b"GGG", b"GCC", b"GCC",
        b"CCCCGCCT"
    ];

    let exp_mod_ids: Vec<_> = vec![
        usize::MAX, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2
    ];

    println!();
    let x = get_module_nomenclature(&seq, &partition, &mod_ids, 0);
    println!("{}", x);
    let x = get_module_nomenclature(&seq, &partition, &mod_ids, 1);
    println!("{}", x);
    let x = get_module_nomenclature(&seq, &partition, &mod_ids, 2);
    println!("{}", x);
    let x = get_module_nomenclature(&seq, &partition, &mod_ids, 3);
    println!("{}", x);
    println!();

    for (i, p) in partition.into_iter().enumerate() {
        println!("{}", str::from_utf8(&seq[p]).unwrap());
        println!("{}", str::from_utf8(exp_split[i]).unwrap());
        println!("{} {}", mod_ids[i], exp_mod_ids[i])
    }
}
