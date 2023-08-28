use noodles::bam as bam;
use noodles::fasta as fasta;
use noodles::sam::record::quality_scores::Score;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::str;
use clap::Parser;

use crate::hmm::HMM;
use crate::hmm::Module;
use crate::repeats::TandemRepeat;
use crate::cli::Args;

mod hmm;
mod repeats;
mod cli;

fn main() {
    let args = Args::parse();
    let references = read_reference(&args.ref_file);
    let repeats = read_nomenclature(&args.hgvs_file);
    let bam_refs = read_bam_refs(&args.bam_file);

    let mut valid_repeats = Vec::new();
    for repeat in repeats {
        if is_present(&repeat, &references) {
            valid_repeats.push(repeat);
        }
    }

    valid_repeats.par_iter().for_each(|repeat| {
        let mut reader = bam::indexed_reader::Builder::default()
            .build_from_path(&args.bam_file).unwrap();
        let header = reader.read_header().unwrap();

        let modules = get_modules(&repeat, &references, 20);
        let model = HMM::from(&modules).log();

        let region = format!("{}:{}-{}", repeat.reference, repeat.start+1, repeat.end).parse().unwrap();
        let reads = reader.query(&header, &region).unwrap();

        for read in reads {
            let read = read.expect("Incorrect read.");
            let seq: Vec<_> = read.sequence().as_ref().iter().map(|&x| x.into()).collect();
            let qual: Vec<_> = read.quality_scores().as_ref().iter().map(|&x| remap(x)).collect();
            let (likelihood, annotation) = model.log_predict(&seq, &qual);

            let reconstructed_reference = model.reconstruct_sequence(&annotation);
            let reconstructed_read = model.realign_read(&annotation, &seq); 
            let mods = model.reconstruct_mod_ids(&annotation);

            println!(">{} {} {}\n{}\n{}\n{}", 
                read.read_name().unwrap(), repeat, likelihood,
                str::from_utf8(&reconstructed_read).unwrap(),
                str::from_utf8(&reconstructed_reference).unwrap(),
                str::from_utf8(&mods).unwrap()
            );
        }
    })
}

fn remap(x: Score) -> u8 {
    let c: char = x.into();
    return c as u8;
}

fn read_reference(filename: &str) -> HashMap<String, Vec<u8>> {
    let mut reader = fasta::reader::Builder
        .build_from_path(filename).unwrap();

    let mut result = HashMap::new();
    for record in reader.records() {
        let record = record.unwrap();

        result.insert(
            record.name().to_string(),
            (&record.sequence()[..]).to_vec()
            // ^- Is there a better way to get Vec<u8>
            // Do I need Vec<u8>? Cannot I leave it as Sequence?
        );
    }
    return result;
}

fn read_nomenclature(filename: &str) -> Vec<TandemRepeat> {
    let mut repeats = Vec::new();
    
    let file = File::open(filename).expect("Cannot find nomenclature file.");
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line
            .expect("Cannot read line from nomenclature file.")
            .trim().to_owned();
        let repeat = line.parse().expect("Cannot parse nomenclature.");
        repeats.push(repeat);
    }
    return repeats;
}

fn read_bam_refs(filename: &str) -> HashMap<String, usize> {
    let mut result = HashMap::new();
    
    let file = File::open(filename).unwrap();
    let header = bam::Reader::new(file).read_header().unwrap();

    for s in header.reference_sequences().iter() {
        let name = s.0.to_string();
        let length = s.1.length().get();
        result.insert(name.clone(), length);
    }
    return result;
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

fn get_modules(
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
    for i in 0..repeat.copy_unit.len() {
        modules.push((&repeat.copy_unit[i][..], repeat.copy_number[i]).into());
    }
    modules.push(right_flank.into());
    return modules;
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

fn bam_wrt_ref_is_valid(bam_refs: &HashMap<String, usize>, references: &HashMap<String, Vec<u8>>) -> bool {
    for (id, &len) in bam_refs {
        let seq = match references.get(id) {
            None => {
                println!("{} not found in reference.", id);
                return false;
            },
            Some(s) => { s }
        };
        if len != seq.len() {
            println!("{} lengths differ in bam and fasta.", id);
            return false;
        }
    } 
    return true;
}

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
    // for item in ref_map {
    //     let new_id = "".to_string();
    //     let seq = "N".repeat(1000);
    //     result.insert(new_id, seq);
    // }
    return result;
}

fn correct_ref(references: HashMap<String, Vec<u8>>) -> HashMap<String, Vec<u8>> {
    let mut result = HashMap::new();
    // chr1    248956422
    // chr2    242193529
    // chr3    198295559
    // chr4    190214555
    // chr5    181538259
    // chr6    170805979
    // chr7    159345973
    // chr8    145138636
    // chr9    138394717
    // chr10   133797422
    // chr11   135086622
    // chr12   133275309
    // chr13   114364328
    // chr14   107043718
    // chr15   101991189
    // chr16   90338345
    // chr17   83257441
    // chr18   80373285
    // chr19   58617616
    // chr20   64444167
    // chr21   46709983
    // chr22   50818468
    // chrX    156040895
    // chrY    57227415
    // chrM    16569

//     for (id, seq) in references {
//         let new_id = "".to_string();
//         result.insert(new_id, seq);
//     }
//     for item in ref_map {
//         let new_id = "".to_string();
//         let seq = "N".repeat(1000);
//         result.insert(new_id, seq);
//     }
    return result;
}

fn correct_repeats(repeats: Vec<TandemRepeat>) -> Vec<TandemRepeat> {
    return repeats;
}

#[cfg(test)]
mod tests {
    use hgvs::parser::HgvsVariant;
    use std::fs::File;
    use std::io::{prelude::*, BufReader};
    use super::*;

    #[test]
    fn can_get_sequence_id_from_bam() {
        let filename: &str = "data/mini2.bam";
        let file = File::open(filename).unwrap();
        let index = bam::bai::read(filename.to_owned() + ".bai").unwrap();
        let mut reader = bam::IndexedReader::new(file.try_clone().unwrap(), index.clone());

        let header = reader.read_header().unwrap();
        let seqs = header.reference_sequences();
        for s in seqs.iter() {
            let name = s.0.to_string();
            let length = s.1.length().get();
            println!("{} {}", name, length);
        }
    }

    #[test]
    fn can_get_sequence_id_from_fasta() {
        let filename = "data/chromosomeX.fna";
        let file = File::open(filename).unwrap();
        let mut reader = fasta::Reader::new(BufReader::new(file));

        for record in reader.records() {
            let record = record.unwrap();
            let name = record.name().to_string();
            let length = record.sequence().len();
            println!("{} {}", name, length);
        }
    }

    #[test]
    fn test_reference_checking() {
        let bam_refs = read_bam_refs("/home/balaz/projects/STRs/remaSTR/data/real/twist_S22-157-01_S1.bam");
        let mut references = read_reference("data/chromosomeX.fna");
        let mut repeats = read_nomenclature("data/mini_HGVS.txt");
        for (k, v) in &bam_refs { println!("{}\t{}", k, v); }
        println!("{:?}", references.keys());
        println!("{:?}", repeats);

        if ! bam_wrt_ref_is_valid(&bam_refs, &references) {
            println!("IDs in bam and reference differ. Attempting correction of reference... ");
            references = correct_ref(references);
            match bam_wrt_ref_is_valid(&bam_refs, &references) {
                true => { println!("Success!"); }
                false => { panic!("Unable to correct reference!"); }
            }
        }
        if ! hgvs_wrt_ref_is_valid(&repeats, &references) {
            println!("IDs in hgvs and reference differ. Attempting correction of hgvs... ");
            repeats = correct_repeats(repeats);
            match hgvs_wrt_ref_is_valid(&repeats, &references) {
                true => { println!("Success!"); },
                false => { panic!("Unable to correct repeats!"); }
            }
        }
    }

    #[test]
    fn input_can_be_remapped() {
        let args = Args::try_parse_from([
            "remastr",
            "-f", "data/chromosomeX.fna",
            "-b", "data/real/twist_S22-157-01_S1.bam",
            "-n", "data/chromosomeX.fna",
            "--map1", "data/chromosomeX_to_bam_map.txt"
        ]).unwrap();

        let bam_refs = read_bam_refs(&args.bam_file);
        let mut references = read_reference(&args.ref_file);

        if ! bam_wrt_ref_is_valid(&bam_refs, &references) {
            println!("IDs in BAM and reference differ. Attempting correction of reference... ");
            if let Some(ref2bam) = args.ref2bam {
                println!("Correcting with map provided in {}.", ref2bam);
                references = fix_reference(references, &bam_refs, &ref2bam);
            } else {
                println!("Correcting with best effort heuristic.");
            }

            match bam_wrt_ref_is_valid(&bam_refs, &references) {
                true => { println!("Success!"); }
                false => { panic!("Unable to correct reference!"); }
            }
        }

    }

    #[test]
    fn can_load_fasta() {
        let sequences = read_reference("data/chromosomeX.fna");
        let hgvs = File::open("data/mini_HGVS.txt").unwrap();
        let reader = BufReader::new(hgvs);

        let expected = vec![
            false, false, true, false, false, false, false, true, true, false
        ];
        for (i, line) in reader.lines().enumerate() {
            let line = line.unwrap();
            let line = line.trim();
            let tr: TandemRepeat = line.parse().unwrap();
            let is_correct = is_present(&tr, &sequences);
            assert_eq!(is_correct, expected[i]);
        }
    }

    #[test]
    fn count_present() {
        let references = read_reference("data/chromosomeX.fna");
        let hgvs = File::open("data/HGVS.txt").unwrap();
        let reader = BufReader::new(hgvs);

        let mut present_count = 0;
        let mut max_count = 0;
        for line in reader.lines() {
            let line = line.unwrap().trim().to_owned();
            let tr: TandemRepeat = line.parse().unwrap();
            if is_present(&tr, &references) {
                present_count += 1;
            } else {
                println!("{}", tr);
                print_diff(&tr, &references);
                // println!();
            }
            max_count += 1;
        }
        println!("Present repeats: {}/{}", present_count, max_count);
    }

    fn print_diff(tr: &TandemRepeat, refs: &HashMap<String, Vec<u8>>) {
        let n = 10;
        let rflank = ref_region(refs, &tr.reference, tr.start-n, tr.start).unwrap();
        let ref_repeat = ref_region(refs, &tr.reference, tr.start, tr.end).unwrap();
        let lflank = ref_region(refs, &tr.reference, tr.end, tr.end+n).unwrap();
        println!("{} {} {}", 
            str::from_utf8(rflank).unwrap(),
            str::from_utf8(ref_repeat).unwrap(),
            str::from_utf8(lflank).unwrap()
        );
        println!("{} {} {}",
            " ".repeat(n),
            str::from_utf8(&tr.sequence()).unwrap(),
            " ".repeat(n)
        );
    }

    #[test]
    fn can_parse_hgvs() {
        let _record = "NM_01234.5:c.456-6_*22A>T";
        let _record = "NC_000017.11:g.43091687del";
        let tmp: HgvsVariant = _record.parse().unwrap();
        println!("{:?}", tmp);

        println!("{}", tmp.accession().value);
    }

    #[test]
    fn can_read_and_parse_hgvs_file() {
        let file = "data/mini_HGVS.txt";
        let file = File::open(file).unwrap();
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line.unwrap();
            let line = line.trim();
            println!("{}", line);
            let tr: TandemRepeat = line.parse().unwrap();
            println!("{:?}", tr);
        }
    }
}

