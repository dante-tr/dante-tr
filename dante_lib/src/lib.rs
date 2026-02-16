mod annotation;
mod bam_index;
mod hmm;
mod io;
mod motif_correction;
mod repeats;

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
use std::path::PathBuf;
use std::str;

use rayon::prelude::*;
use noodles::bam::io::Writer;
use noodles::bam;
use noodles::bgzf as bgzf;
use noodles::sam::Header;

use crate::annotation::{annotate_reads, print_dbg_file, print_tsv_file};
use crate::bam_index::check_bai;
use crate::bam_ops::RelevantReads;
use crate::hmm::{Module, Hmm};
use crate::io::get_modules;
use crate::repeats::TandemRepeat;

pub fn run_v2(bam_file: &Path, motif_file: &Path, output: &Path, out_bam_flag: bool) {
    check_bai(bam_file);

    let motif_records = read_motifs(motif_file);
    motif_records.par_iter().for_each(|motif_record| {

        let (left_flank, repeat, right_flank) = motif_record;
        let name = repeat.name.as_ref().unwrap().clone();
        let region_str = format!("{}:{}-{}", repeat.reference, repeat.start + 1, repeat.end);
        let region = region_str.parse().unwrap();

        let mut relevant_reads = RelevantReads::from(bam_file, region);
        if out_bam_flag {
            let h = relevant_reads.header();
            let out_bam_file = output.join(name.to_owned() + ".annotated.bam");
            let mut out_bam = init_bam(&out_bam_file.to_string_lossy(), &h);
            for record in relevant_reads.iter() {
                out_bam.write_record(&h, &record).expect("Cannot write to out bam.");
            }
            // TODO: sort bam + create bai index
        }

        // build HMM and annotate reads - polars alternative
        let model = Hmm::from(&get_modules(left_flank, repeat, right_flank)).log();
        let mut annotation_df = annotate_reads(relevant_reads.iter(), model, repeat);

        // write results to tsv
        let out_tsv_file = output.join(name.to_owned() + ".annotations.tsv");
        print_tsv_file(&mut annotation_df, &out_tsv_file).expect("Failed writing tsv file.");

        let out_tsv_file = output.join(name.to_owned() + ".annotations.dbg.txt");
        print_dbg_file(&annotation_df, &out_tsv_file).expect("Failed writing dbg file.");
    });

    println!("Annotation finished successfully.");
}

fn read_motifs(filename: &Path) -> Vec<(Vec<u8>, TandemRepeat, Vec<u8>)> {
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
 
fn init_bam(tsv_file: &str, header: &Header) -> Writer<bgzf::Writer<File>> {
    let mut filename = PathBuf::from(tsv_file);
    filename.set_extension("bam");
    let new_bam = File::create(filename).expect("Cannot open file for writing.");
    let mut writer = bam::io::Writer::new(new_bam);
    writer.write_header(header).unwrap();
    return writer;
}

mod bam_ops {
    use std::fs::File;
    use std::path::Path;

    use noodles::bam;
    use noodles::bgzf::Reader;
    use noodles::core::Region;
    use noodles::sam::Header;

    pub struct RelevantReads {
        reader: bam::io::indexed_reader::IndexedReader<Reader<File>>,
        header: Header,
        region: Region,
    }

    impl RelevantReads {
        pub fn from(bam_file: &Path, region: Region) -> RelevantReads {
            let mut reader = bam::io::indexed_reader::Builder::default()
                .build_from_path(bam_file)
                .expect("Unable to read the associated index (.bai).");
            let header = reader.read_header().expect("Error. TODO");
            RelevantReads { reader, header, region }
        }

        pub fn header(&self) -> Header { self.header.clone() }

        pub fn iter(&mut self) -> impl Iterator<Item = bam::Record> + '_ {
            self.reader
                .query(&self.header, &self.region).expect("Error. TODO")
                .map(|x| x.expect("Error. TODO"))
        }
    }
}
