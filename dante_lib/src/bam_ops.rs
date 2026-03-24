use std::fs::File;
use std::path::Path;
use std::path::PathBuf;

use noodles::bam;
use noodles::bgzf::Reader;
use noodles::core::Region;
use noodles::sam::Header;

use noodles::bam::io::Writer;
use noodles::bgzf as bgzf;

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

pub(crate) fn init_bam(tsv_file: &str, header: &Header) -> Writer<bgzf::Writer<File>> {
    let mut filename = PathBuf::from(tsv_file);
    filename.set_extension("bam");
    let new_bam = File::create(filename).expect("Cannot open file for writing.");
    let mut writer = bam::io::Writer::new(new_bam);
    writer.write_header(header).unwrap();
    return writer;
}
