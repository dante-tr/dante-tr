use std::fs::File;
use std::path::Path;
use std::path::PathBuf;

use noodles::bam;
use noodles::bgzf::Reader;
use noodles::core::Region;
use noodles::sam::Header;

use noodles::bam::io::Writer;
use noodles::bgzf as bgzf;

use noodles::bam::bai;
use noodles::bam::io::reader;
use noodles::csi::binning_index::{index::reference_sequence::bin::Chunk, Indexer};
use noodles::sam::{self as sam, alignment::Record};
use std::io;

pub(crate) struct RelevantReads {
    reader: bam::io::indexed_reader::IndexedReader<Reader<File>>,
    header: Header,
    region: Region,
}

impl RelevantReads {
    pub(crate) fn from(bam_file: &Path, region_str: &str) -> RelevantReads {
        let region = region_str.parse().unwrap();
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_file)
            .expect("Unable to read the associated index (.bai).");
        let header = reader.read_header().expect("Error. TODO");
        RelevantReads { reader, header, region }
    }

    pub(crate) fn header(&self) -> Header { self.header.clone() }

    pub(crate) fn iter(&mut self) -> impl Iterator<Item = bam::Record> + '_ {
        self.reader
            .query(&self.header, &self.region).expect("Error. TODO")
            .map(|x| x.expect("Error. TODO"))
    }

    pub(crate) fn write_to_file(&mut self, out_bam_file: &Path) {
        let h = self.header();
        let mut out_bam = init_bam(&out_bam_file.to_string_lossy(), &h);
        for record in self.iter() {
            out_bam.write_record(&h, &record).expect("Cannot write to out bam.");
        }
        // TODO: sort bam + create bai index
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

pub(crate) fn check_bai<P: AsRef<Path>>(bam: P) {
    let bai_path = bai(&bam);
    if !bai_path.exists() {
        println!("BAM index (.bai) does not exist. Creating...");
        create_bai_file(&bam).expect("Cannot create index file.");
    }
}

fn create_bai_file<P: AsRef<Path>>(bam_file: P) -> io::Result<()> {
    let mut reader = reader::Builder.build_from_path(&bam_file)?;
    let header = reader.read_header()?;

    if !is_coordinate_sorted(&header) {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "BAM must be coordinate-sorted"));
    }

    let mut builder = Indexer::default();
    let mut start_position = reader.get_ref().virtual_position();

    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let alignment_context = match (
            record.reference_sequence_id().transpose()?,
            record.alignment_start().transpose()?,
            record.alignment_end().transpose()?,
        ) {
            (Some(id), Some(start), Some(end)) => {
                let is_mapped = !record.flags().is_unmapped();
                Some((id, start, end, is_mapped))
            },
            _ => None,
        };

        builder.add_record(alignment_context, chunk)?;
        start_position = end_position;
    }

    let index = builder.build(header.reference_sequences().len());

    let bai_path = bai(&bam_file);
    let mut writer = File::create(bai_path).map(bai::Writer::new)?;
    writer.write_index(&index)?;

    Ok(())
}

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    use sam::header::record::value::map::header::{sort_order, tag};
    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&tag::SORT_ORDER))
        .map(|sort_order| sort_order == sort_order::COORDINATE)
        .unwrap_or_default()
}

fn bai<P: AsRef<Path>>(bam: P) -> PathBuf {
    let mut bai = bam.as_ref().to_path_buf();
    bai.set_extension("bam.bai");
    return bai;
}
