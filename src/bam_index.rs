use noodles::bam::{self, bai};
use noodles::bam::io::reader;
use noodles::csi::binning_index::{index::reference_sequence::bin::Chunk, Indexer};
use noodles::sam::{self as sam, alignment::Record};
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

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

#[test]
fn test_coordinate_sorted() {
    let bam_path = "data/real/ilr_lib3.bam";
    let mut reader = bam::io::reader::Builder.build_from_path(bam_path).unwrap();
    let header = reader.read_header().unwrap();

    assert!(is_coordinate_sorted(&header));
}

fn bai<P: AsRef<Path>>(bam: P) -> PathBuf {
    let mut bai = bam.as_ref().to_path_buf();
    bai.set_extension("bam.bai");
    return bai;
}

pub fn check_bai<P: AsRef<Path>>(bam: P) {
    let bai_path = bai(&bam);
    if !bai_path.exists() {
        println!("BAM index (.bai) does not exist. Creating...");
        create_bai_file(&bam).expect("Cannot create index file.");
    }
}

#[test]
fn test_bai_construction() {
    let bam_path = "data/real/ilr_lib3.bam";
    let bai_path = bai(bam_path);
    if !bai_path.exists() {
        println!("Create .bai for bam file.");
        create_bai_file(bam_path).expect("Cannot create .bai file.");
    }

    use noodles::bam::io::indexed_reader::Builder;
    let mut reader = Builder::default()
        .build_from_path(bam_path)
        .expect("Unable to read the associated index (.bai).");

    let _header = reader.read_header().unwrap();
}
