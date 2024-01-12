use noodles::bam::{self as bam, bai};
use noodles::csi::binning_index::{index::reference_sequence::bin::Chunk, Indexer};
use noodles::sam::{self as sam, alignment::Record};
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    use sam::header::record::value::map::header::SortOrder;

    if let Some(hd) = header.header() {
        if let Some(sort_order) = hd.sort_order() {
            return sort_order == SortOrder::Coordinate;
        }
    }
    return false;
}

fn create_bai_file<P: AsRef<Path>>(bam_file: P) -> io::Result<()> {
    let mut reader = bam::reader::Builder.build_from_path(&bam_file)?;
    let header = reader.read_header()?;

    if !is_coordinate_sorted(&header) {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "BAM must be coordinate-sorted"));
    }

    let mut builder = Indexer::default();
    let mut start_position = reader.virtual_position();

    let mut record = Record::default();
    while reader.read_record(&header, &mut record)? != 0 {
        let end_position = reader.virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let alignment_context =
            match (record.reference_sequence_id(), record.alignment_start(), record.alignment_end()) {
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

fn bai<P: AsRef<Path>>(bam: P) -> PathBuf {
    let mut bai = bam.as_ref().to_path_buf();
    bai.set_extension("bam.bai");
    return bai;
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_bai_construction() {
        let bam_path = "data/real/ilr_lib3.bam";
        let bai_path = bai(bam_path);
        if !bai_path.exists() {
            println!("Create .bai for bam file.");
            create_bai_file(bam_path).expect("Cannot create .bai file.");
        }

        let mut reader = bam::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .expect("Unable to read the associated index (.bai).");

        let _header = reader.read_header().unwrap();
    }
}
