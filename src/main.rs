use std::str::FromStr;
use nom::{
    IResult,
    bytes::complete::tag,
    character::complete::alphanumeric0,
};
use nom::bytes::complete::take_until;
use nom::sequence::delimited;
use nom::character::complete::digit1;
use nom::multi::many0;
use nom::character::complete::alpha1;

fn main() {
    // read nomenclature
    // read reference
    // check nomenclature w.r.t. reference
    // load bam
    // for all nomenclatures:
    //     build HMM
    //     reads = bam.query()
    //     for read in reads:
    //         prob, annotation = HMM.annotate(read)
    //         postfilter
    //         report()
    //     report_row()
    println!("Bu!");
}

#[derive(Default, Debug, PartialEq)]
struct TandemRepeat {
    reference: String,
    start: usize,
    end: usize,
    copy_unit: Vec<Vec<u8>>,
    copy_number: Vec<usize>
}

#[derive(Debug, PartialEq, Eq)]
struct ParseTandemRepeatError;

impl FromStr for TandemRepeat {
    type Err = ParseTandemRepeatError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        let (input, tr) = tandem_repeat(input).map_err(|_| ParseTandemRepeatError)?;
        if !input.is_empty() { return Err(ParseTandemRepeatError); }
        return Ok(tr);
    }
}

fn parse_repeat(input: &str) -> IResult<&str, (Vec<u8>, usize)> {
    let (remaining, unit) = alpha1(input)?;
    let unit = unit.as_bytes().to_vec();
    let (remaining, number) = delimited(tag("["), digit1, tag("]"))(remaining)?;
    let number = number.parse().unwrap();   // this is safe due to previous line
    return Ok((remaining, (unit, number)));
}

fn tandem_repeat(input: &str) -> IResult<&str, TandemRepeat> {
    let (input, reference) = take_until(":")(input)?;
    let (input, _) = delimited(tag(":"), alphanumeric0, tag("."))(input)?;
    let (input, start) = digit1(input)?;
    let (input, _) = tag("_")(input)?;
    let (input, end) = digit1(input)?;
    let (input, repeats) = many0(parse_repeat)(input)?;

    let mut copy_unit = Vec::new();
    let mut copy_number = Vec::new();
    for i in 0..repeats.len() {
        copy_unit.push(repeats[i].0.clone());
        copy_number.push(repeats[i].1);
    }

    Ok((input, TandemRepeat {
        reference: reference.to_string(),
        start: start.parse().unwrap(),
        end: end.parse().unwrap(),
        copy_unit, copy_number
    }))
}

#[cfg(test)]
mod tests {
    use noodles::bam as bam;
    use noodles::fasta as fasta;
    use hgvs::parser::HgvsVariant;
    use std::fs::File;
    use std::io::{prelude::*, BufReader};
    use super::TandemRepeat;

    #[test]
    fn can_load_bam() {
        let mut reader = bam::indexed_reader::Builder::default()
            .build_from_path("data/mini.bam").unwrap();

        let header = reader.read_header().unwrap();

        // let region = "sq0:5-8".parse().unwrap();
        // let query = reader.query(&header, &region).unwrap();

        for result in reader.records(&header) {
            let record = result.unwrap();
            println!("{:?}", record);
        }
    }

    #[test]
    fn can_load_fasta() {
        let mut reader = fasta::reader::Builder
            .build_from_path("data/chromosomeX.fna").unwrap();

        for result in reader.records() {
            let record = result.unwrap();

            println!("{}\t{}", record.name(), record.sequence().len());
        }
    }

//    #[test]
//    fn can_parse_hgvs() {
//        let _record = b"NC_000023.11:g.2789717_2789870ATTTT[30]";
//        let _record = "NM_000044.3:g.123_191CAG[25]";
//        // let _record = "NM_01234.5:c.456-6_*22A>T";
//        // let _record = "NC_000017.11:g.43091687del";
//        let tmp: HgvsVariant = _record.parse().unwrap();
//        println!("{:?}", tmp);
//
//        println!("{}", tmp.accession().value);
//        // println!("{}", tmp.loc_edit().loc);
//    }

    #[test]
    fn repeat_can_be_parsed() {
        let s = "NM_000044.3:g.123_191CAG[25]";
        let tr: TandemRepeat = s.parse().unwrap();
        println!("{:?}", tr);
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

