use std::str::FromStr;
use nom::IResult;
use nom::bytes::complete::{tag, take_until};
use nom::character::complete::{alpha1, alphanumeric0, digit1};
use nom::multi::many0;
use nom::sequence::delimited;

#[derive(Default, Debug, PartialEq)]
pub struct TandemRepeat {
    pub reference: String,
    pub start: usize,
    pub end: usize,
    pub copy_unit: Vec<Vec<u8>>,
    pub copy_number: Vec<usize>
}

#[derive(Debug, PartialEq, Eq)]
pub struct ParseTandemRepeatError;

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
    use super::*;

    #[test]
    fn repeat_can_be_parsed() {
        let s = "NM_000044.3:g.123_191CAG[25]";
        let tr: TandemRepeat = s.parse().unwrap();
        println!("{:?}", tr);
    }

    #[test]
    fn complex_repeat_can_be_parsed() {
        let s = "NC_000008.11:g.118366816_118366918TAAAA[13]TAA[1]TAAAA[7]";
        let tr: TandemRepeat = s.parse().unwrap();
        println!("{:?}", tr);
    }
}
