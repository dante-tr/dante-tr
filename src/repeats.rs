use nom::bytes::complete::{tag, take_until};
use nom::character::complete::{alpha1, alphanumeric0, digit1};
use nom::multi::many0;
use nom::sequence::delimited;
use nom::IResult;
use std::fmt;
use std::str;
use std::str::FromStr;

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
    fn test_inner_motif_view() {
        let motif: TandemRepeat = "S1:g.1_10A[10]".parse().unwrap();
        let v = motif.view(0, 10);
        println!("{}", str::from_utf8(&v).unwrap());
    }

    #[test]
    fn complex_repeat_can_be_parsed() {
        let s = "NC_000008.11:g.118366816_118366918TAAAA[13]TAA[1]TAAAA[7]";
        let tr: TandemRepeat = s.parse().unwrap();
        println!("{:?}", tr);
    }
}

#[derive(Default, Debug, PartialEq, Clone)]
pub struct TandemRepeat {
    pub name: Option<String>,
    pub reference: String,
    pub start: usize,
    pub end: usize,
    pub copy_unit: Vec<Vec<u8>>,
    pub copy_number: Vec<usize>,
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

impl fmt::Display for TandemRepeat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:g.", self.reference)?;
        write!(f, "{}_{}", self.start + 1, self.end)?;
        for i in 0..self.copy_number.len() {
            write!(f, "{}[{}]", str::from_utf8(&self.copy_unit[i]).unwrap(), self.copy_number[i])?;
        }
        Ok(())
    }
}

#[test]
fn fmt_is_inverse_to_parse() {
    let motif1 = "NC_000008.11:g.118366816_118366918TAAAA[20]";
    let motif2 = format!("{}", motif1.parse::<TandemRepeat>().unwrap());
    assert_eq!(motif1, motif2);
}

fn parse_repeat(input: &str) -> IResult<&str, (Vec<u8>, usize)> {
    let (remaining, unit) = alpha1(input)?;
    let unit = unit.as_bytes().to_vec();
    let (remaining, number) = delimited(tag("["), digit1, tag("]"))(remaining)?;
    let number = number.parse().unwrap(); // this is safe due to previous line
    return Ok((remaining, (unit, number)));
}

fn tandem_repeat(input: &str) -> IResult<&str, TandemRepeat> {
    let (input, reference) = take_until(":")(input)?;
    let (input, _) = delimited(tag(":"), alphanumeric0, tag("."))(input)?;
    let (input, start) = digit1(input)?;
    let start: usize = start.parse().unwrap(); // this is safe
    let (input, _) = tag("_")(input)?;
    let (input, end) = digit1(input)?;
    let end: usize = end.parse().unwrap(); // this is safe
    let (input, repeats) = many0(parse_repeat)(input)?;

    let mut copy_unit = Vec::new();
    let mut copy_number = Vec::new();
    for r in &repeats {
        copy_unit.push(r.0.clone());
        copy_number.push(r.1);
    }

    Ok((input, TandemRepeat {
        name: None,
        reference: reference.to_string(),
                            // HGVS is 1-based
        start: start - 1,   // 1-based -> 0-based
        end,                // 1-based -> 0-based+1 (half-open) is nop
        copy_unit, copy_number
    }))
}

impl TandemRepeat {
    pub fn sequence(&self) -> Vec<u8> {
        let mut res = Vec::new();
        for i in 0..self.copy_number.len() {
            for _ in 0..self.copy_number[i] {
                res.extend_from_slice(&self.copy_unit[i]);
            }
        }
        return res;
    }

    pub fn correct_boundary(&self) -> Self {
        let mut new = self.clone();
        let l = new.sequence().len();
        if new.start + l != new.end {
            eprintln!("{} has incorrect end.", new);
            new.end = new.start + l;
            eprintln!("Corrected to {}", new);
        }
        return new;
    }

    #[cfg(test)]
    pub fn view(&self, from: usize, to: usize) -> Vec<u8> {
        if from <= self.start && self.end <= to {
            let mut v = b"-".repeat(to - from);
            let seq = self.sequence();
            for i in 0..seq.len() { v[self.start-from+i] = seq[i]; }
            return v;
        } else {
            eprintln!("View not yet implemented.");
            return Vec::new();
        }
    }
}

