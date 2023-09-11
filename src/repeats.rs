use std::str::FromStr;
use nom::IResult;
use nom::bytes::complete::{tag, take_until};
use nom::character::complete::{alpha1, alphanumeric0, digit1};
use nom::multi::many0;
use nom::sequence::delimited;
use std::fmt;
use std::str;
use std::collections::HashMap;
use ndarray::{Array, Array2};

#[derive(Default, Debug, PartialEq, Clone)]
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

impl fmt::Display for TandemRepeat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // NC_000008.11:g.118366816_118366918TAAAA[13]TAA[1]TAAAA[7]
        write!(f, "{}:g.", self.reference)?;
        write!(f, "{}_{}", self.start+1, self.end)?;
        for i in 0..self.copy_number.len() {
            write!(f, "{}[{}]",
                str::from_utf8(&self.copy_unit[i]).unwrap(),
                self.copy_number[i]
            )?;
        }
        Ok(())
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
    let start: usize = start.parse().unwrap();  // this is safe
    let (input, _) = tag("_")(input)?;
    let (input, end) = digit1(input)?;
    let end: usize = end.parse().unwrap();      // this is safe
    let (input, repeats) = many0(parse_repeat)(input)?;

    let mut copy_unit = Vec::new();
    let mut copy_number = Vec::new();
    for i in 0..repeats.len() {
        copy_unit.push(repeats[i].0.clone());
        copy_number.push(repeats[i].1);
    }

    Ok((input, TandemRepeat {
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
}

pub fn ref_region<'a>(
    refseq: &'a HashMap<String, Vec<u8>>, id: &str, start: usize, end: usize
) -> Option<&'a[u8]> {
    let seq = match refseq.get(id) {
        None => { return None; },
        Some(x) => { x },
    };
    return Some(&seq[start..end]);
}

pub fn is_present(tr: &TandemRepeat, seq: &HashMap<String, Vec<u8>>) -> bool {
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

const FLANK_SIZE: usize = 4;

pub fn modify_repeat(repeat: &TandemRepeat, refs: &HashMap<String, Vec<u8>>)
    -> TandemRepeat
{
    let repeat_seq = repeat.sequence();
    let refs_seq = ref_region(
        refs, &repeat.reference, repeat.start-FLANK_SIZE, repeat.end+FLANK_SIZE
    ).unwrap().to_owned();

    let dp = fill_dp_table(&repeat_seq, &refs_seq);
    // find end
    // find start
    // report number of edits

    return repeat.clone();
}

use ndarray::array;

fn fill_dp_table(ref_seq: &[u8], mot_seq: &[u8]) -> Array2<u8> {
    // TODO:
    // let mut dp = Array::zeros((ref_seq.len(), mot_seq.len()));
    // println!("{:?}", dp.strides());
    // println!("{:?}", dp.shape());
    return array![[1,2,3], [4,5,6]]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fromstr_and_display_traits_are_inverse() {
        let inputs = [
            "s1:g.10_14A[5]",
            "s2:g.1_16AG[8]",
            "NC_000008.11:g.118366816_118366918TAAAA[13]TAA[1]TAAAA[7]",
            "NM_000044.3:g.123_191CAG[25]"
        ];

        for i in 0..inputs.len() {
            let repeat = inputs[i].parse::<TandemRepeat>().unwrap();
            let repr = format!("{}", repeat);
            assert_eq!(inputs[i], repr);
        }
    }

    #[test]
    fn test_dp_fill() {
        let ref_seq = b"ACCCA";
        let mot_seq = b"CCC";
        let dp = array![
            [0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 1],
            [2, 2, 1, 0, 0, 1],
            [3, 3, 2, 1, 0, 1]
        ];

        let dp2 = fill_dp_table(&ref_seq[..], &mot_seq[..]);
        assert_eq!(dp, dp2);
    }

    #[test]
    fn modify_repeat_can_move_repeat() {
        let motif = "s1:g.10_15A[5]";
        let fasta = [
            ("s1", "CCCCCCCAAAAACCCCCCC")
        ];

        let repeat = motif.parse().unwrap();
        let refs: HashMap<String, Vec<u8>> = HashMap::from_iter(
            fasta.iter().map(
                |(id, seq)| ((*id).to_owned(), (*seq).as_bytes().to_owned())
            )
        );

        let mod_rep = modify_repeat(&repeat, &refs);
        println!("{}", mod_rep);
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
}
