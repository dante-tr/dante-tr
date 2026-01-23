# remaSTR

## Description
Dante is an algorithm designed for genotyping STR alleles based on NGS reads originating from the STR locus of interest.
This method takes into account natural deviations from the expected sequence, including variations in repeat count,
sequencing errors, ambiguous bases, and complex loci containing various motifs.

The reported figures provide evidence for expanded alleles which are too long to be captured by a single NGS read,
as well as allelic single point mutations, small insertions, and deletions that may be relevant for diagnostic evaluations.

## Installation

1. Install [rust](https://rust-lang.org/tools/install/)
2. Install [uv](https://docs.astral.sh/uv/getting-started/installation/)
3. Download, compile, and install python libraries:
```
curl -L https://gitlab.com/andrejbalaz/remastr/-/archive/v0.12.0/remastr-v0.12.0.tar.gz | tar xz
cd remastr-v0.12.0/dante_cli/
cargo build --release
cd ../dante_py/
uv sync
source .venv/bin/activate  # or .venv/bin/activate.fish if you are using fish 
cd ../
```
4. Check that programs can be run
```
./target/release/dante_cli -h
python ./dante_py/dante_remastr_simple.py -h
```

The rust binary is statically linked, therefore it does not need any dynamic libraries
and can be moved to any directory without dependency issues.

## Usage
To annotate reads:
```
remastr -f data/chromosomeX.fna -n data/HGVS.txt -b data/mini2.bam -o output.tsv
```

