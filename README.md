# dante-tr

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
wget https://github.com/dante-tr/dante-tr/archive/refs/heads/main.zip
unzip main.zip
cd dante-tr-main/

cd dante_cli/
cargo build --release
cd ../

cd dante_py/
uv sync
cd ../
```
4. Check that programs can be run
```
./dante-tr ./example_data/01_in_dante_nomenclatures_predominant.tsv ./example_data/in_HG002.GRCh38.selected_w_pairs.bam ./dante_output
```

