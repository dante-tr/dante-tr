# remaSTR

## Description
Given a set of reads aligned to reference and a set of short tandem repeat (STR) motifs, remaSTR annotates each read with an information about the length of repeated regions. This information is used in subsequent analysis in [Dante](https://github.com/marcelTBI/dante-remaSTR). RemaSTR utilizes fast compiled language (Rust) and multithreading to calculate the annotations quickly, enabling Dante to scale to milions of reads and thousands of motifs.

## Installation
Assuming rust is [installed](https://www.rust-lang.org/tools/install), run:
```
cd remastr
cargo build --release
```
The resulting binary can be found in `target/release/remastr`.
The binary is statically linked, therefore it does not need any dynamic libraries and can be moved to any directory without dependency issues.

## Usage
To annotate reads:
```
remastr -f data/chromosomeX.fna -n data/HGVS.txt -b data/mini2.bam -o output.tsv
```

