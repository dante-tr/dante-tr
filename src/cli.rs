use clap::Parser;

// Predict short tandem repeat annotation
#[derive(Parser)]
pub struct Args {
    /// Reference in FASTA format
    #[arg(short='f')]
    pub ref_file: String,
    
    /// Reads mapped to reference in BAM format, index (.bai) has to be present
    #[arg(short='b')]
    pub bam_file: String,

    /// Repeats in HGVS nomenclature, one per line
    #[arg(short='n')]
    pub hgvs_file: String,
    
    /// ID mapping from reference to bam file.
    #[arg(long="map1")]
    pub ref2bam: Option<String>,
    
    /// ID mapping from nomenclature to bam file.
    #[arg(long="map2")]
    pub hgvs2bam: Option<String>,
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn empty_args_prints_help() {
        let args = Args::try_parse_from(["remastr"].iter()).err().unwrap();
        println!("{}", args.to_string()); 
    }
    
    #[test]
    fn prints_help() {
        let args = Args::try_parse_from(["remastr", "-h"].iter()).err().unwrap();
        println!("{}", args.to_string());
    }

    #[test]
    fn cli_small_example() {
        let args = Args::try_parse_from([
            "remastr",
            "-f", "data/chromosomeX.fna",
            "-n", "data/mini_HGVS.txt",
            "-b", "data/mini2.bam"
        ].iter()).unwrap();
    }
}

