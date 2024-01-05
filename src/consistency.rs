use std::collections::HashMap;
use crate::TandemRepeat;

pub fn ensure_consistency(
    bam_refs: HashMap<String, usize>,
    mut fna_refs: HashMap<String, Vec<u8>>,
    mut repeats: Vec<TandemRepeat>
) -> (HashMap<String, Vec<u8>>, Vec<TandemRepeat>) {

    let id = bam_refs.iter().next().expect("BAM should have at least one @SQ record.").0;
    let bam_variant = Variant::detect(id);
    let id = fna_refs.iter().next().expect("FASTA should have at least one record.").0;
    let fna_variant = Variant::detect(id);
    let id = &repeats.first().expect("HGVS should have at least one record.").reference;
    let hgvs_variant = Variant::detect(id);

    if (bam_variant.is_chr() || fna_variant.is_chr() || hgvs_variant.is_chr())
        && (bam_variant.is_nc() || fna_variant.is_nc() || hgvs_variant.is_nc())
    {
        eprintln!("WARNING: data is inconsistent. Results might be incorrect.");
        eprintln!("Attempting ID remapping.");
        if bam_variant.is_nc() { panic!("Cannot remap IDs in BAM."); }
        if fna_variant.is_nc() {
            eprintln!("Remapping fasta records.");
            fna_refs = remap_fna(fna_refs);
        }
        if hgvs_variant.is_nc() {
            eprintln!("Remapping nomenclatures.");
            repeats = remap_hgvs(repeats);
        }
    }
    return (fna_refs, repeats);
}

#[derive(PartialEq, Eq)]
enum Variant { Chr, NC, Other }

impl Variant {
    fn is_chr(&self) -> bool { *self == Variant::Chr }
    fn is_nc(&self) -> bool { *self == Variant::NC }
    fn detect(id: &String) -> Variant {
        if id.starts_with("chr") { return Variant::Chr; }
        if id.starts_with("NC") { return Variant::NC; }
        return Variant::Other;
    }
}

const MAP: [(&str, &str); 24] = [
    ("NC_000001", "chr1"),  ("NC_000002", "chr2"),  ("NC_000003", "chr3"),  ("NC_000004", "chr4"),
    ("NC_000005", "chr5"),  ("NC_000006", "chr6"),  ("NC_000007", "chr7"),  ("NC_000008", "chr8"),
    ("NC_000009", "chr9"),  ("NC_000010", "chr10"), ("NC_000011", "chr11"), ("NC_000012", "chr12"),
    ("NC_000013", "chr13"), ("NC_000014", "chr14"), ("NC_000015", "chr15"), ("NC_000016", "chr16"),
    ("NC_000017", "chr17"), ("NC_000018", "chr18"), ("NC_000019", "chr19"), ("NC_000020", "chr20"),
    ("NC_000021", "chr21"), ("NC_000022", "chr22"), ("NC_000023", "chrX"),  ("NC_000024", "chrY"),
];

fn remap_fna(fna_refs: HashMap<String, Vec<u8>>) -> HashMap<String, Vec<u8>> {
    let mut result = HashMap::new();
    for (id, seq) in fna_refs.iter() {
        let id = id.split('.').next()
            .expect("Fasta ID is not in form NC_000000.0");
        for i in 0..MAP.len() {
            if id == MAP[i].0 {
                result.insert(MAP[i].1.to_string(), seq.clone());
            }
        }
    }
    return result;
}

fn remap_hgvs(repeats: Vec<TandemRepeat>) -> Vec<TandemRepeat> {
    let mut result = Vec::new();
    for repeat in repeats {
        let name = repeat.reference;
        let id = name.split('.').next()
            .expect(&format!("{name} is not in form NC_000000.0"));
        for i in 0..MAP.len() {
            if id == MAP[i].0 {
                result.push(TandemRepeat{
                    reference: MAP[i].1.to_string(),
                    start: repeat.start,
                    end: repeat.end,
                    copy_unit: repeat.copy_unit.clone(),
                    copy_number: repeat.copy_number.clone()
                })
            }
        }
    }
    return result;
}

