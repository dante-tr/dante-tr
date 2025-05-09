# reads
- info about read sets: https://github.com/genome-in-a-bottle/giab_data_indexes
- raw reads have 560GB, 300x cov, and were used to create VCFs from giabtr_comparison
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/Leicester_HipSTR_GangSTR_05182022/HipSTR/
- downloaded from:
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/
    HG002.GRCh38.300x.bam
- filtered by 
`samtools view -b -h -L diseases_clean.bed HG002.GRCh38.300x.bam > HG002.GRCh38.selected.bam`

