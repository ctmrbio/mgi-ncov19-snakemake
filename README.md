# mgi-ncov19-snakemake

This snakemake pipeline can conduct cov-19 virus classification, de novo assembly, coverage assessment and variant calling.

The pipeline is built according to https://github.com/BGI-IORI/nCoV_Meta (preprint: https://doi.org/10.1101/2020.03.16.993584)

Differences between mgi-ncov19-snakemake and nCoV_Meta: 
  1) low complexity reads removal were implemented with fastp (bgi: prinSEQ)
  2) kraken2 was employed in this snakemake pipeline (bgi: kraken1)  
  3) SOAPnuke v2 (bgi: SOAPnuke v1) (better to change to SOAPnuke v1)
  4) not yet finished with the alignment and variant calling steps.
  
  updated: 2020-04-14
  
## Usage:
### 0. Install Conda and Snakemake 


### 1. Clone workflow
    git clone git@github.com:huyue87/mgi-ncov19-snakemake

### 2. Execute workflow
    # 2.1 load input files (paired-end raw reads)
    cd mgi-ncov19-snakemake/input
    ln -s Sample_{1,2}.fq.gz . 
    
    # 2.2 run de novo assembly and generating sam files
    cd mgi-ncov19-snakemake
    snakemake --use-conda -n 
    snakemake --use-conda 
    
    
