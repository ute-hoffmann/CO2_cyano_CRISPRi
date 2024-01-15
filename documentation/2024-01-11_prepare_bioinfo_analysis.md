# Preparing input files

Collect all relevant input files in [input](../input).

* [Synechocystis_v2_trimmed.fasta](../input/Synechocystis_v2_trimmed.fasta) is the library of used sgRNAs
* [samplesheet_CRISPRi_CO2_Elena.csv](../input/samplesheet_CRISPRi_CO2_Elena.csv) gives which sample is which
* the respective .fastq.gz files are collected in the folder [../input/fastq](../input/fastq) and are the output of a sequencing run on a NextSeq 2000 system using NextSeq 1000/2000 P1 reagents (Illumina). The used workflow is "Illumina DRAGEN BCL Convert 3.10.4". Downloaded using the command "bs download run --id 273520278" (compare [Basespace CLI documentation](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview))

# nf-core-crispriscreen parameter decisions

Pipeline commit used for analysis: [c7c6563eb269bc3426e7481a8324788bc1006306](https://github.com/MPUSP/nf-core-crispriscreen/commit/c7c6563eb269bc3426e7481a8324788bc1006306) from 25th May 2023

```
conda activate env_nf
nextflow run ../../pipelines/nf-core-crispriscreen/ -profile singularity --input "input/samplesheet_CRISPRi_CO2_Elena.csv" --fasta "input/Synechocystis_v2_trimmed.fasta" --outdir "results" --three_prime_adapter ^CAGTGATAGAGATACTGGGAGCTA...GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC --filter_mapq=1 --max_cpus 5 --max_memory 12GB --run_mageck false --gene_fitness true
```
