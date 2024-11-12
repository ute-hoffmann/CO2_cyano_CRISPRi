# CO2_cyano_CRISPRi

This repository collects the analysis for the *Synechocystis* sp. PCC 6803 CRISPRi library (compare [Miao and Jahn et al., 2023](https://doi.org/10.1093/plcell/koad208)) when grown at 4% or 30% CO_2. For analysis, sequencing files were processed using the [nf-core-crispriscreen workflow](https://github.com/MPUSP/nf-core-crispriscreen).

Data is uploaded to [ShinyLib](https://m-jahn.shinyapps.io/ShinyLib/) as data set "CRISPRi_library_2024_CO2data" and can be explored interactively, 

## Documentation nf-core-crispriscreen

### Input files

Collect all relevant input files in [input](input).

* [Synechocystis_v2_trimmed.fasta](input/Synechocystis_v2_trimmed.fasta) is the library of used sgRNAs
* [samplesheet_CRISPRi_CO2_Elena.csv](input/samplesheet_CRISPRi_CO2_Elena.csv) gives which sample is which
* the respective .fastq.gz files are collected in the folder [input/fastq](input/fastq) and are the output of a sequencing run on a NextSeq 2000 system using NextSeq 1000/2000 P1 reagents (Illumina). The used workflow is "Illumina DRAGEN BCL Convert 3.10.4". Downloaded using the command "bs download run --id 273520278" (compare [Basespace CLI documentation](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview))

### nf-core-crispriscreen parameter decisions

Use control sgRNAs as point of comparision and updated verison of pipeline [e4aad5be10264d99e632761fc8fc56e68d6357c4](https://github.com/MPUSP/nf-core-crispriscreen/commit/e4aad5be10264d99e632761fc8fc56e68d6357c4):

```
conda activate env_nf
nextflow run ../../../pipelines/nf-core-crispriscreen/ -profile singularity --input "input/samplesheet_CRISPRi_CO2_Elena.csv" --fasta "input/Synechocystis_v2_trimmed.fasta" --outdir "results_controlsgRNAs" --three_prime_adapter ^CAGTGATAGAGATACTGGGAGCTA...GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC --filter_mapq=1 --max_cpus 5 --max_memory 12GB --run_mageck false --gene_controls "ctrl" --gene_fitness true
```

### Comments on results

Results of the nf-core-crispriscreen pipeline were saved in [results_controlssgRNAs](results_controlssgRNAs). The minimal amount of mapped reads was mapped for sample 4% CO2, generation 0, replicate 1. It was 1,967,391 reads for 21,470 barcodes. This gives a coverage of 91, which is a bit below the advised x100 coverage. 

### Exploratory data analysis / data analysis downstream of nf-core-crispriscreen

Further data processing was performed using R scripts provided in [R_code](R_code), results were saved in [R_results_controlssgRNAs](R_results_controlssgRNAs). Further R analysis included plotting and some diagnostics, annotation documented in [Visualization_EDA_with_adjp.pdf](R_code/Visualization_EDA_with_adjp.pdf).

File [Wilcoxon_compareConditions.Rmd](R_code/Wilcoxon_compareConditions.Rmd) was used for calculating adjusted p values comparing fitness values for sgRNAs targeting the same gene under different conditions. File [Create_RDataSet_ShinyLib.Rmd](R_code/Create_RDataSet_ShinyLib.Rmd) was used to format output file of nf-core-crispripipeline to be compatible with [https://m-jahn.shinyapps.io/ShinyLib/](https://m-jahn.shinyapps.io/ShinyLib/), compare [https://github.com/m-jahn/ShinyLib](https://github.com/m-jahn/ShinyLib). Finally, file [CompareToMiaoJahn.Rmd](R_code/CompareToMiaoJahn.Rmd) was used to compare data from this study to results obtained as part of [Miao and Jahn et al., 2023](https://doi.org/10.1093/plcell/koad208). 

## Brief description of remaining content of repository

Folder [content_AdriansRepository](content_AdriansRepository) contains files used by Adrian for analysis of the data set.
