Instructions based mainly on instructions by Michael Jahn

# Submission of sequencing data to ENA - European Nucleotide Archive

This file is a mark down document and can be rendered to HTML, PDF and other formats using for example Rstudio. The entries were used for a file submission of TnSeq files, but it can also serve as template for future submissions.

## Submission How-To

1. upload individual samples using FTP, see details below
2. login to ENA web client at https://www.ebi.ac.uk/ena/submit/sra/#home
3. create new study or select existing one
4. enter study details 
5. enter sample details (alternatively upload a premade table)
6. enter run details  (alternatively upload a premade table)
7. finalize submission


----------

### FTP file submission

Upload via FTP is the best option. Data must be uploaded **before** entering any study/sample details, so that the system can map it. First start the local FTP program in linux terminal:

```
ftp webin2.ebi.ac.uk
```

Login data for ENA is saved in keepass

Enter user name and password. Then enter and execute the following commands (in this order).
**Important note: Switch to passive mode on client side** if the following error message appears:
`200 PORT command successful. Consider using PASV. 425 Failed to establish connection.` after trying to upload a file. If `passive` switches the client to active, just repeat the command to switch back to passive. The command `mput <file-pattern>` (= multiple file upload) should match the file names that are supposed to be uploaded.

```
ftp> bin
ftp> prompt
ftp> passive
ftp> mput *.fastq.gz
```

Comment 2024-09-02: Somehow worked without all the commands above, press "a" for "all" when asked [anpqy?]. "prompt" would suppress re-iterating questions.

### Create new study (optional)

**Release date**

some future day

**Short name of study**

Growth of CRISPRi library of Synechocystis at 4% and 30% CO2 gas feed

**Full name of study**

Stress responses in Synechocystis sp. PCC 6803 to very high CO2 concentrations revealed using a CRISPR interference mutant library

**Keywords**

Synechocystis, CRISPRi, CO2, gas feed, growth conditions

**Abstract**

Optimal cultivation of photosynthetic microbes is often performed under high CO2 conditions to promote growth, but excessively high CO2 levels inhibit growth through mechanisms that are not fully understood. This study investigates the physiological responses of the cyanobacterium Synechocystis sp. PCC 6803 to very high CO2 concentrations under conditions where pH is maintained at around 7.5. The growth rate of the wild type (WT) at a gas phase content of 30% CO2 was 2.6-fold lower compared to 4% CO2. Using a CRISPR interference knockdown mutant library, we identified mutations that either enhanced or impaired growth under these conditions. The results showed that the knockdown of phycocyanin biosynthesis genes (cpc) promoted growth in 30% CO2 suggesting that the cells were under high light stress under very high CO2. Conversely, knockdown of key regulators of CO2 fixation (cp12 and rbcR) and carbon utilization (pmgA) resulted in increased growth inhibition at 30% CO2, suggesting that adaptation to very high CO2 conditions requires alterations in carbon fixation and utilization pathways. Further experiments confirmed that WT cells under 30% CO2 experienced light stress at lower light intensities than WT cells under 4% CO2. In addition, experiments with knockout mutants lacking cpcG1 showed better growth than WT at 30% CO2, while pmgA and cp12 mutants performed worse. These findings suggest that enhanced fitness under very high CO2 involves modifications in light harvest, CO2 fixation, and carbon management pathways. The genetic profiling in this study provides potential targets for engineering cyanobacteria with improved photosynthetic efficiency and stress resilience, which could be valuable for biotechnological applications.

### Sample annotation

A tab-separated table in `*.tsv` format has to be prepared for sample annotation. The file for the 1st cultivation large library is ENA_defaultSampleTable_CRISPRi_CO2.tsv

### Run annotation

A tab-separated table in `*.tsv` format has to be prepared with file (run) names and descriptions. For the first cultivation large library, the file is called fastq1_template_1725275624938_CRISPRi_CO2.tsv

### Finalize submission

Push the submit button. For citing data in pulbications, use project accession number PRJEB79744
