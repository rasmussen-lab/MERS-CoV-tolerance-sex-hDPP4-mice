# MERS-CoV-tolerance-sex-hDPP4-mice
## Differential expression (DE) analysis of hDPP4 mice infected with a tolerant-inducing low dose or lethal high dose of MERS-CoV. DE analysis was performed both with and without dissagregation by sex, for days 1-5 of infection.

#### Step 1: Download the FASTQ files from NCBI BioProject PRJNA1276520
#### Step 2: Check the quality of the FASTQ files
##### FastQC version 0.12.1
```
fastqc -o RESULTS/QC -f fastq -t 20 FASTQ/*.gz

RESULTS/QC = Folder to save the output
FASTQ = Folder with FASTQ files in .gz format
```
#### Step 3: Download mouse genome and annotation
```
wget https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz
gunzip Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz
gunzip Mus_musculus.GRCm39.112.gtf.gz
```

#### Step 4: Build mouse genome index
##### STAR version 2.7.11b
```
mkdir GenomeIndex
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir GenomeIndex --genomeFastaFiles Mus_musculus.GRCm39.dna_sm.toplevel.fa --sjdbGTFfile Mus_musculus.GRCm39.112.gtf  --sjdbOverhang 99

GenomeIndex = Location where generated output files will be saved
```
#### Step 5: Aligning FASTQ reads to the mouse genome index
```
python Scripts/StarAlignS.py Data/Input.txt GenomeIndex Raw ALIGNMENTS

Input.txt = Tab delimited file listing paired-end FASTQ files
GenomeIndex =  Folder with mouse genome index
Raw = Folder with FASTQ files in .gz format
AlIGNMENTS = Folder to save the output files
```
#### Step 6: Generate read counts
##### R version 4.4.0 and Rsubread_2.18.0
```
Rscript Scripts/ReadCountSE.R
```
#### Step 7: Perform differential expression (DE) analysis for tolerance only
##### R version 4.4.0 and DESeq2_1.44.0 
```
Rscript Scripts/DiffExpr_Tolerance.R
```

#### Step 8: Perform DE analysis for tolerance and sex
##### R version 4.4.0 and DESeq2_1.44.0 
```
Rscript Scripts/DiffExpr_ToleranceSex.R
```

#### Step 9: Generate volcano plots of DE genes
```

```

#### Step 10: Perform multidimensional scaling (MDS) plot of DE genes
```

```
#### Step 11: Upload DE results to IPA

Uploaded DE lists generated in Step 8. to Ingenuity Pathway Analysis for functional enrichment analysis.

## Workflow Engineering & Implementation

### Core Systems Architect & Lead Developer
Dr. Reema Singh (res498@usask.ca) 

Responsible for the foundational technical architecture, logic-gated framework, and core analytical wrappers.

### Technical Infrastructure
This repository is built upon the standardized Workflow Engineering logic developed and validated by Dr. Reema Singh within the Rasmussen Lab. This architecture provides the modular framework and automated logic-gates required for high-throughput reproducibility of the MERS-CoV dataset.

### Study-Specific Implementation & Functional Mapping
#### Pipeline Execution: 
Marin Habbick (Responsible for the execution of the MERS-CoV sex-tolerance study pipeline, including parameter configuration and study-level data mapping.)
#### Functional Enrichment & IPA: 
Dr. Angela Rasmussen and Marin Habbick (Responsible for the downstream biological interpretation and Ingenuity Pathway Analysis mapping).

### Code Adaptations & Visualization
Implementation-level modifications: Marin Habbick

Note on AI-Assistance: In alignment with the manuscript's disclosure, specific study-level adaptations (Steps 9 & 10) and visualization components utilized AI-assisted code generation for implementation-specific mapping.

#### Validation & Reproducibility

Framework Stability: This pipeline has been validated across multiple datasets (Ebola, MERS-CoV) within the Rasmussen Lab to ensure cross-species logic consistency.

Logic-Gating: All filtering thresholds and batch-correction wrappers were benchmarked for high-fidelity transcriptomic recovery.
