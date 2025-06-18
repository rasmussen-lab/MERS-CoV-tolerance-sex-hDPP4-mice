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
```
#### Step 5: Aligning FASTQ reads to the mouse genome index
```
python StarAlignS.py Input1.txt GenomeIndex Raw ALIGNMENTS

Input1.txt = Tab delimited file listing paired-end FASTQ files
GenomeIndex =  Folder with mouse genome index
Raw = Folder with FASTQ files in .gz format
AlIGNMENTS = Folder to save the output files
```
#### Step 6: Generate read counts
##### R version 4.4.0 and Rsubread_2.18.0
```

```
#### Step 7: Perform differential expression (DE) analysis for tolerance only
##### R version 4.4.0 and DESeq2_1.44.0 
```

```

#### Step 8: Perform DE analysis for tolerance and sex
##### R version 4.4.0 and DESeq2_1.44.0 
```

```
#### Step 8: Upload DE results to IPA
```

```
#### Step 9: Generate volcano plots of DE genes
```

```

#### Step 10: Perform multidimensional scaling (MDS) plot of DE genes
```

```
