import csv
import os
import sys

Input = sys.argv[1]
GenomeIndex = sys.argv[2]
Fastq = sys.argv[3]
Output = sys.argv[4]

Fastq1 = ","+Fastq+"/"

with open(Input,encoding="utf8") as single_end:
    tsv_readers = csv.DictReader(single_end, delimiter="\t")
    for single in tsv_readers:
        name = single["Name"]
        fastq = single["Fastq"]
        out1 = str(name).rstrip()
        out2 = out1+'Aligned.sortedByCoord.out.bam'

        if(',' in fastq):
            fq = fastq.split(",")
            cmd = "STAR --runThreadN 40 --genomeDir" + ' ' + GenomeIndex + ' '+"--readFilesIn" + ' ' + Fastq + "/"+ Fastq1.join(map(str,fq)) + ' ' +  "--readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000 --outBAMsortingThreadN 40 --outFileNamePrefix" + ' ' +Output+"/"+ out1+' '+"--genomeLoad NoSharedMemory --outSAMunmapped Within --outSAMstrandField intronMotif --outSJtype Standard"
            cmd1 = "samtools index"+ ' '+Output+"/"+out2
            os.system(cmd)
            os.system(cmd1)
        else:
            cmd2 = "STAR --runThreadN 40 --genomeDir" + ' ' + GenomeIndex + ' '+"--readFilesIn" + ' ' + Fastq + "/"+ fastq + ' ' +  "--readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000 --outBAMsortingThreadN 40 --outFileNamePrefix" + ' ' +Output+"/"+ out1+' '+"--genomeLoad NoSharedMemory --outSAMunmapped Within --outSAMstrandField intronMotif --outSJtype Standard"
            cmd21 = "samtools index"+ ' '+Output+"/"+out2
            os.system(cmd2)
            os.system(cmd21)
