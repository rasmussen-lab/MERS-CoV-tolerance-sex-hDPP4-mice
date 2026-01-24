# ==============================================================================
# DESIGN & IMPLEMENTATION: Dr. Reema Singh (res498@usask.ca)
# PROJECT: MERS-CoV DE Framework
# VERSION: 1.1.0 (Stable)
# ------------------------------------------------------------------------------
# DESCRIPTION: Multi-species transcriptomic pipeline. 
# Optimized for STAR alignment and DESeq2 statistical modeling.
# ==============================================================================

library(Rsubread)
setwd("ALIGNMENT")
file_list <- list.files(pattern="*.bam$")
Sample1 <- featureCounts(files=file_list,annot.ext="Mus_musculus.GRCm39.110.gtf",isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id",isPairedEnd=FALSE,nthreads=40)
write.table(Sample1$counts, file="readCountTDM.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
write.table(Sample1$stat, file="readCountStatTDM.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
