#module load r/4.4.0

library(DESeq2)
library(vsn)
library(ggplot2)

Counts <- read.table(file="readCountTDM_MERS.txt",sep="\t", header = TRUE, check.names = FALSE)
names(Counts) <- gsub("Aligned.sortedByCoord.out.bam", "", names(Counts))
StudyDes <- read.table(file="StudyDesignComp2.txt",sep="\t", header = TRUE)

ddsMat <- DESeqDataSetFromMatrix(countData = Counts, colData = StudyDes, design = ~Collection_Date)
nrow(ddsMat)
keep <- rowSums(counts(ddsMat) >= 15) >= 6
ddsMat <- ddsMat[keep,]
nrow(ddsMat)

vsd <- vst(ddsMat, blind = FALSE)
head(assay(vsd))
sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(StudyDes$Dose_Group, StudyDes$Collection_Date, sep="-")


mds <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mds, colData(vsd))
pdf("Figure1_MDSPlot_DosegroupVsDate.pdf")
p <- qplot(X1,X2,color=Dose_Group,shape=Sex,data=as.data.frame(mds),xlab = "Dimension1", ylab = "Dimension2")+geom_point(size =3)
p + theme_classic()+theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=9), axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=13))
dev.off()

ddsMat <- DESeq(ddsMat)
write.csv(as.data.frame(counts(ddsMat)),file="MERSRawCount.csv")

##### Differential expression comparisons

res1L_F <- results(ddsMat, alpha = 0.05,contrast=c("Collection_Date", "Day1_LOW_F","Day1_MOCK_F"))
mcols(res1L_F, use.names = TRUE)
summary(res1L_F)
res1L_FOrdered <- res1L_F[order(res1L_F$pvalue),]
write.csv(as.data.frame(res1L_FOrdered), file="Day1_LOWvsMOCK_F.csv")

res1H_F <- results(ddsMat, alpha = 0.05, contrast=c("Collection_Date", "Day1_HIGH_F", "Day1_MOCK_F"))
mcols(res1H_F, use.names = TRUE)
summary(res1H_F)
res1H_FOrdered <- res1H_F[order(res1H_F$pvalue),]
write.csv(as.data.frame(res1H_FOrdered), file="Day1_HIGHvsMOCK_F.csv")

res1L_M <- results(ddsMat, alpha = 0.05,contrast=c("Collection_Date", "Day1_LOW_M","Day1_MOCK_M"))
mcols(res1L_M, use.names = TRUE)
summary(res1L_M)
res1L_MOrdered <- res1L_M[order(res1L_M$pvalue),]
write.csv(as.data.frame(res1L_MOrdered), file="Day1_LOWvsMOCK_M.csv")

res1H_M <- results(ddsMat, alpha = 0.05, contrast=c("Collection_Date", "Day1_HIGH_M", "Day1_MOCK_M"))
mcols(res1H_M, use.names = TRUE)
summary(res1H_M)
res1H_MOrdered <- res1H_M[order(res1H_M$pvalue),]
write.csv(as.data.frame(res1H_MOrdered), file="Day1_HIGHvsMOCK_M.csv")

### Day 2
res2L_F <- results(ddsMat, alpha = 0.05,contrast=c("Collection_Date", "Day2_LOW_F","Day2_MOCK_F"))
mcols(res2L_F, use.names = TRUE)
summary(res2L_F)
res2L_FOrdered <- res2L_F[order(res2L_F$pvalue),]
write.csv(as.data.frame(res2L_FOrdered), file="Day2_LOWvsMOCK_F.csv")

res2H_F <- results(ddsMat, alpha = 0.05, contrast=c("Collection_Date", "Day2_HIGH_F", "Day2_MOCK_F"))
mcols(res2H_F, use.names = TRUE)
summary(res2H_F)
res2H_FOrdered <- res2H_F[order(res2H_F$pvalue),]
write.csv(as.data.frame(res2H_FOrdered), file="Day2_HIGHvsMOCK_F.csv")

res2L_M <- results(ddsMat, alpha = 0.05,contrast=c("Collection_Date", "Day2_LOW_M","Day2_MOCK_M"))
mcols(res2L_M, use.names = TRUE)
summary(res2L_M)
res2L_MOrdered <- res2L_M[order(res2L_M$pvalue),]
write.csv(as.data.frame(res2L_MOrdered), file="Day2_LOWvsMOCK_M.csv")

res2H_M <- results(ddsMat, alpha = 0.05, contrast=c("Collection_Date", "Day2_HIGH_M", "Day2_MOCK_M"))
mcols(res2H_M, use.names = TRUE)
summary(res2H_M)
res2H_MOrdered <- res2H_M[order(res2H_M$pvalue),]
write.csv(as.data.frame(res2H_MOrdered), file="Day2_HIGHvsMOCK_M.csv")

# Day 3

res3L_F <- results(ddsMat, alpha = 0.05,contrast=c("Collection_Date", "Day3_LOW_F","Day3_MOCK_F"))
mcols(res3L_F, use.names = TRUE)
summary(res3L_F)
res3L_FOrdered <- res3L_F[order(res3L_F$pvalue),]
write.csv(as.data.frame(res3L_FOrdered), file="Day3_LOWvsMOCK_F.csv")

res3H_F <- results(ddsMat, alpha = 0.05, contrast=c("Collection_Date", "Day3_HIGH_F", "Day3_MOCK_F"))
mcols(res3H_F, use.names = TRUE)
summary(res3H_F)
res3H_FOrdered <- res3H_F[order(res3H_F$pvalue),]
write.csv(as.data.frame(res3H_FOrdered), file="Day3_HIGHvsMOCK_F.csv")

res3L_M <- results(ddsMat, alpha = 0.05,contrast=c("Collection_Date", "Day3_LOW_M","Day3_MOCK_M"))
mcols(res3L_M, use.names = TRUE)
summary(res3L_M)
res3L_MOrdered <- res3L_M[order(res3L_M$pvalue),]
write.csv(as.data.frame(res3L_MOrdered), file="Day3_LOWvsMOCK_M.csv")

res3H_M <- results(ddsMat, alpha = 0.05, contrast=c("Collection_Date", "Day3_HIGH_M", "Day3_MOCK_M"))
mcols(res3H_M, use.names = TRUE)
summary(res3H_M)
res3H_MOrdered <- res3H_M[order(res3H_M$pvalue),]
write.csv(as.data.frame(res3H_MOrdered), file="Day3_HIGHvsMOCK_M.csv")

# Day 4

res4L_F <- results(ddsMat, alpha = 0.05,contrast=c("Collection_Date", "Day4_LOW_F","Day4_MOCK_F"))
mcols(res4L_F, use.names = TRUE)
summary(res4L_F)
res4L_FOrdered <- res4L_F[order(res4L_F$pvalue),]
write.csv(as.data.frame(res4L_FOrdered), file="Day4_LOWvsMOCK_F.csv")

res4H_F <- results(ddsMat, alpha = 0.05, contrast=c("Collection_Date", "Day4_HIGH_F", "Day4_MOCK_F"))
mcols(res4H_F, use.names = TRUE)
summary(res4H_F)
res4H_FOrdered <- res4H_F[order(res4H_F$pvalue),]
write.csv(as.data.frame(res4H_FOrdered), file="Day4_HIGHvsMOCK_F.csv")

res4L_M <- results(ddsMat, alpha = 0.05,contrast=c("Collection_Date", "Day4_LOW_M","Day4_MOCK_M"))
mcols(res4L_M, use.names = TRUE)
summary(res4L_M)
res4L_MOrdered <- res4L_M[order(res4L_M$pvalue),]
write.csv(as.data.frame(res4L_MOrdered), file="Day4_LOWvsMOCK_M.csv")

res4H_M <- results(ddsMat, alpha = 0.05, contrast=c("Collection_Date", "Day4_HIGH_M", "Day4_MOCK_M"))
mcols(res4H_M, use.names = TRUE)
summary(res4H_M)
res4H_MOrdered <- res4H_M[order(res4H_M$pvalue),]
write.csv(as.data.frame(res4H_MOrdered), file="Day4_HIGHvsMOCK_M.csv")

# Day 5

res5L_F <- results(ddsMat, alpha = 0.05,contrast=c("Collection_Date", "Day5_LOW_F","Day5_MOCK_F"))
mcols(res5L_F, use.names = TRUE)
summary(res5L_F)
res5L_FOrdered <- res5L_F[order(res5L_F$pvalue),]
write.csv(as.data.frame(res5L_FOrdered), file="Day5_LOWvsMOCK_F.csv")

res5H_F <- results(ddsMat, alpha = 0.05, contrast=c("Collection_Date", "Day5_HIGH_F", "Day5_MOCK_F"))
mcols(res5H_F, use.names = TRUE)
summary(res5H_F)
res5H_FOrdered <- res5H_F[order(res5H_F$pvalue),]
write.csv(as.data.frame(res5H_FOrdered), file="Day5_HIGHvsMOCK_F.csv")

res5L_M <- results(ddsMat, alpha = 0.05,contrast=c("Collection_Date", "Day5_LOW_M","Day5_MOCK_M"))
mcols(res5L_M, use.names = TRUE)
summary(res5L_M)
res5L_MOrdered <- res5L_M[order(res5L_M$pvalue),]
write.csv(as.data.frame(res5L_MOrdered), file="Day5_LOWvsMOCK_M.csv")

res5H_M <- results(ddsMat, alpha = 0.05, contrast=c("Collection_Date", "Day5_HIGH_M", "Day5_MOCK_M"))
mcols(res5H_M, use.names = TRUE)
summary(res5H_M)
res5H_MOrdered <- res5H_M[order(res5H_M$pvalue),]
write.csv(as.data.frame(res5H_MOrdered), file="Day5_HIGHvsMOCK_M.csv")

