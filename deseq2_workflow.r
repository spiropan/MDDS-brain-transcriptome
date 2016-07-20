# Reading in data 
suppressPackageStartupMessages( library( "DESeq2" ) )
gtfDir="/nfs/zorba/jjmdata/spiro/"
inDir="/nfs/zorba/jjmdata/DESeq/"
create_se = 0 #if 1 will create the counts matrix for all samples if not already done

if (create_se == 1) { 
suppressPackageStartupMessages( library( "Rsamtools" ) )
# Here construct the TranscriptDB
suppressPackageStartupMessages( library("GenomicFeatures"))
gtffile <- list.files(gtfDir,pattern="genes.gtf$",full.names=TRUE)
print("Now creating transcript DB from GTF file")
(txtdb <- makeTranscriptDbFromGFF(gtffile, format="gtf"))
print("Now group by gene")
(genes <- exonsBy(txtdb, by="gene"))

print("Now performing the counting")
library("GenomicAlignments")
library('BiocParallel')
register(MulticoreParam(workers=8))
filenames=list.files(inDir,pattern="*.bam$",full.names=TRUE)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
se <- summarizeOverlaps(genes, bamfiles, mode="Union", singleEnd=FALSE,ignore.strand=TRUE)
save.image(file = "DeSeq2.Rdata" ) 
}

load(file="DeSeq2.Rdata")
sampleTable=read.csv('Covs.csv') # This to compute LRT testing for overall differences between 3 groups
sampleTable1=read.csv('Covs_NonSuivsSui.csv') # Thes to compute depression and suicide contrasts, 2 groups
sampleTable2=read.csv('Covs_MDDvsNonMDD.csv')
se1=se
se2=se

(colData(se) <- DataFrame(sampleTable))
(colData(se1) <- DataFrame(sampleTable1))
(colData(se2) <- DataFrame(sampleTable2))

dds <- DESeqDataSet(se = se, design = ~ age + sex + condition)
dds1 <- DESeqDataSet(se = se1, design = ~ age + sex + condition)
dds2 <- DESeqDataSet(se = se2, design = ~ age + sex + condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)

# Here compute results while including RIN as covariate
seRIN <- se[, which(colData(se)$RIN > 0)]
ddsRIN <- DESeqDataSet(se = seRIN, design = ~ age + sex + RIN + condition)
ddsRIN <- DESeq(ddsRIN, reduced = ~ age + sex + RIN, test = "LRT")
resRINOrdered <- resRIN[order(resRIN$padj),]

# Here compute the depression and suicide contrasts
tmp <- results(dds1, contrast=c("condition","Sui","NonSui"), altHypothesis="greater")
SuivsNonSui <- tmp[order(tmp$padj),]
tmp <- results(dds1, contrast=c("condition","Sui","NonSui"), altHypothesis="less")
NonSuivsSui <- tmp[order(tmp$padj),]

tmp <- results(dds2, contrast=c("condition","MDD","NonMDD"), altHypothesis="greater")
MDDvsNonMDD <- tmp[order(tmp$padj),]
tmp <- results(dds2, contrast=c("condition","MDD","NonMDD"), altHypothesis="less")
NonMDDvsMDD <- tmp[order(tmp$padj),]

# Here compute the pairwise contrasts with 10% min fold change threshld
dds <- DESeq(dds)
tmp <- results(dds, contrast=c("condition","MDD-S","CON"), lfcThreshold=0.138, altHypothesis="greater")
MDD_SvsCON <- tmp[order(tmp$padj),]
tmp <- results(dds, contrast=c("condition","MDD-S","CON"), lfcThreshold=0.138, altHypothesis="less")
CONvsMDD_S <- tmp[order(tmp$padj),]

tmp <- results(dds, contrast=c("condition","MDD-S","MDD"), lfcThreshold=0.138, altHypothesis="greater")
MDD_SvsMDD <- tmp[order(tmp$padj),]
tmp <- results(dds, contrast=c("condition","MDD-S","MDD"), lfcThreshold=0.138, altHypothesis="less")
MDDvsMDD_S <- tmp[order(tmp$padj),]

tmp <- results(dds, contrast=c("condition","CON","MDD"), lfcThreshold=0.138, altHypothesis="greater")
CONvsMDD <- tmp[order(tmp$padj),]
tmp <- results(dds, contrast=c("condition","CON","MDD"), lfcThreshold=0.138, altHypothesis="less")
MDDvsCON <- tmp[order(tmp$padj),]
 

# Here write out .csv files
write.csv(as.data.frame(resRINOrdered),file="DESeq2_LRT_RIN_results.csv")
write.csv(as.data.frame(SuivsNonSui),file="DESeq2_SuivsNonSui_results.csv")
write.csv(as.data.frame(NonSuivsSui),file="DESeq2_NonSuivsSui_results.csv")
write.csv(as.data.frame(MDDvsNonMDD),file="DESeq2_MDDvsNonMDD_results.csv")
write.csv(as.data.frame(NonMDDvsMDD),file="DESeq2_NonMDDvsMDD_results.csv")

write.csv(as.data.frame(MDD_SvsCON),file="DESeq2_MDD_SgtCON_results.csv")
# write.csv(as.data.frame(MDD_SvsCON_RIN),file="DESeq2_MDD_SgtCON_RIN_results.csv")
write.csv(as.data.frame(CONvsMDD_S),file="DESeq2_CONgtMDD_S_results.csv")
# write.csv(as.data.frame(CONvsMDD_S_RIN),file="DESeq2_CONgtMDD_S_RIN_results.csv")

write.csv(as.data.frame(MDD_SvsMDD),file="DESeq2_MDD_SgtMDD_results.csv")
# write.csv(as.data.frame(MDD_SvsMDD_RIN),file="DESeq2_MDD_SgtMDD_RIN_results.csv")
write.csv(as.data.frame(MDDvsMDD_S),file="DESeq2_MDDgtMDD_S_results.csv")
# write.csv(as.data.frame(MDDvsMDD_S_RIN),file="DESeq2_MDDgtMDD_S_RIN_results.csv")

write.csv(as.data.frame(CONvsMDD),file="DESeq2_CONgtMDD_results.csv")
# write.csv(as.data.frame(CONvsMDD_RIN),file="DESeq2_CONgtMDD_RIN_results.csv")
write.csv(as.data.frame(MDDvsCON),file="DESeq2_MDDgtCON_results.csv")
# write.csv(as.data.frame(MDDvsCON_RIN),file="DESeq2_MDDgtCON_RIN_results.csv")


