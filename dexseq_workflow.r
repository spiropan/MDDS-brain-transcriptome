# Reading in data 
inDir="/nfs/zorba/jjmdata/DESeq/"
countfiles=list.files(inDir,pattern="*dexseq_count.txt$",full.names=TRUE)
flattenedFile=list.files(inDir,pattern="genes.gff$",full.names=TRUE)
#sampleTable=read.csv('covs_gen_simple.csv')
sampleTable=read.csv('covs_age_gen_noSampleID.csv')
subset=1 # Set this to 1 if going to feed in a text file listing subset of genes

# create data analysis object
suppressPackageStartupMessages( library( "DEXSeq" ) )
print("Now creating data object")
dxd=DEXSeqDataSetFromHTSeq(
	countfiles,
	sampleData=sampleTable,
	design= ~ sample + exon + condition:exon,
	flattenedfile=flattenedFile)

sampleAnnotation( dxd )

print("Now estimating library sizes")
dxd=estimateSizeFactors(dxd)

# here to take a subset of Genes
if (subset > 0) {
print("Ok taking only a subset of genes")
genesForSubset = read.table(file.path(inDir,"Top_ID.txt"),stringsAsFactors=FALSE)[[1]]
#exonsForSubset = read.table(file.path(inDir,"TopExons_plot.txt"),stringsAsFactors=FALSE)[[1]]
exonsForSubset = read.table("TopExons2.txt",stringsAsFactors=FALSE)[[1]]
#genesForSubset = read.table(file.path(inDir,"SAT1_ID.txt"),stringsAsFactors=FALSE)[[1]]
#genesForSubset = read.table(file.path(inDir,"SAT1_ID.txt"),stringsAsFactors=FALSE)[[1]]
dxd=dxd[geneIDs( dxd ) %in% genesForSubset,]
dxd=dxd[exonIDs( dxd ) %in% exonsForSubset,]
}

print("Now estimating dispersions")
BPPARAM = MulticoreParam(workers=6)
formulaFullModel = ~ sample + exon + exon:condition + age:exon + gen:exon
formulaReducedModel = ~ sample + exon + age:exon + gen:exon

#dxd = estimateDispersions( dxd, formula = formulaFullModel, BPPARAM=BPPARAM)
dxd = estimateDispersions( dxd, formula = formulaFullModel, BPPARAM=BPPARAM)
#dxd = estimateDispersionsGeneEst( dxd )
#dxd = estimateDispersions( dxd, formula = formulaFullModel, fitType='mean' )

print("now testing for DEU")
# here include extra covs
dxd = testForDEU( dxd, reducedModel = formulaReducedModel, fullModel = formulaFullModel,BPPARAM=BPPARAM )
#dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)
dxr = DEXSeqResults( dxd )

DEXSeqHTML( dxr, FDR=0.1, BPPARAM=BPPARAM)
