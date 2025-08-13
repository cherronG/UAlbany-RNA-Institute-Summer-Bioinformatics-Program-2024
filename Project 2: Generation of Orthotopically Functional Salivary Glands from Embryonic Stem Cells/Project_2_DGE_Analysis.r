#Set working directory to my folder in bioinformatics lab

#Load DESeq2 package
library("DESeq2")

#DESeq2 Pipeline

#Count File Generation
file.list <- list.files( path = "./", pattern = "*ReadsPerGene.out.tab$")

counts.files <- lapply(file.list, read.table, skip = 4)

counts <- as.data.frame( sapply( counts.files, function(x) x[ ,2] ) )

colnames(counts) <- file.list

row.names(counts) <- counts.files[[1]]$V1

#Conditions and Metadata
condition <- c(rep("iSG",3), rep("iSG_mes",3))

sampleTable <- data.frame(sampleName = file.list, condition = condition)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~condition)

#Running and Writing Data
output <- DESeq(dds)

results_Control_DM <- results(output, contrast=c("condition","iSG_mes","iSG"))

results_Control_DM_PValue <- results_Control_DM[order(results_Control_DM$padj),]

head(results_Control_DM_PValue)

write.csv(as.data.frame(results_Control_DM_PValue), file="iSG_iSG_mes_DE.csv")

#Converting ENSEMBL IDs to Gene IDs

#load annotation package
library(AnnotationDbi)

#install ad load mouse gene ID reference package
BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")

iSG_iSG_mes <- read.table("iSG_iSG_mes_DE.csv", header = TRUE, sep =',')

IDs <- c(iSG_iSG_mes$X)

iSG_iSG_mes$Symbol <- mapIds(org.Mm.eg.db, IDs, 'SYMBOL', 'ENSEMBL')

rownames(iSG_iSG_mes) <- iSG_iSG_mes$X

write.csv(as.data.frame(iSG_iSG_mes), file="iSG_iSG_mes_DE_GeneIDs.csv", header = TRUE, sep =',')
write.csv(iSG_iSG_mes, file = "iSG_iSG_mes_DE_GeneIDs.csv")

# PCA Plot
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup=c("condition"))

#Filtered PCA Plot
library(magrittr) #load pipeline operator package
significant_events <- iSG_iSG_mes %>%
  filter(abs(iSG_iSG_mes$log2FoldChange)>1.5 & iSG_iSG_mes$pvalue < 0.05)

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd[significant_events],intgroup=c("condition"))

# MA Plot
plotMA(results_Control_DM, ylim=c(-10,10))

# Volcano Plot
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(results_Control_DM, lab = iSG_iSG_mes$Symbol, x = "log2FoldChange", y = "padj", pCutoff = 0.05, FCcutoff = 1.5, pointSize = 1.5, labSize = 3.0, title = "iSG+mes vs iSG", xlim = c(-10,10), ylim = c(0,3))

#Heatmap Plot
install.packages("pheatmap")
library("pheatmap")
rld <- rlog(dds, blind=FALSE)
topgenes <- head(rownames(results_Control_DM_PValue),10)
topgenes_symbols <- iSG_iSG_mes[topgenes,]
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)
rownames(mat) <- topgenes_symbols$Symbol
df <- as.data.frame(colData(dds)[,c("condition")])
colnames(df) <- "Condition"
row.names(df) <- sampleTable$sampleName
pheatmap(mat, annotation_col=df, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE)

#Normalized counts plot
plotCounts(dds, gene=which.min(results_Control_DM$padj), intgroup="condition") #Chst4
plotCounts(dds, gene="ENSMUSG00000060550", intgroup="condition") #H2-Q7
plotCounts(dds, gene="ENSMUSG00000027559", intgroup="condition") #Car3

