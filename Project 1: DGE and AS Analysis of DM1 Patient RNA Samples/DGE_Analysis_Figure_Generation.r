#Set working directory: select "Set Working Directory" in "Sessions" tab then "Choose Directory"

#Download DESeq2 package (only done once)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

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
condition <- c(rep("Control",5), rep("DM1",5))

sampleTable <- data.frame(sampleName = file.list, condition = condition)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~condition)

#Running and Writing Data
output <- DESeq(dds)

results_Control_DM <- results(output, contrast=c("condition","DM1","Control"))

results_Control_DM_PValue <- results_Control_DM[order(results_Control_DM$padj),]

head(results_Control_DM_PValue)

write.csv(as.data.frame(results_Control_DM_PValue), file="Control_DM_DE.csv")

#Converting ENSEMBL IDs to Gene IDs
BiocManager::install("AnnotationDbi") 

BiocManager::install("org.Hs.eg.db")

library(AnnotationDbi)

library(org.Hs.eg.db)

DM1_Control <- read.table("Control_DM_DE.csv", header = TRUE, sep =',')

IDs <- c(DM1_Control$X)

DM1_Control$Symbol <- mapIds(org.Hs.eg.db, IDs, 'SYMBOL', 'ENSEMBL')

rownames(DM1_Control) <- DM1_Control$X

write.csv(as.data.frame(DM1_Control), file="Control_DM_DE_GeneIDs.csv", header = TRUE, sep =',')
write.csv(DM1_Control, file = "Control_DM_DE_GeneIDs.csv") 

# PCA Plot
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup=c("condition"))

#Filtered PCA Plot
significant_events <- Control_DM_DE %>%
  filter(abs(Control_DM_DE$log2FoldChange) > 1.5 & Control_DM_DE$pvalue < 0.05)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd[significant_events$NameOfEnsembleIDColumn,],intgroup=c("condition"))

# MA Plot
plotMA(results_Control_DM, ylim=c(-10,10))

# Volcano Plot
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(results_Control_DM, lab = DM1_Control$Symbol, x = "log2FoldChange", y = "padj", pCutoff = 0.05, FCcutoff = 1.5, pointSize = 1.5, labSize = 3.0, title = "DM1 vs Control", xlim = c(-10,10), ylim = c(0,10))

#Heatmap Plot
install.packages("pheatmap")
library("pheatmap")
rld <- rlog(dds, blind=FALSE)
topgenes <- head(rownames(results_Control_DM_PValue),20)
topgenes_symbols <- DM1_Control[topgenes,]
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)
rownames(mat) <- topgenes_symbols$Symbol
df <- as.data.frame(colData(dds)[,c("condition")])
colnames(df) <- "Condition"
row.names(df) <- sampleTable$sampleName
pheatmap(mat, annotation_col=df, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE)

#Normalized counts plot
plotCounts(dds, gene=which.min(results_Control_DM$padj), intgroup="condition")

plotCounts(dds, gene=which.min(results_Control_DM$padj), intgroup="condition")
