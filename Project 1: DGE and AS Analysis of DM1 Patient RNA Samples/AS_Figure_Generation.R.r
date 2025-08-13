#Install and load maser package
BiocManager::install('maser')
library(maser)

#Download and load the rMATS data (using rMATS output folder)
DM1 <- maser("./", c("DM1", "Control"), ftype ="JCEC")

#Summary of the first 8 columns
head(summary(DM1, type = "SE")[, 1:8])

#Filters by read coverage using an average number of reads of 5
FilteredDM1 <- filterByCoverage(DM1, avg_reads = 5)

#Further filters using FDR and delta PSI
DM_top <- topEvents(FilteredDM1, fdr = 0.05, deltaPSI = 0.1)

#Unfiltered Events PCA Plot
pca(DM1, type = "SE")

#Filtered Events PCA Plot
pca(DM_top, type = "SE")

#Violin Plot
mbnl1_events <- geneEvents(DM1, geneS = "MBNL1", fdr = 0.05, deltaPSI = 0.1)
plotGenePSI(mbnl1_events, type = "SE", show_replicates = TRUE)

#Modified Sashimi Plot
gtf_path <- system.file("extdata", file.path("GTF","Ensembl85_examples.gtf.gz"), package = "maser")
ens_gtf <- rtracklayer::import.gff(gtf_path)
plotTranscripts(mbnl1_events, type = "SE", event_id = 59859, gtf = ens_gtf, zoom = FALSE, show_PSI = TRUE)
                