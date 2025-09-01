library("DESeq2")
library("RColorBrewer")
library("gplots")
library("vsn")
library("pheatmap")
library("ggplot2")
library("BiocParallel")

register(SnowParam(4))

#Download count tables and put in folder "AbATE_raw": https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA338797

#prepare the data
condition1 <- "Control_visit_0_6_12_24"
condition2 <- "Treated_responders_non_responders_visit_0_6_12_24"

print(paste("--R--DESeq2-- Condition ", condition1, " vs ", condition2))
print("--R--DESeq2-- Loading samples")
print(getwd())
all_samples = read.csv("samples.csv")
samples = subset(all_samples, condition %in% c(condition1, condition2))
allmisscols <- sapply(samples, function(x) all(is.na(x) | x == '' ))

if (allmisscols[4]==FALSE){
  dds = DESeqDataSetFromHTSeqCount(sampleTable=samples, directory = getwd(), design = ~ pairing + condition)} else{dds = DESeqDataSetFromHTSeqCount(sampleTable=samples, directory = paste0(getwd(),"/AbATE_raw"), design = ~ condition)}

#prefiltering
print("--R--DESeq2-- Prefiltering")
countsdds=counts(dds)
write.csv(as.data.frame(countsdds), file=paste("counts_results.csv", sep=""))
dds = dds[rowSums(counts(dds)) > 1,]
nsub=2

print("--R--DESeq2-- End Loading samples")

print("--R--DESeq-- Estimate Size Factors")
dds = estimateSizeFactors(dds)

write.csv(counts(dds, normalized=TRUE), paste("global", ".normalized_counts.csv", sep = ""))
