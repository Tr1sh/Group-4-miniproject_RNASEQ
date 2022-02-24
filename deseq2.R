# load libraries ----------------------------------------------------------
library(DESeq2)
library(tximport)
library(rhdf5)
library(ggplot2)
library(apeglm)
library(org.Hs.eg.db)


# import counts data -------------------------------------------
countdata <- read.table("./hisat2_counts.txt", header = TRUE, skip = 1, row.names = 1)


# Remove length/char columns
countdata <- countdata[ ,c(-1:-5)]

# rename columns of countdata with sample names
colnames(countdata) <- paste0("sample", 37:42)

# load and process metadata -----------------------------------------------
metadata <- read.delim("practice.dataset.metadata.tsv", row.names = 1)

# deseq2 ------------------------------------------------------------------
# compare colnames of count data to rownames of metadata
colnames(countdata) == rownames(metadata)

# create deseq2 data object
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = metadata,
                                 design = ~Condition)


# Find differential expressed genes
ddsMat <- DESeq(ddsMat)

# obtain statistically results 
results <- results(ddsMat)

# head of results
head(results)

# Generate summary of testing. 
summary(results)

# store results in text file
# store all results
write.table(x = as.data.frame(results), 
            file = "results.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# store only statistically significant results(padj < 0.05)
results_sig <- subset(results, padj < 0.05)

write.table(x = as.data.frame(results_sig), 
            file = "results_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)


# plots of gene expression data -------------------------------------------------------------------
# Convert all samples to rlog(regularized logarithm) for visualization
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# head
head(assay(ddsMat_rlog))

# pca
plotPCA(ddsMat_rlog, intgroup = "Condition") +
  ggtitle(label = "Principal Component Analysis (PCA)") +
  scale_y_continuous(limits = c(-40, 40)) +
  scale_x_continuous(limits = c(-70, 40)) +
  theme_light()

# ma plot
# remove noise using apeglm
resultslfc <- lfcShrink(ddsMat, coef="Condition_normal_vs_disease", type="apeglm")

#  plot ma 
plotMA(resultslfc, ylim=c(-5, 5))


# import kallisto data ---------------------------------------------------------
# create path to the abundance files
# sample names
samples <- paste0("sample", 37:42)

# file path
files <- file.path(".", samples, "abundance.h5")

# file names
names(files) <- paste0("sample", 37:42)

# import abundance data
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE, ignoreAfterBar = TRUE)

# view first few lines of count data
head(txi.kallisto$counts)

# create deseq2 data object
kallisto_dds <- DESeqDataSetFromTximport(txi.kallisto,
                                         colData = metadata,
                                         design = ~ Condition)

# Find differential expressed genes
kallisto_diff <- DESeq(kallisto_dds)

# obtain results
kallisto_results <- results(kallisto_diff)

# results
head(kallisto_results)

# summary
summary(kallisto_results)
