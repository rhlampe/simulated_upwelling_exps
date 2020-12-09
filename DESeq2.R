#
# DESeq2 boilerplate
# Manual: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#

library("DESeq2")

# Read in tab-delimited counts file
counts <- read.table('counts.tsv', row.names = 1, header = T)

# Read in experimental design
# Has columns for time points, iron level, and both combined
design <- read.csv('design.csv', row.names = 1)

# Create DDS object
# Create dds object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = design,
                              design = ~ combined)
dds <- DESeq(dds)

# Generate variance stablized counts
vsd <- vst(dds, blind=FALSE)
vsd_counts <- assay(vsd)
write.csv(vsd_counts, "vsd_counts.csv")


# Examine PCA plot
plotPCA(vsd, intgroup=c("timepoint", "fe_level"))

# Example results comparison for across time
res_high_2vs1 <- results(dds, contrast=c("combined", "tp2_high", "tp1_high"))

# Example results comparison for across Fe level
res_tp1_fe <- results(dds, contrast = c('combined', 'tp1_high', 'tp1_low'))
