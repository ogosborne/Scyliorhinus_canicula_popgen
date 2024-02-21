library(pcadapt) # version: 4.3.5
library(snpStats) # version 1.44.0
# reformat data
path_to_file <- "results/var/filt.var.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
# choose optimal K
x <- pcadapt(input = filename, K = 10, LD.clumping = list(size = 500, thr = 0.1)) 
## diagnostic plots
# scree plot to choose optimal K
plot(x, option = "screeplot")
## optimal K = 2
K <- 2
# scores plot to see which PC differentiates populations
plot(x, option = "scores", pop = c(rep("A", 5), rep("M", 9)))
## PC1 only
PC = 1
# quantile-quantile plot to check presence of outliers 
plot(x, option = "qqplot", K = 1)
# plot loadings for all SNPs in order to check for clumping, if LD is responsible for outliers, we would expect them to be clumped. They aren't.
par(mfrow = c(1, 2))
for (i in 1:2)
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
# run with optimal K, for each PC separately
x <- pcadapt(filename, K = K, LD.clumping = list(size = 500, thr = 0.1), method = "componentwise")
# correct P-values, only take PC1 as this differentiates Atl vs Med
padj <- p.adjust(x$pvalues[,PC],method="bonferroni")
# load bed to get SNP positions
bed <- read.plink("results/var/filt.var.bed")
# load per site fst
fst <- read.table(file = "results/vcftools_fst/persite_fst.weir.fst", header = T)
# combine SNP positions with pcadapt
df <- data.frame(CHROM = bed$map$chromosome,
                 POS = bed$map$position,
                 pcadapt_chi2.stat = x$chi2.stat[,PC], 
                 pcadapt_p = x$pvalues[,PC], 
                 pcadapt_padj = padj)
# combine with per-site FST
df <- merge(df, fst, all.x = T, all.y = F) 
# set negative FST to 0
df[which(df$WEIR_AND_COCKERHAM_FST < 0), "WEIR_AND_COCKERHAM_FST"] <- 0
# save results
write.csv(df, file = "results/pcadapt/pcadapt.csv", row.names = F)
# make bed of significant SNPs to intersect with genes
dfsig <- df[which(df$pcadapt_padj < 0.05),]
write.table(data.frame(dfsig$CHROM, dfsig$POS, dfsig$POS), file = "results/pcadapt/sig_snps.bed", sep = "\t", row.names = F, col.names = F, quote = F)
