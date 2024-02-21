library(stringr) # version 1.5.1
library(qqman) # version 0.1.9
library(calibrate) # version 1.7.7
# modified version of qqman::manhattan
source("code/manhattan2.R")
# load data
# pcadapt
pcadapt <- read.csv("results/pcadapt/pcadapt.csv")

# pcadapt gene matches
pcadapt_genes <- read.table("results/pcadapt/sig_snp_genes.tsv", sep = "\t")
pcadapt_genes <- pcadapt_genes[, c(1,2,5,6,9,13)]
colnames(pcadapt_genes) <- c("chrom", "snp_pos", "gene_start", "gene_end", "strand", "info")
pcadapt_genes$gene_name <- str_match(pcadapt_genes$info, ";Name=\\s*(.*?)\\s*;")[,2]
# remove genes without proper names
pcadapt_genes <- pcadapt_genes[grep(pattern = "LOC", pcadapt_genes$gene_name, invert = T),]
pcadapt_genes <- pcadapt_genes[grep(pattern = ":", pcadapt_genes$gene_name, invert = T),]
# add gene names to pcadapt table
pcadapt$gene_name <- pcadapt_genes[match(paste0(pcadapt$CHROM, "_", pcadapt$POS), paste0(pcadapt_genes$chrom, "_", pcadapt_genes$snp_pos)),"gene_name"]
# keep only unique genes
pcadapt_genes <- pcadapt_genes[!duplicated(pcadapt_genes$gene_name),]
# add gene names as snp names to the first snp from each gene for labelling
pcadapt$snp_name <- pcadapt_genes[match(paste0(pcadapt$CHROM, "_", pcadapt$POS), paste0(pcadapt_genes$chrom, "_", pcadapt_genes$snp_pos)),"gene_name"]
# make snp names for genes without one
pcadapt[which(is.na(pcadapt$snp_name)), "snp_name"] <- paste("SNP",pcadapt[which(is.na(pcadapt$snp_name)),"CHROM"], pcadapt[which(is.na(pcadapt$snp_name)),"POS"], sep = "_")
# save pcadapt results
write.csv(pcadapt, file = "results/pcadapt/pcadapt_res.csv", row.names = F, quote = F)
# keep only sig snps
pcadapt <- pcadapt[which(pcadapt$pcadapt_padj < 0.05),]
# remake df to bind with main df
pcadapt <- data.frame(chromosome = pcadapt$CHROM, window_pos_1 = NA, window_pos_2 = NA, avg_pi_A = NA, avg_pi_M = NA, avg_dxy = NA, avg_wc_fst = pcadapt$WEIR_AND_COCKERHAM_FST, no_snps = NA, GC = NA, Ngenes = NA, mid = pcadapt$POS, snp_name = pcadapt$snp_name)
# pixy
dxy <- read.table(file = "results/pixy/1Mb_window_dxy.txt", header = T)[,3:6]
pi <- read.table(file = "results/pixy/1Mb_window_pi.txt", header = T)[,1:5]
pi_A <- pi[which(pi$pop == "A"),2:5]
colnames(pi_A)[4] <- "avg_pi_A"
pi_M <- pi[which(pi$pop == "M"),2:5]
colnames(pi_M)[4] <- "avg_pi_M"
fst <- read.table(file = "results/pixy/1Mb_window_fst.txt", header = T)[,3:7]
df <- Reduce(function(x, y) merge(x, y, all.x=TRUE), list(pi_A, pi_M, dxy, fst))
# set Fst for windows with negative FST to 0
df[which(df$avg_wc_fst < 0), "avg_wc_fst"] <- 0
# gcbias
gc <- read.table("results/gc/GC_1Mb.tsv", header = F)[,c(1:3,5)]
colnames(gc) <- c("chromosome", "window_pos_1", "window_pos_2", "GC")
gc$GC <- gc$GC*100
df <- merge(df, gc, all.x = T)
# gene density
Ngenes <- read.table("results/gene_density/ngenes_1Mb.tsv", header = F) 
colnames(Ngenes) <- c("chromosome", "window_pos_1", "window_pos_2", "Ngenes")
df <- merge(df, Ngenes, all.x = T)
# save window stats
write.csv(df, file = "results/manhattan_plots/window_stats.csv", quote = F, row.names = F)
# add snp_name and mid position
df$mid <- round((df$window_pos_1 + df$window_pos_2)/2)
df$snp_name <- paste("WIND",df$chromosome, df$mid, sep = "_")
# add pcadapt results
df <- rbind(df, pcadapt)
# sort
df <- df[with(df, order(df$chromosome, df$mid)),]
# chr short name and number
chrnames <- paste0("chr",sprintf("%02d", 01:31))
names(chrnames) <- sort(unique(df$chromosome))
df$chrnames <- chrnames[df$chromosome]
df$chrnum <- as.numeric(gsub("chr", "", df$chrname))

# with gene labelling
layout(matrix(c(1,2,3,4,5,6,7), ncol = 1))
par(mar = c(1,6,0,1))
manhattan2(df, chr = "chrnum", p = "avg_wc_fst", bp = "mid", snp = "snp_name", logp=FALSE, 
           ylab=expression(italic('F'["ST"])), cex.lab = 1.5, 
          chrlabs = rep("",31), ylim = c(0,1.1), las = 2,
          xlab = "",mgp=c(3.5,1,0), xaxt = "n", 
          highlight = pcadapt$snp_name, annotateSet = pcadapt_genes$gene_name, ann_cex = 0.75, highlightcol = "lightcoral")
manhattan2(df, chr = "chrnum", p = "avg_dxy", bp = "mid", snp = "snp_name", logp=FALSE, 
          ylab=expression(italic('d'["XY"])), cex.lab = 1.5,  
          chrlabs = rep("",31), ylim = c(0,max(df$avg_dxy, na.rm = T)), 
          xlab = "",mgp=c(3.5,1,0), xaxt = "n")
manhattan2(df, chr = "chrnum", p = "avg_pi_A", bp = "mid", snp = "snp_name", logp=FALSE, 
          ylab = expression(italic(pi["Atl."])), cex.lab = 1.5,  
          chrlabs = rep("",31), ylim = c(0,max(df$avg_pi_A, na.rm = T)),
          col = c("darkblue", "cornflowerblue"),
          xlab = "",mgp=c(3.5,1,0), xaxt = "n")
manhattan2(df, chr = "chrnum", p = "avg_pi_M", bp = "mid", snp = "snp_name", logp=FALSE, 
          ylab = expression(italic(pi["Med."])), cex.lab = 1.5,  
          chrlabs = rep("",31), ylim = c(0,max(df$avg_pi_M, na.rm = T)),
          col = c("darkred", "lightcoral"),
          xlab = "",mgp=c(3.5,1,0), xaxt = "n")
manhattan2(df, chr = "chrnum", p = "Ngenes", bp = "mid", snp = "snp_name", logp=FALSE, 
          ylab="Genes per Mb", cex.lab = 1.5,
          chrlabs = rep("",31), ylim = c(0,max(df$Ngenes, na.rm = T)),
          col = c( "palegreen4", "palegreen2"),
          xlab = "",mgp=c(3.5,1,0), suggestiveline = F, genomewideline = F, xaxt = "n")
manhattan2(df, chr = "chrnum", p = "GC", bp = "mid", snp = "snp_name", logp=FALSE, 
          ylab="GC %", cex.lab = 1.5, xlab = "Chromosome",
          col =c("darkgoldenrod", "gold"),
          chrlabs = paste0("Chr",1:31), ylim = c(min(df$GC, na.rm = T),max(df$GC, na.rm = T)),
          mgp=c(3.5,1,0))

# no gene labelling
layout(matrix(c(1,2,3,4,5,6,7), ncol = 1))
par(mar = c(1,6,0,1))
manhattan2(df, chr = "chrnum", p = "avg_wc_fst", bp = "mid", snp = "snp_name", logp=FALSE, 
           ylab=expression(italic('F'["ST"])), cex.lab = 1.5, 
           chrlabs = rep("",31), ylim = c(0,1.1), las = 2,
           xlab = "",mgp=c(3.5,1,0), xaxt = "n", 
           highlight = pcadapt$snp_name, annotateSet = NULL, ann_cex = 0.75, highlightcol = "lightcoral")
manhattan2(df, chr = "chrnum", p = "avg_dxy", bp = "mid", snp = "snp_name", logp=FALSE, 
           ylab=expression(italic('d'["XY"])), cex.lab = 1.5,  
           chrlabs = rep("",31), ylim = c(0,max(df$avg_dxy, na.rm = T)), 
           xlab = "",mgp=c(3.5,1,0), xaxt = "n")
manhattan2(df, chr = "chrnum", p = "avg_pi_A", bp = "mid", snp = "snp_name", logp=FALSE, 
           ylab = expression(italic(pi["Atl."])), cex.lab = 1.5,  
           chrlabs = rep("",31), ylim = c(0,max(df$avg_pi_A, na.rm = T)),
           col = c("darkblue", "cornflowerblue"),
           xlab = "",mgp=c(3.5,1,0), xaxt = "n")
manhattan2(df, chr = "chrnum", p = "avg_pi_M", bp = "mid", snp = "snp_name", logp=FALSE, 
           ylab = expression(italic(pi["Med."])), cex.lab = 1.5,  
           chrlabs = rep("",31), ylim = c(0,max(df$avg_pi_M, na.rm = T)),
           col = c("darkred", "lightcoral"),
           xlab = "",mgp=c(3.5,1,0), xaxt = "n")
manhattan2(df, chr = "chrnum", p = "Ngenes", bp = "mid", snp = "snp_name", logp=FALSE, 
           ylab="Genes per Mb", cex.lab = 1.5,
           chrlabs = rep("",31), ylim = c(0,max(df$Ngenes, na.rm = T)),
           col = c( "palegreen4", "palegreen2"),
           xlab = "",mgp=c(3.5,1,0), suggestiveline = F, genomewideline = F, xaxt = "n")
manhattan2(df, chr = "chrnum", p = "GC", bp = "mid", snp = "snp_name", logp=FALSE, 
           ylab="GC %", cex.lab = 1.5, xlab = "Chromosome",
           col =c("darkgoldenrod", "gold"),
           chrlabs = paste0("Chr",1:31), ylim = c(min(df$GC, na.rm = T),max(df$GC, na.rm = T)),
           mgp=c(3.5,1,0))
