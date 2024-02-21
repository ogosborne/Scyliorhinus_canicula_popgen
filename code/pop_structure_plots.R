###### PCA
# load pca
pca <- read.table("results/pca/pca.eigenvec", header = FALSE)
eigenval <- scan("results/pca/pca.eigenval")
# remove useless col
pca <- pca[,-1]
# set col names
colnames(pca) <- c("ind",paste0("PC", 1:(ncol(pca)-1)))
# load metadata
metadata <- read.csv("data/sample_sheet.csv")
rownames(metadata) <- metadata$sample
# change subpop names
subpop_names <- c(Brittany = "English Channel", North_Wales = "N. Wales", Med = "Mediterranean")
subpop <- subpop_names[metadata[pca$ind, "subpop"]]
# remake pca
pca <- data.frame(pca, subpop)
# convert eigenvalues to percentage variance explained
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

###### ADMIXTURE
# crossvalidation
cv <- read.table("results/admixture/cv.error")
colnames(cv) <- c("K", "CVE")

# K2
K2 <- as.matrix(read.table("results/admixture/thin_100kb.2.Q"))
#rownames(K2) <- names
K2 <- t(K2)

###### Plots
# combined plot
layout(matrix(c(1,3,3,
                2,3,3), 
              byrow = T,
              ncol = 3
              ))
par(mar = c(6,8,5,0), pty="m", xpd = T)
# plot admixture plot
barplot(K2, las = 2, col = c("red","blue"), main = "",  cex.axis = 1.3,  cex.lab = 1.5, names.arg = rep("", 14), ylab = "Ancestry")
lines(x = c(0.2,2*1.2+1.2), y = c(-0.03, -0.03), col = "darkblue", lwd = 4)
lines(x = c(3.8,4*1.2+1.2), y = c(-0.03, -0.03), col = "cornflowerblue", lwd = 4)
lines(x = c(6.2,13*1.2+1.2), y = c(-0.03, -0.03), col = "red", lwd = 4)
lines(x = c(0.2,2*1.2+1.2), y = c(1.03, 1.03), col = "darkblue", lwd = 4)
lines(x = c(3.8,4*1.2+1.2), y = c(1.03, 1.03), col = "cornflowerblue", lwd = 4)
lines(x = c(6.2,13*1.2+1.2), y = c(1.03, 1.03), col = "red", lwd = 4)
mtext("a", side = 2, line = 4, at = 1, font = 2, las = 2, cex = 3)
# plot ADMIXTURE cross-validation error
par(mar = c(5,8,0,0))
plot(cv$K, cv$CVE, type = "b", xlab = "K", ylab = "CV error",  cex.axis = 1.3, cex.lab = 1.5, bty = "n", las = 1,  mgp = c(3,1,0))
mtext("b", side = 2, line = 4, at = 3, font = 2, las = 2, cex = 3)
# plot PCA
col = c("Mediterranean" = "red", "N. Wales" = "darkblue", "English Channel" = "cornflowerblue")
par(mar = c(5,7,2,3))
plot(pca$PC1, pca$PC2, col = col[pca$subpop],  
     xlim = c(-0.4, 0.55), ylim = c(-0.65, 0.65), 
     xlab = paste0("PC1 (", signif(pve$pve[1], 3), "%)"),
     ylab = paste0("PC2 (", signif(pve$pve[2], 3), "%)"),
     cex = 4, pch = 19, cex.axis = 1.3, cex.lab = 1.5
     )
# legend
legend("bottomleft", legend = names(col), col = col, pch = 19, cex = 1.5, bg = NA, bty = "n")
mtext("c", side = 4, line = 0.6, at = 0.7, font = 2, las = 2, cex = 3)

#Other plots
par(mfrow=c(1,1))
# PCA 3 - 4
plot(pca$PC3, pca$PC4, col = col[pca$subpop],  
     xlab = paste0("PC3 (", signif(pve$pve[3], 3), "%)"),
     ylab = paste0("PC4 (", signif(pve$pve[4], 3), "%)"),
     cex = 2, pch = 19, cex.axis = 1.3, cex.lab = 1.5, pty = "s"
)
# PCA PVE
plot(pve$PC, pve$pve, type = "b", xlab = "PC", ylab = "Percent variance explained")
