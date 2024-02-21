# this is a slightly modified version of the manhattan function from Stephen Turner's qqman package, adding the ability to set highlight and line colours and label an arbitrary set of SNPs

manhattan2 <- function(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", "gray60"), highlightcol = "green", suggestiveline_col = "blue", genomewideline_col = "red", chrlabs = NULL, suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = FALSE, annotateSet = NULL, ann_cex = 0.45, ...) 
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  if (!is.null(x[[snp]])) 
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                   pos = NA, index = NA, SNP = x[[snp]], stringsAsFactors = FALSE)
  else d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                      pos = NA, index = NA)
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, 
                                                             d$CHR, length))
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    xlabel = paste("Chromosome", unique(d$CHR), "position")
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + max(d[d$index == (i - 
                                                  1), "BP"])
        d[d$index == i, "BP"] = d[d$index == i, "BP"] - 
          min(d[d$index == i, "BP"]) + 1
        d[d$index == i, "pos"] = d[d$index == i, "BP"] + 
          lastbase
      }
    }
    ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep_len(col, max(d$index))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      points(d[d$index == i, "pos"], d[d$index == i, "logp"], 
             col = col[icol], pch = 20, ...)
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = suggestiveline_col)
  if (genomewideline) 
    abline(h = genomewideline, col = genomewideline_col)
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = highlightcol, 
                             pch = 20, ...))
  }
  if(!is.null(annotateSet)){
    topHits = subset(d, SNP %in% annotateSet)
    if (logp) {
      with(subset(d, SNP %in% annotateSet), calibrate::textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = ann_cex), ...)
    }
    else with(subset(d, SNP %in% annotateSet), calibrate::textxy(pos, P, offset = 0.625, labs = topHits$SNP, cex = ann_cex), ...)
  }
  if (!is.null(annotatePval)) {
    if (logp) {
      topHits = subset(d, P <= annotatePval)
    }
    else topHits = subset(d, P >= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      if (logp) {
        with(subset(d, P <= annotatePval), calibrate::textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = ann_cex), ...)
      }
      else with(subset(d, P >= annotatePval), calibrate::textxy(pos, P, offset = 0.625, labs = topHits$SNP, cex = ann_cex), ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      if (logp) {
        calibrate::textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
               labs = topSNPs$SNP, cex = ann_cex, ...)
      }
      else calibrate::textxy(topSNPs$pos, topSNPs$P, offset = 0.625, 
                  labs = topSNPs$SNP, cex = ann_cex, ...)
    }
  }
  par(xpd = FALSE)
}
