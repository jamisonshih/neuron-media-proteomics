library(readxl)

input <- commandArgs(TRUE)[1]
log_base <- 10        
distance <- "pearson"     # "pearson"  or "euclidean" 
linkage <- "complete"     # "complete" or "ward.D2"
prefix <- "output/02"     # output file name 
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

if (grepl("\\.xlsx?$", input, ignore.case = TRUE)) {
  d <- as.data.frame(read_excel(input, .name_repair = "minimal"))
} else {
  d <- read.delim(input, sep = ifelse(grepl("\\.csv$", input), ",", "\t"), check.names = FALSE)
}
smp <- names(d)[grepl("\\d+$", names(d)) & sapply(d, is.numeric)]
grp <- trimws(sub("[_.-]?\\d+$", "", smp))
grp <- factor(grp, levels = unique(grp))
m <- log(as.matrix(d[, smp]), log_base)

# spearman 
gm <- sapply(levels(grp), function(g) rowMeans(m[, grp == g, drop = FALSE], na.rm = TRUE))
gm[!is.finite(gm)] <- NA
rho <- cor(gm, method = "spearman", use = "pairwise.complete.obs")
write.csv(round(rho, 4), paste0(prefix, "_group_spearman.csv"))

pdf(paste0(prefix, "_group_correlation.pdf"), width = 7, height = 7)
pairs(gm, upper.panel = NULL, lower.panel = function(x, y) {
  points(x, y, pch = 16, cex = 0.3, col = densCols(x, y, colramp = colorRampPalette(hcl.colors(9, "Inferno"))))
  legend("topleft", legend = bquote(rho == .(sprintf("%.3f", cor(x, y, method = "spearman", use = "complete.obs")))), bty = "n")
})
invisible(dev.off())

# replicates
mc <- m[complete.cases(m), ]
r <- cor(mc, method = "pearson")
write.csv(round(r, 4), paste0(prefix, "_sample_pearson.csv"))
write.csv(round(cor(mc, method = "spearman"), 4), paste0(prefix, "_sample_spearman.csv"))
hc <- hclust(if (distance == "pearson") as.dist(1 - r) else dist(t(mc)), method = linkage)

pal <- c("#F05133", "#75C043", "#46A5DC", "#A674B2")   # label colors
dend <- dendrapply(as.dendrogram(hc), function(n) {
  if (is.leaf(n)) attr(n, "nodePar") <- list(pch = NA, lab.col = pal[grp[match(attr(n, "label"), smp)]])
  n
})
for (f in c("pdf", "png")) {
  if (f == "pdf") pdf(paste0(prefix, "_replicate_dendrogram.pdf"), width = 5, height = 5)
  else png(paste0(prefix, "_replicate_dendrogram.png"), width = 1000, height = 1000, res = 200)
  par(mar = c(2, 1, 1, 10))
  plot(dend, horiz = TRUE)
  invisible(dev.off())
}
