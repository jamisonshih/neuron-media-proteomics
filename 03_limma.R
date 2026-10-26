library(readxl)
library(limma)
library(ggplot2)

input <- commandArgs(TRUE)[1]
fc_cut <- 1       # |log2FC|
padj_cut <- 0.05
prefix <- "output/03"   # output file name
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

if (grepl("\\.xlsx?$", input, ignore.case = TRUE)) {
  d <- as.data.frame(read_excel(input, .name_repair = "minimal"))
} else {
  d <- read.delim(input, sep = ifelse(grepl("\\.csv$", input), ",", "\t"), check.names = FALSE)
}
smp <- names(d)[grepl("\\d+$", names(d)) & sapply(d, is.numeric)]
grp <- trimws(sub("[_.-]?\\d+$", "", smp))
grp <- factor(grp, levels = unique(grp))
pal <- setNames(c("#FAAFAA", "#B8E0B1", "#BFE8F5", "#ACB1DA"), levels(grp))   # group color

res <- list()
for (k in combn(levels(grp), 2, simplify = FALSE)) {
  # comparison groups
  num <- smp[grp == k[1]]; den <- smp[grp == k[2]]
  y <- log2(as.matrix(d[, c(den, num)]) + 1)
  f <- factor(rep(c("den", "num"), c(length(den), length(num))))
  design <- model.matrix(~ 0 + f)
  colnames(design) <- levels(f)
  fit <- eBayes(contrasts.fit(lmFit(y, design), makeContrasts(num - den, levels = design)))
  cmp <- paste0(k[1], "/", k[2])
  r <- data.frame(gene = d$PG.Genes, comparison = cmp, log2FC = fit$coefficients[, 1], p = fit$p.value[, 1])
  r$padj <- p.adjust(r$p, "BH")
  r$sig <- ifelse(!is.na(r$padj) & r$padj < padj_cut & r$log2FC > fc_cut, "up",
                  ifelse(!is.na(r$padj) & r$padj < padj_cut & r$log2FC < -fc_cut, "down", "ns"))
  cat(sprintf("%-40s up %4d down %4d\n", cmp, sum(r$sig == "up"), sum(r$sig == "down")))
  res[[cmp]] <- r

  p <- ggplot(r, aes(log2FC, -log10(padj), colour = sig)) + geom_point(size = 1, na.rm = TRUE) +
    scale_colour_manual(values = c(up = pal[[k[1]]], down = pal[[k[2]]], ns = "grey87"), guide = "none") +
    labs(x = sprintf("log2(%s / %s)", k[1], k[2])) + theme_classic()
  fn <- paste0(prefix, "_volcano_", gsub("[^A-Za-z0-9]+", "_", cmp))
  ggsave(paste0(fn, ".pdf"), p, width = 4, height = 4)
  ggsave(paste0(fn, ".png"), p, width = 4, height = 4, dpi = 200)
}
all <- do.call(rbind, res)
write.csv(all, paste0(prefix, "_limma_results.csv"), row.names = FALSE)

ref <- levels(grp)[1]
rr <- all[startsWith(all$comparison, paste0(ref, "/")), ]
cmp <- factor(rr$comparison)  
for (dir in c("down", "up")) {
  hits <- Reduce(intersect, lapply(split(rr$gene[rr$sig == dir], cmp[rr$sig == dir]), unique))
  cat(dir, "in", ref, "consistent:", length(hits), "\n")
  writeLines(sort(hits), sprintf("%s_consistent_%s_in_%s.txt", prefix, dir, gsub("[^A-Za-z0-9]+", "_", ref)))
}
