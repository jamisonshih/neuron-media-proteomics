library(readxl)
library(pheatmap)

args <- commandArgs(TRUE)
input <- args[1]
genes <- readLines(args[2])   # gene list
prefix <- "output/05_heatmap"   # output file name
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

if (grepl("\\.xlsx?$", input, ignore.case = TRUE)) {
  d <- as.data.frame(read_excel(input, .name_repair = "minimal"))
} else {
  d <- read.delim(input, sep = ifelse(grepl("\\.csv$", input), ",", "\t"), check.names = FALSE)
}
smp <- names(d)[grepl("\\d+$", names(d)) & sapply(d, is.numeric)]
grp <- trimws(sub("[_.-]?\\d+$", "", smp))
grp <- factor(grp, levels = unique(grp))

x <- log10(as.matrix(d[match(genes, d$PG.Genes), smp]))
rownames(x) <- genes

z <- t(apply(x, 1, function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)))
write.csv(round(z, 4), paste0(prefix, "_zscore.csv"))
for (f in c("pdf", "png")) {
  pheatmap(t(z), cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = head(cumsum(table(grp)), -1),
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100), breaks = seq(-2, 2, length.out = 101),
           annotation_row = data.frame(Group = grp, row.names = smp),
           annotation_colors = list(Group = setNames(c("#FAAFAA", "#B8E0B1", "#BFE8F5", "#ACB1DA"), levels(grp))),
           show_rownames = FALSE, filename = paste0(prefix, ".", f), width = 8, height = 3.5)
}
