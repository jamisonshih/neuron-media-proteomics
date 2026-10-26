library(readxl)
library(ggplot2)

input <- commandArgs(TRUE)[1]
prefix <- "output/01_pca"   # output file name
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

if (grepl("\\.xlsx?$", input, ignore.case = TRUE)) {
  d <- as.data.frame(read_excel(input, .name_repair = "minimal"))
} else {
  d <- read.delim(input, sep = ifelse(grepl("\\.csv$", input), ",", "\t"), check.names = FALSE)
}
# sample columns end in a replicate number (Ctrl1, Ctrl2, etc.)
smp <- names(d)[grepl("\\d+$", names(d)) & sapply(d, is.numeric)]
grp <- trimws(sub("[_.-]?\\d+$", "", smp))
grp <- factor(grp, levels = unique(grp))

m <- log10(as.matrix(d[, smp]))
m <- m[complete.cases(m), ]
pca <- prcomp(t(m), scale. = TRUE)

# biggest loading on each PC is positive
pc <- sweep(pca$x, 2, apply(pca$rotation, 2, function(v) sign(v[which.max(abs(v))])), "*")
ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
cat(sprintf("PC%d %.1f%%\n", 1:4, ve[1:4]), sep = "")

sc <- data.frame(sample = smp, group = grp, PC1 = pc[, 1], PC2 = pc[, 2], PC3 = pc[, 3])
write.csv(sc, paste0(prefix, "_scores.csv"), row.names = FALSE)
write.csv(data.frame(PC = paste0("PC", seq_along(ve)), variance_pct = ve), paste0(prefix, "_variance.csv"), row.names = FALSE)

# ellipses
pad <- 0.0265 * max(diff(range(sc$PC1)), diff(range(sc$PC2)))
th <- seq(0, 2 * pi, length.out = 100)
ell <- do.call(rbind, lapply(split(sc, sc$group), function(g) {
  cx <- mean(g$PC1); cy <- mean(g$PC2); dx <- g$PC1 - cx; dy <- g$PC2 - cy
  ev <- eigen(matrix(c(sum(dx^2), sum(dx * dy), sum(dx * dy), sum(dy^2)), 2))$vectors
  ang <- atan2(ev[2, 1], ev[1, 1])
  a <- max(abs(dx * cos(ang) + dy * sin(ang))) + pad
  b <- max(max(abs(-dx * sin(ang) + dy * cos(ang))) + pad, 0.15 * a)
  data.frame(group = g$group[1],
             x = cx + a * cos(th) * cos(ang) - b * sin(th) * sin(ang),
             y = cy + a * cos(th) * sin(ang) + b * sin(th) * cos(ang))
}))

pal <- c("#FAAFAA", "#B8E0B1", "#BFE8F5", "#ACB1DA")   # group colors
p <- ggplot(sc, aes(PC1, PC2, fill = group)) +
  geom_polygon(data = ell, aes(x, y), alpha = 0.6) +
  geom_point(size = 3, shape = 21) +
  scale_fill_manual(values = pal) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) + theme_classic()
ggsave(paste0(prefix, ".pdf"), p, width = 5, height = 4)
ggsave(paste0(prefix, ".png"), p, width = 5, height = 4, dpi = 200)
