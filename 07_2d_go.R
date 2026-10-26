library(readxl)
library(ggplot2)

args <- commandArgs(TRUE)
ref <- "BrainPhys"       # numerator; all other groups are pooled as denominator
ont <- "BP"              # display ontology; BH includes BP, CC and MF
min_term_size <- 5L       # minimum scored genes inside and outside a term
fdr_cut <- 0.03           # FDR cutoff
max_terms <- 20L          # number of terms to show
prefix <- "output/07_2d_go" # output file name

read_table <- function(path) {
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    as.data.frame(read_excel(path, .name_repair = "minimal"))
  } else {
    read.delim(path, sep = if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t",
               check.names = FALSE, quote = "\"")
  }
}

make_scores <- function(d) {
  smp <- names(d)[grepl("\\d+$", names(d)) & vapply(d, is.numeric, logical(1))]
  grp <- trimws(sub("[_.-]?\\d+$", "", smp))
  m <- as.matrix(d[, smp, drop = FALSE])
  fc <- log2(rowMeans(m[, grp == ref, drop = FALSE], na.rm = TRUE) + 1) -
        log2(rowMeans(m[, grp != ref, drop = FALSE], na.rm = TRUE) + 1)
  names(fc) <- trimws(as.character(d$PG.Genes))
  valid_symbol <- !is.na(names(fc)) & nzchar(names(fc)) & names(fc) != "0" &
                  !grepl("[;[:space:]]", names(fc))
  fc <- fc[is.finite(fc) & valid_symbol]
  if (anyDuplicated(names(fc))) {
    fc <- fc[order(abs(fc), decreasing = TRUE)]
    fc <- fc[!duplicated(names(fc))]
  }
  fc
}

sx <- make_scores(read_table(args[1]))
sy <- make_scores(read_table(args[2]))
genes <- intersect(names(sx), names(sy)) 
scores <- data.frame(gene = genes, log2FC_abundance = unname(sx[genes]),
                     log2FC_half_life = unname(sy[genes]))
N <- nrow(scores)
R <- cbind(rank(scores$log2FC_abundance, ties.method = "average"),
           rank(scores$log2FC_half_life, ties.method = "average"))

up <- read_table(args[3])
gene_col <- if ("Gene Names (primary)" %in% names(up)) "Gene Names (primary)" else "Gene Names"
go_cols <- c(BP = "Gene Ontology (biological process)",
             CC = "Gene Ontology (cellular component)",
             MF = "Gene Ontology (molecular function)")
symbols <- sub("[; ].*", "", trimws(as.character(up[[gene_col]])))
ann <- do.call(rbind, lapply(names(go_cols), function(o) {
  text <- as.character(up[[go_cols[[o]]]])
  text[is.na(text)] <- ""
  items <- lapply(strsplit(text, ";", fixed = TRUE), trimws)
  data.frame(gene = rep(symbols, lengths(items)), ontology = rep(o, sum(lengths(items))),
             term = unlist(items))
}))
ann <- ann[!is.na(ann$gene) & nzchar(ann$gene) & ann$gene != "0" & grepl("\\[GO:\\d+\\]$", ann$term), ]
ann$term_id <- sub(".*\\[(GO:\\d+)\\]$", "\\1", ann$term)
ann$term_name <- sub("\\s*\\[GO:\\d+\\]$", "", ann$term)
term_genes <- lapply(split(ann$gene, ann$term_id), intersect, y = genes)
term_info <- ann[!duplicated(ann$term_id), c("term_id", "term_name", "ontology")]

results <- lapply(names(term_genes), function(id) {
  inside <- genes %in% term_genes[[id]]
  n_in <- sum(inside)
  if (n_in < min_term_size || N - n_in < min_term_size) return(NULL)
  group <- factor(ifelse(inside, "in", "out"))
  p <- summary(manova(R ~ group), test = "Wilks")$stats[1, "Pr(>F)"]
  data.frame(term_info[match(id, term_info$term_id), ], n_genes = n_in,
             score_x = 2 * (mean(R[inside, 1]) - mean(R[!inside, 1])) / N,
             score_y = 2 * (mean(R[inside, 2]) - mean(R[!inside, 2])) / N,
             mean_log2FC_abundance = mean(scores$log2FC_abundance[inside]),
             mean_log2FC_half_life = mean(scores$log2FC_half_life[inside]),
             pval = p,
             genes = paste(sort(genes[inside]), collapse = ";"))
})
results <- do.call(rbind, results)
if (is.null(results)) stop("No eligible GO terms.")
results$fdr <- p.adjust(results$pval, method = "BH")
results$neglog10_fdr <- -log10(pmax(results$fdr, 1e-300))
results <- results[order(results$fdr, -sqrt(results$score_x^2 + results$score_y^2)), ]
selected <- results[results$ontology == ont & results$fdr <= fdr_cut, ]

dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
write.csv(scores, paste0(prefix, "_input_scores.csv"), row.names = FALSE)
write.csv(results, paste0(prefix, "_all_terms.csv"), row.names = FALSE)
write.csv(selected, paste0(prefix, "_", ont, "_significant.csv"), row.names = FALSE)
if (nrow(selected)) {
  top <- head(selected, max_terms)
  p <- ggplot(top, aes(score_x, score_y)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(aes(size = n_genes, colour = neglog10_fdr)) +
    geom_text(aes(label = term_name), vjust = -0.8, size = 2.5, check_overlap = TRUE) +
    scale_colour_gradient(low = "#FEE8C8", high = "#B30000", name = "-log10(FDR)") +
    labs(x = "Abundance enrichment score", y = "Half-life enrichment score", size = "Genes") +
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) + theme_classic()
  ggsave(paste0(prefix, "_", ont, ".pdf"), p, width = 7, height = 5)
  ggsave(paste0(prefix, "_", ont, ".png"), p, width = 7, height = 5, dpi = 200)
}
