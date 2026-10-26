library(readxl)
library(limma)
library(ggplot2)

args <- commandArgs(TRUE)
input <- args[1]
go_file <- args[2]
ont <- "CC"       # ontology type BP, CC or MF
fc_cut <- 1
padj_cut <- 0.05
prefix <- "output/04_go"   # output file name
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

# UniProt export: Gene Names (primary) (or Gene Names) and all three GO columns.
if (grepl("\\.xlsx?$", go_file, ignore.case = TRUE)) {
  up <- as.data.frame(read_excel(go_file, col_types = "text"))
} else {
  up <- read.delim(go_file, quote = "", check.names = FALSE, colClasses = "character")
}
up[is.na(up)] <- ""
gene_col <- if ("Gene Names (primary)" %in% names(up)) "Gene Names (primary)" else "Gene Names"
if (!gene_col %in% names(up)) stop("missing Gene Names (primary)")
up$symbol <- sub("[; ].*", "", trimws(up[[gene_col]]))
go_cols <- c(BP = "biological process", CC = "cellular component", MF = "molecular function")
ann <- do.call(rbind, Map(function(x, o) {
  s <- lapply(strsplit(up[[grep(x, names(up), fixed = TRUE)]], ";", fixed = TRUE), trimws)
  data.frame(symbol = rep(up$symbol, lengths(s)), ont = rep(o, sum(lengths(s))), term = unlist(s))
}, go_cols, names(go_cols)))
ann <- ann[nzchar(ann$symbol) & ann$symbol != "0" & grepl("\\[GO:\\d+\\]$", ann$term), ]
if (!nrow(ann)) stop("No valid GO annotations.")
ann$id <- sub(".*\\[(GO:\\d+)\\]$", "\\1", ann$term)
ann$name <- sub("\\s*\\[GO:\\d+\\]$", "", ann$term)

key <- function(s) setdiff(unique(s), c("", "0", NA_character_))
term_ids <- lapply(split(ann$symbol, ann$id), key)
term_info <- ann[!duplicated(ann$id), c("id", "ont", "name")]
universe <- key(up$symbol)

if (grepl("\\.xlsx?$", input, ignore.case = TRUE)) {
  d <- as.data.frame(read_excel(input, .name_repair = "minimal"))
} else {
  d <- read.delim(input, sep = ifelse(grepl("\\.csv$", input), ",", "\t"), check.names = FALSE)
}
smp <- names(d)[grepl("\\d+$", names(d)) & sapply(d, is.numeric)]
grp <- trimws(sub("[_.-]?\\d+$", "", smp))
grp <- factor(grp, levels = unique(grp))
ref <- levels(grp)[1]

calls <- vapply(levels(grp)[-1], function(o) {
  y <- log2(as.matrix(d[, c(smp[grp == o], smp[grp == ref])]) + 1)
  f <- factor(rep(c("other", "ref"), c(sum(grp == o), sum(grp == ref))), levels = c("other", "ref"))
  design <- model.matrix(~ 0 + f)
  colnames(design) <- levels(f)
  fit <- eBayes(contrasts.fit(lmFit(y, design), makeContrasts(ref - other, levels = design)))
  padj <- p.adjust(fit$p.value[, 1], "BH")
  ifelse(!is.na(padj) & padj < padj_cut, sign(fit$coefficients[, 1]) * (abs(fit$coefficients[, 1]) > fc_cut), 0)
}, numeric(nrow(d)))
queries <- list(decreased = unique(d$PG.Genes[rowSums(calls == -1) == ncol(calls)]),
                increased = unique(d$PG.Genes[rowSums(calls == 1) == ncol(calls)]))
queries <- lapply(queries, function(s) key(trimws(as.character(s))))

N <- length(universe)
M <- lengths(term_ids)
for (q in names(queries)) {
  hits <- intersect(queries[[q]], universe)
  n <- length(hits)
  k <- vapply(term_ids, function(t) sum(hits %in% t), 0)
  keep <- M >= 5 & k > 0
  # BH
  r <- data.frame(term_info[match(names(term_ids)[keep], term_info$id), ], n_genes = k[keep], n_term = M[keep],
                  pval = phyper(k[keep] - 1, M[keep], N - M[keep], n, lower.tail = FALSE))
  r$fold_enrichment <- r$n_genes / (n * r$n_term / N)
  r$fdr <- p.adjust(r$pval, "BH")
  r$genes <- vapply(term_ids[r$id], function(t) paste(sort(intersect(hits, t)), collapse = ","), "")
  r <- r[r$ont == ont & r$fdr <= 0.05 & r$n_genes >= 3, ]
  r <- r[order(r$fdr, -r$fold_enrichment), ]
  fn <- sprintf("%s_%s_%s", prefix, ont, q)
  write.csv(r, paste0(fn, ".csv"), row.names = FALSE)
  if (!nrow(r)) next

  top <- head(r, 5)
  lab <- sprintf("%s (n = %d)", top$name, top$n_genes)
  p <- ggplot(top, aes(-log10(fdr), factor(lab, rev(lab)))) + geom_col(fill = "grey80") +
    labs(y = NULL) + theme_classic()
  ggsave(paste0(fn, ".pdf"), p, width = 5, height = 3)
  ggsave(paste0(fn, ".png"), p, width = 5, height = 3, dpi = 200)
}
