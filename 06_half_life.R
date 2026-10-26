library(readxl)

args <- commandArgs(TRUE)
input <- args[1]
td <- if (length(args) > 1) as.numeric(args[2]) else 4          # labeling time 
lower <- if (length(args) > 2) as.numeric(args[3]) else 0.02    # H/L ratio cutoffs
upper <- if (length(args) > 3) as.numeric(args[4]) else 100
protein_col <- "PG.ProteinGroups"      
gene_col <- "PG.Genes"
peptide_col <- "PEP.StrippedSequence"
prefix <- "output/06"   # output file name 
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

if (grepl("\\.xlsx?$", input, ignore.case = TRUE)) {
  d <- as.data.frame(read_excel(input, .name_repair = "minimal"))
} else {
  d <- read.delim(input, sep = ifelse(grepl("\\.csv$", input), ",", "\t"), check.names = FALSE)
}
smp <- names(d)[grepl("\\d+$", names(d)) & sapply(d, is.numeric)]
hl <- as.matrix(d[, smp])

# ratio filter
hl[!is.finite(hl) | hl <= 0 | hl < lower | hl > upper] <- NA

t_half <- td * log(2) / log(1 + hl)
pep <- data.frame(d[, c(protein_col, gene_col, peptide_col)], t_half, check.names = FALSE)
write.csv(pep, paste0(prefix, "_peptide_half_life.csv"), row.names = FALSE)

hmean <- function(x) { x <- x[!is.na(x)]; if (length(x)) length(x) / sum(1 / x) else NA }
key <- paste(d[[protein_col]], d[[gene_col]])
f <- factor(key, levels = unique(key))
prot <- apply(t_half, 2, function(v) tapply(v, f, hmean))
out <- data.frame(d[!duplicated(key), c(protein_col, gene_col)], n_peptides = as.vector(table(f)), prot,
                  check.names = FALSE)
out <- out[rowSums(!is.na(prot)) > 0, ]
write.csv(out, paste0(prefix, "_protein_half_life.csv"), row.names = FALSE)
