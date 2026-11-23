# Supplementary Figure S1: correlations (r) between the 16 PCSK9 instruments and the Setting 2 LDL-C
# instruments within 10 Mb of PCSK9 (24 variants, 1000 Genomes EUR), the matrix used in the estimation.
# Inputs: summary_results/. Output: results/figures/FigureS1.png and FigureS1.pdf.
# Run from the repository root: Rscript scripts/figureS1.R
suppressMessages(library(ggplot2))

root <- Sys.getenv("PCSK9_ROOT", getwd())
summ <- file.path(root, "summary_results")
out_dir <- file.path(Sys.getenv("PCSK9_RESULTS", file.path(root, "results")), "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gene_start <- 55505221
cis <- read.delim(file.path(summ, "pcsk9_instruments.tsv"), stringsAsFactors = FALSE)
cis_ids <- cis$rsid
s2 <- read.delim(file.path(summ, "ldl_setting2.tsv"), stringsAsFactors = FALSE)
ldl_ids <- s2$rsid[s2$chromosome == 1 & abs(s2$position - gene_start) <= 10e6]
bim <- read.table(file.path(summ, "region_snps.bim"), sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("chr", "rsid", "cm", "bp", "a1", "a2"))
stopifnot(setequal(bim$rsid, union(cis_ids, ldl_ids)))
# PLINK reports r for the A1 allele; A1 is the LDL-C effect allele for all 24 variants.
eff <- c(setNames(cis$effect_allele, cis$rsid), setNames(s2$effect_allele, s2$rsid))
stopifnot(all(toupper(eff[bim$rsid]) == toupper(bim$a1)))
mat <- as.matrix(read.table(file.path(summ, "region_corr.ld"), sep = "\t"))
mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
n <- nrow(mat)
labels <- ifelse(bim$rsid %in% ldl_ids, paste0(bim$rsid, "*"), bim$rsid)
grid <- expand.grid(i = seq_len(n), j = seq_len(n))
grid$r <- mat[cbind(grid$i, grid$j)]
grid$y <- n + 1 - grid$i
grid$label <- sprintf("%.2f", round(grid$r, 2) + 0)
p3 <- ggplot(grid, aes(x = j, y = y, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = label, colour = abs(r) >= 0.65), size = 3.4, show.legend = FALSE) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#aa0202", midpoint = 0, limits = c(-1, 1), name = "r") +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
  scale_x_continuous(breaks = seq_len(n), labels = labels, expand = c(0, 0)) +
  scale_y_continuous(breaks = n:1, labels = labels, expand = c(0, 0)) +
  coord_fixed() +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 11, colour = "black", angle = 55, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.title = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 12), plot.margin = margin(10, 14, 10, 10))
ggsave(file.path(out_dir, "FigureS1.png"), p3, width = 12.4, height = 10.6, units = "in", dpi = 220, bg = "white")
ggsave(file.path(out_dir, "FigureS1.pdf"), p3, width = 12.4, height = 10.6, units = "in", device = "pdf", bg = "white")
cat("wrote", file.path(out_dir, "FigureS1.png"), "and FigureS1.pdf\n")
