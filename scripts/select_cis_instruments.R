# Cis PCSK9 instruments from the three-study window overlap.
# Requires results/window/pcsk9_window_overlap.tsv and the 1000 Genomes
# correlation matrix results/window/eur_window.ld described in data/README.md.

root <- Sys.getenv("PCSK9_RESULTS", unset = "results")
ov_path <- Sys.getenv("PCSK9_OVERLAP", unset = file.path(root, "window", "pcsk9_window_overlap.tsv"))
bim_path <- Sys.getenv("PCSK9_BIM", unset = file.path(root, "window", "eur_window.bim"))
ld_path <- Sys.getenv("PCSK9_LD", unset = file.path(root, "window", "eur_window.ld"))
out_path <- Sys.getenv("PCSK9_CIS_OUT", unset = file.path(root, "pcsk9_instruments.tsv"))

comp_map <- c(A = "T", T = "A", C = "G", G = "C")
harmonize_beta <- function(beta, ea, oa, ref_ea, ref_oa) {
  ea <- toupper(ea)
  oa <- toupper(oa)
  ref_ea <- toupper(ref_ea)
  ref_oa <- toupper(ref_oa)
  beta <- as.numeric(beta)
  out <- rep(NA_real_, length(beta))
  known <- !is.na(ea) & !is.na(oa) & !is.na(ref_ea) & !is.na(ref_oa)
  direct <- known & ea == ref_ea & oa == ref_oa
  swapped <- known & ea == ref_oa & oa == ref_ea
  out[direct] <- beta[direct]
  out[swapped] <- -beta[swapped]
  rest <- !(direct | swapped) & is.finite(beta)
  if (any(rest)) {
    cea <- unname(comp_map[ea])
    coa <- unname(comp_map[oa])
    pal <- (ea %in% c("A", "T") & oa %in% c("A", "T")) |
      (ea %in% c("C", "G") & oa %in% c("C", "G"))
    ok <- rest & !pal & !is.na(cea) & !is.na(coa)
    out[ok & cea == ref_ea & coa == ref_oa] <- beta[ok & cea == ref_ea & coa == ref_oa]
    out[ok & cea == ref_oa & coa == ref_ea] <- -beta[ok & cea == ref_oa & coa == ref_ea]
  }
  out
}

raw <- read.delim(ov_path, stringsAsFactors = FALSE)
bim <- read.table(bim_path, stringsAsFactors = FALSE)
snps <- bim$V2
ld <- as.matrix(read.table(ld_path))
stopifnot(nrow(ld) == length(snps), ncol(ld) == length(snps))
rownames(ld) <- snps
colnames(ld) <- snps

df <- raw[raw$rsid %in% snps, ]
df$beta_CAD <- harmonize_beta(
  df$beta_CAD, df$effect_allele_cad, df$other_allele_cad,
  df$effect_allele_ldl, df$other_allele_ldl
)
df$beta_PCSK9_pQTL <- harmonize_beta(
  df$beta_PCSK9_pQTL, df$effect_allele_pqtl, df$other_allele_pqtl,
  df$effect_allele_ldl, df$other_allele_ldl
)
ok <- is.finite(df$beta_LDL_C) & is.finite(df$se_LDL_C) & df$se_LDL_C > 0 &
  is.finite(df$beta_CAD) & is.finite(df$se_CAD) & df$se_CAD > 0 &
  is.finite(df$beta_PCSK9_pQTL) & is.finite(df$se_PCSK9_pQTL) & df$se_PCSK9_pQTL > 0
df <- df[ok, ]
rownames(df) <- df$rsid

pval_of <- function(p, b, se) {
  p <- as.numeric(p)
  z <- abs(as.numeric(b) / as.numeric(se))
  if (all(!is.finite(p))) p <- 2 * pnorm(-z)
  p[!is.finite(p) | p == 0] <- 2 * pnorm(-z[!is.finite(p) | p == 0])
  p
}
clump_ids <- function(ids, pvals, z, r2_thresh, p_thresh) {
  ok <- is.finite(pvals) & pvals < p_thresh & is.finite(z)
  ord <- order(pvals[ok], -z[ok])
  cand <- ids[ok][ord]
  kept <- character(0)
  for (s in cand) {
    if (length(kept) == 0) {
      kept <- s
      next
    }
    r2 <- as.numeric(ld[s, kept])^2
    if (all(is.finite(r2)) && all(r2 < r2_thresh)) kept <- c(kept, s)
  }
  kept
}

df$p_pq <- pval_of(df$pval_PCSK9_pQTL, df$beta_PCSK9_pQTL, df$se_PCSK9_pQTL)
df$z_pq <- abs(df$beta_PCSK9_pQTL / df$se_PCSK9_pQTL)
pq_ids <- clump_ids(df$rsid, df$p_pq, df$z_pq, 0.1, 5e-8)
pq <- df[pq_ids, ]
out <- data.frame(
  rsid = pq$rsid,
  position_b37 = pq$position_b37,
  effect_allele = pq$effect_allele_ldl,
  other_allele = pq$other_allele_ldl,
  eaf_ldl = pq$maf_LDL_C,
  beta_pqtl = pq$beta_PCSK9_pQTL,
  se_pqtl = pq$se_PCSK9_pQTL,
  p_pqtl = pq$pval_PCSK9_pQTL,
  beta_ldl = pq$beta_LDL_C,
  se_ldl = pq$se_LDL_C,
  p_ldl = pq$pval_LDL_C,
  beta_cad = pq$beta_CAD,
  se_cad = pq$se_CAD,
  n_ldl = pq$n_LDL_C,
  n_pqtl = pq$n_PCSK9_pQTL,
  stringsAsFactors = FALSE
)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("cis instruments ", nrow(out), " written to ", out_path)
