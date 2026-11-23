# Figure 1: (a) total effect of PCSK9 on CAD and its indirect and direct effects from the two-step
# method and MVMR; (b) effect of LDL-C on CAD. Main analysis: LDL-C instruments in the APOE and LPA
# regions excluded, LD among instruments within 10 Mb of PCSK9 modelled.
# Input: results/estimates/estimates.tsv (scripts/estimate_effects.py). Output: results/figures/.
# Run from the repository root: Rscript scripts/figure1.R
suppressMessages(library(ggplot2))

root <- Sys.getenv("PCSK9_ROOT", getwd())
results <- Sys.getenv("PCSK9_RESULTS", file.path(root, "results"))
out_dir <- file.path(results, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

est <- read.delim(file.path(results, "estimates", "estimates.tsv"), stringsAsFactors = FALSE)
est <- est[est$analysis == "main" & est$exclusion == "APOE+LPA excluded", ]
stopifnot(nrow(est) == 2)
ci <- function(setting, quantity, method, e, s) {
  data.frame(setting = setting, quantity = quantity, method = method, estimate = e,
             lower = e - 1.96 * s, upper = e + 1.96 * s)
}
rows <- list()
for (i in seq_len(nrow(est))) {
  r <- est[i, ]
  rows[[length(rows) + 1]] <- rbind(
    ci(r$setting, "Total", "Total", r$total, r$total_se),
    ci(r$setting, "Indirect", "Two-step", r$indirect_twostep, r$indirect_twostep_se),
    ci(r$setting, "Indirect", "MVMR", r$indirect_mvmr, r$indirect_mvmr_se),
    ci(r$setting, "Direct", "Two-step", r$direct_twostep, r$direct_twostep_se),
    ci(r$setting, "Direct", "MVMR", r$direct_mvmr, r$direct_mvmr_se),
    ci(r$setting, "LDL", "Two-step", r$ldl_cad_uni, r$ldl_cad_uni_se),
    ci(r$setting, "LDL", "MVMR", r$ldl_cad_mvmr, r$ldl_cad_mvmr_se))
}
d <- do.call(rbind, rows)

ink <- c("Total" = "black", "Two-step" = "#2166ac", "MVMR" = "#aa0202")
theme_fp <- theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    axis.text = element_text(colour = "black", size = 11.5),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11.5),
    panel.spacing = unit(1.4, "lines"),
    plot.tag = element_text(face = "bold", size = 15)
  )

# Panel a. Total effect of PCSK9 on CAD and its indirect and direct effects.
dec <- d[d$quantity != "LDL", ]
dec$y <- c("Total" = 3, "Indirect" = 2, "Direct" = 1)[dec$quantity] +
  ifelse(dec$method == "Two-step", 0.17, ifelse(dec$method == "MVMR", -0.17, 0))
pa <- ggplot(dec, aes(x = estimate, y = y, colour = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.12, linewidth = 0.7, orientation = "y") +
  geom_point(shape = 15, size = 3.2) +
  facet_wrap(~setting, nrow = 1) +
  scale_colour_manual(values = ink, guide = "none") +
  scale_y_continuous(breaks = 3:1, labels = c("Total", "Indirect", "Direct"), limits = c(0.6, 3.4)) +
  scale_x_continuous(limits = c(-0.05, 0.3), breaks = seq(0, 0.3, 0.1)) +
  labs(x = "Log odds ratio for CAD per SD higher PCSK9", tag = "a") +
  theme_fp

# Panel b. Effect of LDL-C on CAD. The legend sits here, at the foot of the figure.
ldl <- d[d$quantity == "LDL", ]
ldl$row <- factor(ifelse(ldl$method == "MVMR", "MVMR,\nadjusted for PCSK9", "Univariable IVW,\nfor the two-step method"),
                  levels = c("MVMR,\nadjusted for PCSK9", "Univariable IVW,\nfor the two-step method"))
# an empty "Total" row so the legend lists all three keys, as in panel a
ldl <- rbind(ldl, transform(ldl[1, ], method = "Total", estimate = NA, lower = NA, upper = NA))
pb <- ggplot(ldl, aes(x = estimate, y = row, colour = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.18, linewidth = 0.7, orientation = "y", na.rm = TRUE) +
  geom_point(shape = 15, size = 3.2, na.rm = TRUE) +
  facet_wrap(~setting, nrow = 1) +
  scale_colour_manual(values = ink, breaks = names(ink), limits = names(ink)) +
  scale_x_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  labs(x = "Log odds ratio for CAD per SD higher LDL-C", tag = "b") +
  theme_fp

ga <- ggplotGrob(pa)
gb <- ggplotGrob(pb)
w <- grid::unit.pmax(ga$widths, gb$widths)
ga$widths <- w
gb$widths <- w
draw <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = grid::unit(c(1.45, 1.2), "null"))))
  grid::pushViewport(grid::viewport(layout.pos.row = 1)); grid::grid.draw(ga); grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 2)); grid::grid.draw(gb); grid::upViewport()
}
png(file.path(out_dir, "Figure1.png"), width = 9.6, height = 8.2, units = "in", res = 300, bg = "white")
draw(); invisible(dev.off())
pdf(file.path(out_dir, "Figure1.pdf"), width = 9.6, height = 8.2, bg = "white")
draw(); invisible(dev.off())
cat("wrote", out_dir, "\n")
