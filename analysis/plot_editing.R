# Net A-to-I editing per construct, read from data/editing_by_rep.csv.
# No editing value is hard-coded here; significance comes from
# results/stats_table.csv (BH-adjusted). Bars are mean, error bars SEM,
# points are individual biological replicates (fixed jitter seed). Net values
# can be negative (e.g. ACTB scramble) and are shown as-is, not truncated.
# Outputs vector PDF + 600-dpi PNG to results/figures/.
suppressMessages(library(ggplot2))
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
set.seed(1)

d  <- read.csv("data/editing_by_rep.csv", stringsAsFactors = FALSE)
d$net_pct <- suppressWarnings(as.numeric(d$net_pct))
d  <- d[!is.na(d$net_pct), ]
st <- read.csv("results/stats_table.csv", stringsAsFactors = FALSE)

lab_from_padj <- function(p, n) {
  if (is.na(p)) return(sprintf("n=%d", n))
  if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else "ns"
}

dose_levels <- c("scramble", "0.3 uM", "1 uM")
dose_pretty <- c("scramble" = "Scramble", "0.3 uM" = "0.3 \u00b5M", "1 uM" = "1 \u00b5M")

plot_one <- function(con, fill) {
  x <- d[d$construct == con, ]
  x$dose <- factor(x$dose, levels = dose_levels[dose_levels %in% x$dose])
  agg <- do.call(rbind, lapply(split(x, x$dose), function(g) data.frame(
    dose = g$dose[1], mean = mean(g$net_pct),
    se = sd(g$net_pct)/sqrt(nrow(g)), n = nrow(g))))
  agg$lab <- mapply(function(dz, n) {
    row <- st[st$endpoint == sprintf("%s editing, %s", con, dz), ]
    lab_from_padj(if (nrow(row)) row$p_adj else NA, n)
  }, as.character(agg$dose), agg$n)

  yr <- range(c(0, x$net_pct, agg$mean + agg$se)); pad <- diff(yr) * 0.12
  p <- ggplot(agg, aes(dose, mean)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
    geom_col(fill = fill, colour = "black", width = 0.6) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.16) +
    geom_point(data = x, aes(dose, net_pct), shape = 21, fill = "white",
               colour = "black", size = 2,
               position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    geom_text(aes(y = mean + se + pad*0.5, label = lab), size = 4) +
    scale_x_discrete(labels = dose_pretty[levels(x$dose)]) +
    scale_y_continuous(limits = c(min(yr[1], 0), yr[2] + pad)) +
    labs(y = "A-to-I editing, net (%)", x = NULL,
         subtitle = sprintf("%s  (mean \u00b1 SEM, n per bar; BH-adjusted)", con)) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.subtitle = element_text(size = 10),
          text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
          axis.title = element_text(colour = "black"), axis.ticks = element_line(colour = "black"))
  ggsave(sprintf("results/figures/%s_editing.pdf", con), p, width = 3.8, height = 3.9)
  ggsave(sprintf("results/figures/%s_editing.png", con), p, width = 3.8, height = 3.9, dpi = 600)
}

plot_one("A74", "grey78")
plot_one("A88", "grey45")
plot_one("ACTB", "grey60")
cat("wrote A74/A88/ACTB editing figures (PDF + PNG) to results/figures/\n")
