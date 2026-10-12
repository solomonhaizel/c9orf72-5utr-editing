# Statistics for single-guide A-to-I editing.
#
# Each net value = guide_raw - matched same-run no-guide, i.e. a paired
# difference. A one-sample t-test of these differences against 0 is therefore
# identical to a paired t-test of guide vs its matched no-guide control; the
# table labels it accordingly. Groups with n < 3 are reported descriptively
# (no test). P values are adjusted across the family of guide-vs-control tests
# by Benjamini-Hochberg (FDR). Effect size is Cohen's dz (mean/SD of the
# paired differences); the 95% CI is on the mean net editing.
#
# Writes results/stats_table.csv.

d <- read.csv("data/editing_by_rep.csv", stringsAsFactors = FALSE)
d$net_pct <- suppressWarnings(as.numeric(d$net_pct))
d <- d[!is.na(d$net_pct), ]

rows <- list()
for (con in unique(d$construct)) {
  for (dose in unique(d$dose[d$construct == con])) {
    v <- d$net_pct[d$construct == con & d$dose == dose]
    n <- length(v); m <- mean(v); s <- sd(v)
    endpoint   <- sprintf("%s editing, %s", con, dose)
    comparison <- "guide vs matched no-guide (paired difference vs 0)"
    if (n < 3 || dose == "scramble") {
      rows[[length(rows)+1]] <- data.frame(
        endpoint, comparison, n, paired = "yes",
        test = if (dose == "scramble") "none (negative control, descriptive only)" else "none (n<3, descriptive only)",
        mean_net = round(m,3), sd = round(s,3),
        statistic = NA, df = NA, p_raw = NA, p_adj = NA,
        cohens_dz = NA, ci95_low = NA, ci95_high = NA)
    } else {
      tt <- t.test(v, mu = 0)
      rows[[length(rows)+1]] <- data.frame(
        endpoint, comparison, n, paired = "yes",
        test = "one-sample t on paired differences (= paired t vs no-guide)",
        mean_net = round(m,3), sd = round(s,3),
        statistic = round(unname(tt$statistic),3),
        df = unname(tt$parameter),
        p_raw = signif(tt$p.value,4), p_adj = NA,   # filled below
        cohens_dz = round(m/s,3),
        ci95_low = round(tt$conf.int[1],3),
        ci95_high = round(tt$conf.int[2],3))
    }
  }
}
res <- do.call(rbind, rows)

# BH across the guide-vs-control tests that actually produced a P value
tested <- !is.na(res$p_raw)
res$p_adj[tested] <- signif(p.adjust(res$p_raw[tested], method = "BH"), 4)

res <- res[order(res$endpoint), ]
write.csv(res, "results/stats_table.csv", row.names = FALSE)
cat("wrote results/stats_table.csv\n")
print(res[, c("endpoint","n","test","statistic","df","p_raw","p_adj","cohens_dz","ci95_low","ci95_high")], row.names = FALSE)
