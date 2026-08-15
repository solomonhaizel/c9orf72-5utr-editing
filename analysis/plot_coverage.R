# Target-site coverage per condition, read from data/editing_by_rep.csv.
# Replicate identity is explicit (rep column), not implied by row order.
suppressMessages(library(ggplot2)); set.seed(1)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
d <- read.csv("data/editing_by_rep.csv", stringsAsFactors = FALSE)
d$guide_cov <- suppressWarnings(as.numeric(d$guide_cov))
d <- d[!is.na(d$guide_cov), ]
d$grp <- paste(d$construct, d$dose)
agg <- aggregate(guide_cov ~ grp, d, function(x) c(m = mean(x), se = sd(x)/sqrt(length(x))))
agg <- data.frame(grp = agg$grp, mean = agg$guide_cov[,"m"], se = agg$guide_cov[,"se"])
p <- ggplot(agg, aes(mean, grp)) +
  geom_col(fill = "grey80", colour = "black", width = 0.7) +
  geom_point(data = d, aes(guide_cov, grp), shape = 21, fill = "white", size = 2,
             position = position_jitter(height = 0.12, width = 0, seed = 1)) +
  geom_vline(xintercept = 10000, linetype = "dashed") +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Target-site coverage (reads, guide well)", y = NULL,
       subtitle = "per replicate; dashed line = 10,000-read QC threshold") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), text = element_text(colour = "black"),
        axis.text = element_text(colour = "black"), plot.subtitle = element_text(size = 10))
ggsave("results/figures/coverage.pdf", p, width = 6, height = 3.6)
ggsave("results/figures/coverage.png", p, width = 6, height = 3.6, dpi = 600)
cat("wrote coverage figure\n")
