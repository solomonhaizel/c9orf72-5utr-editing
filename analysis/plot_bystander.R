# Amplicon-wide adenosine scan (bystander), read from data/bystander_by_rep.csv.
# Scores the 16 reference adenosines in the biological ~30-96 window; target
# sites A74/A88 are drawn black. Six conditions: A74 and A88 guide at each dose,
# No guide, Scramble. Raw editing per replicate (points), mean bars.
suppressMessages(library(ggplot2)); set.seed(1)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
d <- read.csv("data/bystander_by_rep.csv", stringsAsFactors = FALSE)
lv <- c("A74 0.3 uM","A74 1 uM","A88 0.3 uM","A88 1 uM","No guide","Scramble")
d$condition <- factor(d$condition, levels = lv)
d$tgt <- factor(ifelse(d$site_type == "target", "Target (A74/A88)", "Other adenosine"),
                levels = c("Other adenosine","Target (A74/A88)"))
agg <- aggregate(edit_pct ~ condition + position + tgt, d,
                 function(x) c(m = mean(x), se = sd(x)/sqrt(length(x))))
agg <- data.frame(condition = agg$condition, position = agg$position, tgt = agg$tgt,
                  mean = agg$edit_pct[,"m"], se = agg$edit_pct[,"se"])
p <- ggplot() +
  geom_vline(xintercept = c(74,88), linetype = "dotted", colour = "grey45", linewidth = 0.4) +
  geom_col(data = agg, aes(position, mean, fill = tgt), width = 1.6) +
  geom_errorbar(data = agg, aes(position, ymin = mean - se, ymax = mean + se),
                width = 1.4, linewidth = 0.3) +
  geom_point(data = d, aes(position, edit_pct), colour = "grey25", size = 0.5, alpha = 0.55,
             position = position_jitter(width = 0.9, height = 0, seed = 1)) +
  facet_wrap(~condition, ncol = 2) +
  scale_fill_manual(values = c("Other adenosine" = "grey70", "Target (A74/A88)" = "black")) +
  scale_x_continuous(breaks = seq(40,100,20), limits = c(26,99)) +
  labs(x = "Biological 5'UTR position (1-98)", y = "A-to-I editing (%)", fill = NULL,
       subtitle = "amplicon-wide adenosine scan; mean bars, points = replicates") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        strip.background = element_rect(fill = "grey85", colour = "black"),
        text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
        plot.subtitle = element_text(size = 10))
ggsave("results/figures/bystander_scan.pdf", p, width = 8.5, height = 6)
ggsave("results/figures/bystander_scan.png", p, width = 8.5, height = 6, dpi = 600)
cat("wrote bystander scan\n")
