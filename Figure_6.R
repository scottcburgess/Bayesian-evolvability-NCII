rm(list=ls())
library('tidyverse')
library('ggridges')
source('0_misc_funcs.R')


# load and prepare ----
load("Model_outputs/Model_IV_posterior_20251001_0640.rdata") 

# Trunk selection differential  
cov_trunk_p <- vcv.obs$vcv.G.obs[1, 3, ]
mean_p <- vcv.obs$mean.obs[, "Settling"]
cov_trunk_p_relative <- cov_trunk_p / mean_p # selection differential

# Tail selection differential
cov_tail_p <- vcv.obs$vcv.G.obs[2, 3, ]
mean_p <- vcv.obs$mean.obs[, "Settling"]
cov_tail_p_relative <- cov_tail_p / mean_p # selection differential

# Extract means
mean_z <- rbind(p$mean.overall[1, ], p$mean.overall[2, ])

# Compute trunk and tail evolvability
sg <- rbind(cov_trunk_p_relative, cov_tail_p_relative)  # 2 x n_iter

n_iter <- dim(sg)[2]

E_matrix <- matrix(
  NA, nrow = 2, ncol = n_iter,
  dimnames = list(c("Trunk", "Tail"), NULL)
)

for (i in seq_len(n_iter)) {
  S_vec <- sg[, i]
  mean_vec <- mean_z[, i]
  E_matrix[, i] <- S_vec / mean_vec
}


# Make Figure ----
plot_fig6 <- function(df, title) {
  smry <- df |> 
    pull('value') |> 
    vecSmry() |> 
    t() |> 
    as.data.frame() |> 
    mutate(y = 0)
  
  df |> 
    ggplot(aes(x = value)) +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_density(
      alpha = 0.4,
      fill = 'gray30',
      linewidth = 0.1
    ) + 
    geom_segment(
      aes(x = lower.hdi, xend = upper.hdi),
      data = smry, 
      y = 0,
      linetype = "solid", 
      linewidth = 0.5
    ) +
    geom_point(
      aes(x = mode), 
      data = smry,
      y = 0, 
      size = 2
    ) +
    labs(title = title) +
    scale_x_continuous(breaks = seq(-0.5, 1, 0.1), limits = c(-0.5, 1)) +
    default_theme +
    theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1))
}


## Panel A ----
panelA <- plot_fig6(
  data.frame(value = E_matrix['Trunk', ] * 100),
  'a) Trunk length'
)


## Panel B ----
panelB <- plot_fig6(
  data.frame(value = E_matrix['Tail', ] * 100),
  'b) Tail length'
)


fig6 <- gridExtra::grid.arrange(
  panelA, panelB, nrow = 1,
  bottom = '% change per generation',
  left = 'Probability density'
)
fig6


ggsave(
  "Figures and Tables/Figure 6.pdf", 
  plot = fig6, 
  height = 3, 
  width = 6
)