rm(list=ls())
options(scipen = 999)
library('tidyverse')
library('ggridges')

# trunk Tail: load and prepare ----
load("Model_outputs/posterior_20260605_0520.rdata") 

source('0_misc_funcs.R')

# Get the posterior samples by averaging across all of the random beta's
E_trunk_tail <- lapply(betas, function(b) {
  data.frame(
    beta = apply(e.params_BetaMCMC$post.dist[[b]], 1, mean) * 100,
    param = b
  )
}) |> 
  bind_rows() |> 
  mutate(param = factor(param, levels = unique(param))) 

fig4 <- E_trunk_tail |> 
  rename(value = 'beta') |> 
  plot_func(
    bw = 0.0015, 
    breaks = seq(0, 0.1, 0.01), 
    max_x = 0.1,,
    param_df = param_df
  ) +
  labs(x = 'Evolvability (%)', y = 'Metric') +
  theme(
    axis.title.x = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
fig4


ggsave(
  "Figures and Tables/Figure 4.pdf", 
  plot = fig4, 
  height = 3, 
  width = 3
)


# Summaries ----
smry <- E_trunk_tail |> 
  group_by(param) |>
  summarize(beta = list(vecSmry(beta))) |> 
  unnest_wider(beta) |> 
  mutate(across(-param, function(x) round(x, 4))) |> 
  select(param, everything())
