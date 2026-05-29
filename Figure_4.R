rm(list=ls())
options(scipen = 999)
library('tidyverse')
library('ggridges')
source('0_misc_funcs.R')


# trunk Tail: load and prepare ----
load("Model_outputs/Model_I_posterior_20250930_0019.rdata") 
# Get the posterior samples by averaging across all of the random beta's
E_trunk_tail <- lapply(betas, function(b) {
  data.frame(
    beta = apply(e.params_BetaMCMC$post.dist[[b]], 1, mean) * 100,
    param = b
  )
}) |> 
  bind_rows() |> 
  mutate(param = factor(param, levels = unique(param))) 

# Check
# e.params_BetaMCMC$summary # median (called e_mean) should be the same as
# vecSmry(eB_posterior_trunk_tail) # median here (but we're using the mode)

panelA <- E_trunk_tail |> 
  rename(value = 'beta') |> 
  plot_func(
    title = 'a) Trunk-Tail length', 
    bw = 0.0015, 
    breaks = seq(0, 0.1, 0.01), 
    max_x = 0.08
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Hatching Settling : load and prepare ----
load("Model_outputs/Model_III_posterior_20250930_1554.rdata") 
E_hatch_settle <- lapply(betas, function(b) {
  data.frame(
    beta = apply(e.params_BetaMCMC$post.dist[[b]], 1, mean) * 100,
    param = b
  )
}) |> 
  bind_rows() |> 
  mutate(param = factor(param, levels = unique(param)))

panelB <- E_hatch_settle |> 
  rename(value = 'beta') |> 
  plot_func(
    title = 'b) Hatch-Settle probability', 
    bw = 0.1, 
    breaks = seq(0, 10, 1), 
    max_x = 6
  )


# Save plot ----
fig4 <- gridExtra::grid.arrange(
  panelA, panelB, ncol = 2,
  bottom = 'Evolvability (%)',
  left = 'Metric'
)
fig4


ggsave(
  "Figures and Tables/Figure 4.pdf", 
  plot = fig4, 
  height = 3, 
  width = 6
)


# Summaries ----
smry <- bind_rows(
  E_hatch_settle  |> 
    group_by(param) |> 
    summarize(beta = list(vecSmry(beta))) |> 
    unnest_wider(beta) |> 
    mutate(type = 'hatch_settle'),
  E_trunk_tail |> 
    group_by(param) |>
    summarize(beta = list(vecSmry(beta))) |> 
    unnest_wider(beta) |> 
    mutate(type = 'trunk_tail')
) |> 
  mutate(across(-c(type, param), function(x) round(x, 4))) |> 
  select(type, param, everything())