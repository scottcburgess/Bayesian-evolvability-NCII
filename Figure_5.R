rm(list=ls())
options(scipen = 999)
library('tidyverse')
library('ggridges')

# trunk Tail: load and prepare ----
load("Model_outputs/posterior_20260605_0120.rdata") 
source('0_misc_funcs.R')

## Panel A ----
panelA <- e.params_beta |>
  filter(Beta_index == "1") |>
  select(e, c, r) |>
  setNames(betas) |> 
  pivot_longer(
    cols = c(eB, cB, rB),
    names_to = "param",
    values_to = "value"
  ) |>
  mutate(
    param = factor(param, levels = betas),
    value = value * 100
  ) |> 
  plot_func(
    title = 'a) Selection for longer tails', 
    bw = 0.005, 
    breaks = seq(0, 1, 0.05), 
    max_x = 0.25
  ) +
  theme(plot.title = element_text(size = 9))


## Panel B ----
panelB <- e.params_beta |>
  filter(Beta_index == "2") |>
  select(e, c, r) |>
  setNames(betas) |> 
  pivot_longer(
    cols = c(eB, cB, rB),
    names_to = "param",
    values_to = "value"
  ) |>
  mutate(
    param = factor(param, levels = betas),
    value = value * 100
  ) |> 
  plot_func(
    title = 'b) Selection for short trunks, short tails', 
    bw = 0.04,
    breaks = seq(0, 1, 0.05), 
    max_x = 0.25
  ) +
  theme(plot.title = element_text(size = 9))


## Panel C ----
panelC <- e.params_beta |>
  filter(Beta_index == "3") |>
  select(e, c, r) |>
  setNames(betas) |> 
  pivot_longer(
    cols = c(eB, cB, rB),
    names_to = "param",
    values_to = "value"
  ) |>
  mutate(
    param = factor(param, levels = betas),
    value = value * 100
  ) |> 
  plot_func(
    title = 'c) Selection for long trunks, short tails', 
    bw = 0.04,
    breaks = seq(0, 1, 0.05), 
    max_x = 0.25
  ) +
  theme(plot.title = element_text(size = 9))


fig5 <- gridExtra::grid.arrange(
  panelA, panelB, panelC, nrow = 1,
  bottom = 'Evolvability (%)',
  left = 'Metric'
)
fig5

ggsave(
  "Figures and Tables/Figure 5.pdf", 
  plot = fig5, 
  height = 3, 
  width = 8
)
