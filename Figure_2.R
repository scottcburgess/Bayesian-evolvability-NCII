rm(list = ls())
library('tidyverse')
library('ggridges')

# Trunk Tail: load and prepare ----
load("Model_outputs/posterior_20260622_0449.rdata")

source('0_misc_funcs.R')

d_Trunktail <- lapply(vcov_params, function(param) {
  data.frame(
    Trunk = p[[param]]['Trunk', 'Trunk', ] / p$VP['Trunk', 'Trunk', ],
    Tail = p[[param]]['Tail', 'Tail', ] / p$VP['Tail', 'Tail', ],
    param = param
  )
}) |> 
  bind_rows() |> 
  mutate(param = factor(param, levels = rev(unique(param))))

panelA <- d_Trunktail |> 
  rename(value = 'Trunk') |> 
  plot_func(
    title = 'a) Trunk length', 
    bw = 0.01, 
    breaks = seq(0, 1, 0.1), 
    max_x = 1,
    param_df = param_df
  )

panelB <- d_Trunktail |> 
  rename(value = 'Tail') |> 
  plot_func(
    title = 'b) Tail length', 
    bw = 0.01, 
    breaks = seq(0, 1, 0.1), 
    max_x = 1,
    param_df = param_df
  )


# Create plot -------------------------------------------------------------
fig2 <- gridExtra::grid.arrange(
  panelA, panelB,
  nrow = 2,
  bottom = 'Proportion of phenotypic variance',
  left = 'Variance component'
)
fig2


# Save plot ----
ggsave(
  "Figures and Tables/Figure 2.pdf", 
  plot = fig2, 
  height = 4, 
  width = 3
)


# Summaries ----
smry <- bind_rows(
  d_Trunktail |> 
    group_by(param) |> 
    summarize(Trunk = list(vecSmry(Trunk))) |> 
    unnest_wider(Trunk) |> 
    mutate(type = 'Trunk'),
  d_Trunktail |> 
    group_by(param) |>
    summarize(Tail = list(vecSmry(Tail))) |> 
    unnest_wider(Tail) |> 
    mutate(type = 'Tail')
) |> 
  mutate(across(-c(type, param), function(x) round(x, 4))) |> 
  select(type, param, everything())
