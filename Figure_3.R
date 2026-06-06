rm(list = ls())
library('tidyverse')
library('ggridges')

## Panel A ----

# Trunk Tail: load and prepare ----
load("Model_outputs/posterior_20260605_0120.rdata")

source('0_misc_funcs.R')

# quick check 
# vecSmry(p$H[1,1,]) # is ~similar to
# vecSmry(p$VA[1,1,] / p$VP[1,1,])
# vecSmry(vcv.obs$VA$vcv.G.obs[1,1,] / vcv.obs$VA$vcv.P.obs[1,1,])

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
    max_x = 0.8
  )

panelB <- d_Trunktail |> 
  rename(value = 'Tail') |> 
  plot_func(
    title = 'b) Tail length', 
    bw = 0.01, 
    breaks = seq(0, 1, 0.1), 
    max_x = 0.8
  )


# Create plot -------------------------------------------------------------
fig3 <- gridExtra::grid.arrange(
  panelA, panelB,
  nrow = 2,
  bottom = 'Proportion of phenotypic variance',
  left = 'Variance component'
)
fig3


# Save plot ----
ggsave(
  "Figures and Tables/Figure 3.pdf", 
  plot = fig3, 
  height = 4, 
  width = 8
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
