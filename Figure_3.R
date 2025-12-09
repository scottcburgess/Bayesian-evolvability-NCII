rm(list = ls())
library('tidyverse')
library('ggridges')
source('0_misc_funcs.R')


## Panel A ----

# Trunk Tail: load and prepare ----
load("Model_outputs/Model_I_posterior_20250930_0019.rdata")

# quick check 
# vecSmry(p$H[1,1,]) # is ~similar to
# vecSmry(p$VA[1,1,] / p$VP[1,1,])
# vecSmry(vcv.obs$VA$vcv.G.obs[1,1,] / vcv.obs$VA$vcv.P.obs[1,1,])

d_Trunktail <- lapply(vcov_params, function(param) {
  data.frame(
    Trunk = p[[param]][1, 1, ] / p$VP[1, 1, ],
    Tail = p[[param]][2, 2, ] / p$VP[2, 2, ],
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


## Panel C ----

# Trunk Tail ratio: load and prepare ----
load("Model_outputs/Model_II_posterior_20251001_1520.rdata")

## Prepare data
d_Trunktail_ratio <- lapply(vcov_params, function(param) {
  data.frame(
    ratio = var.obs[[param]]$var.a.obs / var.obs[[param]]$var.obs,
    param = param
  )
}) |>  
  bind_rows() |> 
  mutate(param = factor(param, levels = rev(unique(param))))

panelC <- d_Trunktail_ratio |> 
  rename(value = 'ratio') |> 
  plot_func(
    title = 'c) Trunk:Tail ratio', 
    bw = 0.01, 
    breaks = seq(0, 1, 0.1), 
    max_x = 0.8
  )


## Panel D ----

# Hatching Settling : load and prepare ----
load("Model_outputs/Model_III_posterior_20250930_1554.rdata") 

# quick check for hatching
# vecSmry(p$H[1,]) # is ~similar to
# vecSmry(vcv.obs$VA$vcv.G.obs[1,1,] / vcv.obs$VA$vcv.P.obs[1,1,])

## Prepare data 
d_hatchsettle <- lapply(c('VA', 'VM', 'VD'), function(param) {
  data.frame(
    hatch = vcv.obs[[param]]$vcv.G.obs[1, 1, ] / vcv.obs[[param]]$vcv.P.obs[1, 1, ],
    settle = vcv.obs[[param]]$vcv.G.obs[2, 2, ] / vcv.obs[[param]]$vcv.P.obs[2, 2, ],
    param = param
  )
}) |> 
  bind_rows() |> 
  mutate(param = factor(param, levels = rev(unique(param))))

panelD <- d_hatchsettle |> 
  rename(value = 'hatch') |> 
  plot_func(
    title = 'd) Hatching probability', 
    bw = 0.01, 
    breaks = seq(0, 1, 0.1), 
    max_x = 0.8
  )


## Panel E ----

panelE <- d_hatchsettle |> 
  rename(value = 'settle') |> 
  plot_func(
    title = 'e) Settlement probability', 
    bw = 0.01, 
    breaks = seq(0, 1, 0.1), 
    max_x = 0.8
  )



# Create plot -------------------------------------------------------------
fig3 <- gridExtra::grid.arrange(
  panelA, panelB, panelC, panelD, panelE,
  nrow = 2, ncol = 3,
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
    mutate(type = 'Tail'),
  d_Trunktail_ratio |> 
    group_by(param) |>
    summarize(ratio = list(vecSmry(ratio))) |> 
    unnest_wider(ratio) |> 
    mutate(type = 'ratio'),  
  d_hatchsettle |> 
    group_by(param) |>
    summarize(hatch = list(vecSmry(hatch))) |> 
    unnest_wider(hatch) |> 
    mutate(type = 'Hatch'),
  d_hatchsettle |> 
    group_by(param) |>
    summarize(settle = list(vecSmry(settle))) |> 
    unnest_wider(settle) |> 
    mutate(type = 'Settle')
) |> 
  mutate(across(-c(type, param), function(x) round(x, 4))) |> 
  select(type, param, everything())