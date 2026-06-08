rm(list=ls())
library('tidyverse')
library('ggridges')

load("Model_outputs/posterior_20260604_2220.rdata") 
source('0_misc_funcs.R')

# sg |> 
#   filter(!is.na(sg)) |> 
#   mutate(
#     z.W = paste0(z, ':', W),
#     block = factor(block)
#   ) |> 
#   ggplot() +
#   ggridges::geom_density_ridges(aes(x = sg, y = block)) +
#   xlim(c(-2, 3)) +
#   facet_wrap(~z.W, ncol = 2, scales = 'free_y')

fig6 <- sg |> 
  filter(!is.na(sg)) |> 
  mutate(z.W = paste0(z, ':', W)) |> 
  ggplot() +
  geom_density(
    aes(x = sg),
    alpha = 0.4,
    fill = 'gray30'
  ) +
  # Plot 'sg_percent', not sg
  # Add median and 95% hdpi
  # Add x = 'Evolvability (%), y = 'Density'
  xlim(c(-4, 4)) +
  facet_wrap(~z.W, ncol = 2, scales = 'free_y') +
  default_theme
fig6

ggsave(
  "Figures and Tables/Figure 6.pdf",
  plot = fig6,
  height = 4,
  width = 6
)
