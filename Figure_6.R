rm(list=ls())
library('tidyverse')
library('ggridges')

load("Model_outputs/posterior_20260604_2220.rdata") 
source('0_misc_funcs.R')

sg |> 
  filter(!is.na(sg)) |> 
  mutate(
    x.W = paste0(x, ':', W),
    block = factor(block)
  ) |> 
  ggplot() +
  ggridges::geom_density_ridges(aes(x = sg, y = block)) +
  xlim(c(-2, 3)) +
  facet_wrap(~x.W, ncol = 2, scales = 'free_y')

sg |> 
  filter(!is.na(sg)) |> 
  mutate(x.W = paste0(x, ':', W)) |> 
  ggplot() +
  geom_density(aes(x = sg)) +
  xlim(c(-4, 4)) +
  facet_wrap(~x.W, ncol = 2, scales = 'free_y')

# ggsave(
#   "Figures and Tables/Figure 6.pdf", 
#   plot = fig6, 
#   height = 3, 
#   width = 6
# )