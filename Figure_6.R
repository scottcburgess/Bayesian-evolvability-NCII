rm(list=ls())
library('tidyverse')
library('ggridges')

load("Model_outputs/posterior_20260617_0854.rdata") 
source('0_misc_funcs.R')

fig6 <- sg |> 
  filter(!is.na(sg.pct)) |> 
  rename(value = 'sg.pct', param = z.W) |> 
  mutate(param = factor(param, levels = rev(sg.df$z.W))) |> 
  select(value, param) |> 
  plot_func(
    bw = 0.015, 
    breaks = seq(-1, 1, 0.05), 
    min_x = -0.7,
    max_x = 0.7,
    param_df = data.frame(
      param = sg.df$z.W,
      color = c('#a6dba0', '#1b7837', '#c2a5cf', '#762a83'),
      param_label = sapply(sg.df$z.W, function(x) paste0('`', x, '`'))
    )
  ) +
  labs(x = 'Genetic response to selection (%)', y = 'Density') +
  theme(
    axis.title.x = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
fig6

ggsave(
  "Figures and Tables/Figure 6.pdf",
  plot = fig6,
  height = 4,
  width = 6
)
