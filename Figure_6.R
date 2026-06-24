rm(list=ls())
library('tidyverse')
library('ggridges')

load("Model_outputs/posterior_20260622_0449.rdata") 
source('0_misc_funcs.R')

panelA <- sg |> 
  filter(!is.na(sg.pct), W == "Hatch") |> 
  rename(value = 'sg.pct', param = z) |> 
  mutate(param = factor(param, levels = unique(sg.df$z))) |> 
  plot_func(
    title = 'a) Hatching',
    bw = 0.008, 
    breaks = seq(-1, 1, 0.1), 
    min_x = -0.7,
    max_x = 0.7,
    param_df = data.frame(
      param = sg.df$z,
      color = c('#a6dba0', '#1b7837'),
      param_label = sapply(sg.df$z, function(x) paste0('`', x, '`'))
    )
  ) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  labs(x = '', y = '') +
  theme(
    axis.title.x = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

panelB <- sg |> 
  filter(!is.na(sg.pct), W == "Settle|Hatch") |> 
  rename(value = 'sg.pct', param = z) |> 
  mutate(param = factor(param, levels = unique(sg.df$z))) |> 
  plot_func(
    title = 'b) Settlement',
    bw = 0.008, 
    breaks = seq(-1, 1, 0.1), 
    min_x = -0.7,
    max_x = 0.7,
    param_df = data.frame(
      param = sg.df$z,
      color = c('#c2a5cf', '#762a83'),
      param_label = sapply(sg.df$z, function(x) paste0('`', x, '`'))
    )
  ) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  labs(x = '', y = 'Density') +
  theme(
    axis.title.x = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

panelC <- sg |> 
  rename(value = 'sg.pct', param = z) |> 
  mutate(param = factor(param, levels = rev(sg.df$z[1:2]))) |> 
  pivot_wider(
    names_from = W,
    values_from = value,
    id_cols = c(block, sample, param)
  ) |> 
  mutate(value = Hatch + `Settle|Hatch`) |> 
  select(value, param) |> 
  plot_func(
    title = 'c) Total response',
    bw = 0.008, 
    breaks = seq(-1, 1, 0.1), 
    min_x = -0.9,
    max_x = 0.9,
    param_df = data.frame(
      param = sg.df$z,
      color = c('#fdb863', '#e66101'),
      param_label = sapply(sg.df$z, function(x) paste0('`', x, '`'))
    )
  ) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  labs(x = '', y = '') +
  theme(
    axis.title.x = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



fig6 <- gridExtra::grid.arrange(
  panelA, panelB, panelC, nrow = 3,
  bottom = grid::textGrob('Genetic response to selection (%)', vjust = -2)
)
fig6

ggsave(
  "Figures and Tables/Figure 6.pdf",
  plot = fig6,
  height = 7,
  width = 5
)
