param_df <- data.frame(
  param = c('VA', 'VM', 'VD', 'VR', 'VP', 'eB', 'cB', 'rB'),
  param_label = c('V[A]', 'V[M]', 'V[D]', 'V[R]', 'V[P]', 'e(beta)', 'c(beta)', 'r(beta)'),
  color = c(
    '#E76F51', '#E9C46A', '#2a9d8f', 'grey', '#118ab2', '#0077b6', '#00b4d8', '#73e8ff'
  ),
  short = c('G', 'M', 'D', 'R', 'P', 'e', 'c', 'r'),
  title = c(
    'Additive~genetic', 'Maternal~effect', 'Dominance', 'Residual', 'Phenotypic',
    'e', 'c', 'r'
  )
)

vcov_params <- c('VA', 'VM', 'VD', 'VR')

betas <- c('eB', 'rB', 'cB')

vc_colors <- param_df |> 
  select(param, color) |> 
  deframe()

default_theme <- ggridges::theme_ridges() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 10, face = "plain")
  )

plot_func <- function(df, title, bw, breaks, max_x) {
  smry <- df |>
    group_by(param) |>
    summarise(summary_values = list(vecSmry(value)), .groups = 'drop') |>
    unnest_wider(summary_values, names_repair = "unique")
  
  df |> 
    ggplot(aes(x = value, y = param, fill = param)) +
    geom_density_ridges(
      scale = 0.9,
      alpha = 0.4,
      bandwidth = bw,
      color = 'lightgrey',
      linewidth = 0.1
    ) + 
    geom_segment(
      aes(
        x = lower.hdi, 
        xend = upper.hdi, 
        y = as.numeric(param), 
        yend = as.numeric(param),
        color = param
      ), 
      data = smry, 
      linetype = "solid", 
      linewidth = 0.5
    ) +
    geom_point(
      aes(x = mode, y = as.numeric(param), color = param), 
      data = smry, 
      size = 2
    ) +
    labs(title = title) +
    scale_fill_manual(values = vc_colors) +
    scale_color_manual(values = vc_colors) +
    scale_x_continuous(breaks = breaks, limits = c(0, max_x)) +
    scale_y_discrete(
      labels = param_df |> 
        select(param, param_label) |> 
        deframe() |> 
        sapply(function(x) parse(text = x))
    ) + 
    default_theme
}