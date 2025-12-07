rm(list = ls())
library(tidyverse)

load("Model_outputs/Model_I_posterior_20250930_0019.rdata")
p$VR <- p$resid.vcov

param.df <- data.frame(
  param = c('VA', 'VM', 'VD', 'VR', 'VP'),
  color = c('#ef476f', '#118ab2', '#ffd166', '#06d6a0', 'grey30'),
  short = c('G', 'M', 'D', 'R', 'P'),
  title = c(
    'Additive~genetic', 'Maternal~effect', 'Dominance', 'Residual', 'Phenotypic'
  )
) |> 
  mutate(
    pr.gt0 = sapply(
      param,
      function(x) round(mean(p[[x]][1, 2, ] > 0), digits = 2)
    ),
    label = paste0(
      'atop(', 
      title, '~(', short, '),Pr(', short, '[12] > 0) == ', pr.gt0, 
      ')'
    ),
    label = factor(label, levels = label)
  )


slope.smry <- lapply(
  param.df$param, 
  function(param) {
    mat <- evolvability::meanStdGMCMC(
      t(apply(p[[param]], 3, as.vector)),
      t(p$mean.overall)
    )
    
    rads <- apply(mat, 1, function(vec) {
      eigenvectors <- vec |> 
        matrix(nrow = 2, byrow = TRUE) |> 
        eigen() |> 
        pluck('vectors')
      eigenvectors <- ifelse(eigenvectors < 0, eigenvectors * -1, eigenvectors)
      atan2(eigenvectors[2, 1], eigenvectors[1, 1])
    })
    
    c(
      median = rads |> 
        circular::circular(units = 'radians') |> 
        circular::median.circular() |> 
        tan(),
      rads |> 
        HDInterval::hdi() |> 
        tan()
    ) |> 
      rbind() |> 
      as.data.frame() |> 
      mutate(param = param)
  }
) |> 
  bind_rows() |> 
  mutate(intercept = 0) |> 
  left_join(select(param.df, param, label), by = 'param')


lims <- c(-0.15, 0.15)

p1 <- slope.smry |> 
  ggplot() +
  geom_hline(yintercept = 0, linewidth = 1, color = "gray") +
  geom_vline(xintercept = 0, linewidth = 1, color = "gray") +
  geom_abline(
    aes(slope = median, intercept = intercept, color = label), 
    linewidth = 1
  ) + 
  geom_abline(
    aes(slope = lower, intercept = intercept, color = label), 
    alpha = 0.75,
    linetype = 'dashed',
    linewidth = 1
  ) +
  geom_abline(
    aes(slope = upper, intercept = intercept, color = label), 
    alpha = 0.75,
    linetype = 'dashed',
    linewidth = 1
  ) +
  coord_equal() +
  labs(x = "Trunk length", y = "Tail length") +
  lims(x = lims, y = lims) +
  scale_color_manual(
    values = param.df |> 
      select(label, color) |> 
      deframe()
  ) +
  facet_wrap(~label, labeller = label_parsed) +
  theme_minimal(base_size = 18) +
  theme(legend.position = 'none')
p1


ggsave(
  "Figures and Tables/major_eigenvector_posterior.pdf",
  plot = p1,
  height = 8,
  width = 8
)
