rm(list = ls())
library(tidyverse)
source('0_misc_funcs.R')

load("Model_outputs/posterior_20260605_0120.rdata")

param.df <- param_df |> 
  mutate(
    pr.gt0 = sapply(
      param,
      function(x) round(mean(p[[x]]['Trunk', 'Tail', ] > 0), digits = 2)
    ),
    label = paste0(
      'atop(', 
      title, '~(', short, '),Pr(', short, '[12] > 0) == ', pr.gt0, 
      ')'
    ),
    label = factor(label, levels = label)
  )


# draw random values from multivariate-normal for a posterior sample
ran_mvnorm <- function(vec) {
  Sigma <- matrix(vec, nrow = 2, byrow = TRUE)
  if(!all(eigen(Sigma)$values > 0)) return(NULL)
  
  MASS::mvrnorm(100, c(0, 0), Sigma) |> 
    as.data.frame() 
}


# compute ellipses at selected intervals from random multivariate normal draws
smrz_matrix <- function(param, p) {
  # random multivariate normal draws from mean-centered covariance matrices
  mat <- evolvability::meanStdGMCMC(
    t(apply(p[[param]][names.2, names.2, ], 3, as.vector)),
    t(p$overall.mean)
  ) 
  
  rads <- apply(mat, 1, function(vec) {
    eigenvectors <- vec |> 
      matrix(nrow = 2, byrow = TRUE) |> 
      eigen() |> 
      pluck('vectors')
    atan2(eigenvectors[2, 1], eigenvectors[1, 1])
  })
  
  pts <- apply(mat, 1, ran_mvnorm, simplify = FALSE) |> 
    bind_rows() 
  mu <- colMeans(pts)
  sigma <- cov(pts)
  
  list(
    slope.smry = c(
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
      mutate(param = param),
    
    # ellipses at selected intervals
    ellipses = lapply(
      seq(0.05, 0.95, length.out = 10),
      function(pr) {
        ellipse::ellipse(sigma, centre = mu, level = pr) |> 
          as.data.frame() |> 
          setNames(c('x', 'y')) |> 
          mutate(pr = pr)
      }
    ) |> 
      bind_rows() |> 
      mutate(param = param)
  )
}

# get ellipses and slope summaries for each parameter
matrix.smry <- lapply(c(vcov_params, 'VP'), smrz_matrix, p = p) 

# extract ellipses
ellipses <- matrix.smry |> 
  lapply(function(x) x$ellipses) |> 
  bind_rows() |> 
  left_join(select(param.df, param, label), by = 'param') |> 
  arrange(desc(pr))

# extract slope summaries
slope.smry <- lapply(matrix.smry, function(x) x$slope.smry) |> 
  bind_rows() |> 
  mutate(intercept = 0) |> 
  left_join(select(param.df, param, label), by = 'param')

# axis limits to make figure square
lims <- unlist(ellipses[, c('x', 'y')]) |> 
  pretty() |> 
  range()

fig2 <- ellipses |> 
  ggplot() +
  geom_hline(yintercept = 0, linewidth = 1, color = "gray", alpha = 0.6) +
  geom_vline(xintercept = 0, linewidth = 1, color = "gray", alpha = 0.6) +
  geom_polygon(
    aes(x, y, color = label, fill = label, alpha = 1 - pr, group = pr),
    linewidth = 0.1
  ) +
  geom_abline(
    aes(slope = lower, intercept = intercept, color = label), 
    data = slope.smry,
    linetype = 'dashed',
    linewidth = 0.4
  ) +
  geom_abline(
    aes(slope = upper, intercept = intercept, color = label), 
    data = slope.smry,
    linetype = 'dashed',
    linewidth = 0.4
  ) +  
  geom_abline(
    aes(slope = median, intercept = intercept, color = label), 
    data = slope.smry,
    linewidth = 0.6
  ) +
  scale_color_manual(
    values = param.df |> 
      select(label, color) |> 
      deframe()
  ) +
  scale_fill_manual(
    values = param.df |> 
      select(label, color) |> 
      deframe()
  ) +
  coord_equal() +
  labs(
    x = "Trunk length (mean standardized)",
    y = "Tail length (mean standardized)"
  ) +
  lims(x = lims, y = lims) +
  facet_wrap(~label, labeller = label_parsed) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = 'none'
  )
fig2

ggsave(
  "Figures and Tables/Figure 2.pdf",
  plot = fig2,
  height = 5,
  width = 5
)
