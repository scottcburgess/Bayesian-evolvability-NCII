rm(list = ls())
library(tidyverse)

load("Model_outputs/Model_I_posterior_20250930_0019.rdata")
p$VR <- p$resid.vcov

param.df <- data.frame(
  param = c('VA', 'VM', 'VD', 'VR', 'VP'),
  color = c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00'),
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


# draw random values from multivariate-normal for a posterior sample
ran_mvnorm <- function(vec) {
  Sigma <- matrix(vec, nrow = 2, byrow = TRUE)
  if(!all(eigen(Sigma)$values > 0)) return(NULL)
  
  MASS::mvrnorm(100, c(0, 0), Sigma) |> 
    as.data.frame() 
}


# compute ellipses at selected intervals from random multivariate normal draws
ran_ellipses <- function(param, p) {
  # random multivariate normal draws from mean-centered covariance matrices
  pts <- evolvability::meanStdGMCMC(
    t(apply(p[[param]], 3, as.vector)),
    t(p$mean.overall)
  ) |> 
    apply(1, ran_mvnorm, simplify = FALSE) |> 
    bind_rows() 
  
  mu <- colMeans(pts)
  sigma <- cov(pts)
  eigenvectors <- eigen(sigma)$vectors
  
  list(
    # major and minor axes of data
    axes = data.frame(
      major = eigenvectors[2, 1] / eigenvectors[1, 1],
      minor = eigenvectors[2, 2] / eigenvectors[1, 2],
      param = param
    ),
    # ellipses at selected intervals
    ellipses =  lapply(
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

# get ellipses and axes for each sample
ellipse_samples <- lapply(param.df$param, ran_ellipses, p = p) 

# extract ellipses
ellipses <- ellipse_samples |> 
  lapply(function(x) x$ellipses) |> 
  bind_rows() |> 
  left_join(select(param.df, param, label), by = 'param') 

# extract axes
axes <- ellipse_samples |> 
  lapply(function(x) x$axes) |> 
  bind_rows() |> 
  left_join(select(param.df, param, label), by = 'param') |> 
  mutate(intercept = 0)

# axis limits to make figure square
lims <- unlist(ellipses[, c('x', 'y')]) |> 
  pretty() |> 
  range()

p1 <- ellipses |> 
  arrange(pr) |> 
  ggplot() +
  geom_hline(yintercept = 0, linewidth = 1, color = "gray") +
  geom_vline(xintercept = 0, linewidth = 1, color = "gray") +
  geom_abline(
    aes(slope = major, intercept = intercept, color = label), 
    data = axes,
    linewidth = 1
  ) +
  geom_abline(
    aes(slope = minor, intercept = intercept, color = label), 
    data = axes,
    linewidth = 1
  ) +
  geom_polygon(
    aes(x, y, color = label, fill = label, alpha = 1 - pr, group = pr),
    linewidth = 0.1
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
  labs(x = "Trunk length", y = "Tail length") +
  lims(x = lims, y = lims) +
  facet_wrap(~label, labeller = label_parsed) +
  theme_minimal(base_size = 18) +
  theme(legend.position = 'none')
p1

ggsave(
  "Figures and Tables/Visualize_matrices_ea.pdf",
  plot = p1,
  height = 8,
  width = 8
)
