rm(list = ls())
library(tidyverse)
library(runjags)


# Model I -----------------------------------------------------------------

load('3_Model_outputs/Model_I_posterior_20250918_0036.rdata')


# CODA --------------------------------------------------------------------

post.smry <- summary(post) 

post.smry2 <- post.smry |>  
  as.data.frame() |> 
  rownames_to_column('metric') |>
  select(metric, SSeff:psrf) |> 
  filter(!stringr::str_detect(metric, 'deviance|ppc')) |> 
  # filter(stringr::str_detect(metric, 'vcov')) |> 
  pivot_longer(-metric, names_to = 'diag', values_to = 'values') 

post.smry2 |> 
  group_by(diag) |> 
  summarize(
    median = median(values),
    lower = unname(quantile(values, 0.025)),
    upper = unname(quantile(values, 0.975)),
    .groups = 'drop'
  )

ggplot(post.smry2) +
  geom_histogram(aes(values), bins = 40) +
  facet_wrap(~diag, scales = 'free_x')



# Posterior Predictive Check ----------------------------------------------

length.obs <- cbind(Trunk = df$head, Tail = df$tail)
dimnames(p$length.ppc)[[2]] <- colnames(length.obs)

ppc <- expand_grid(
  metric = colnames(length.obs),
  id = 1:nrow(length.obs)
) |> 
  mutate(metric = factor(metric, colnames(length.obs)))
  
ppc$pct.gte.obs <- sapply(1:nrow(ppc), function(i) {
  obs <- length.obs[ppc$id[i], ppc$metric[i]]
  ppc <- p$length.ppc[ppc$id[i], ppc$metric[i], ]
  mean(obs >= ppc)
})

ppc$mean.diff <- sapply(1:nrow(ppc), function(i) {
  obs <- length.obs[ppc$id[i], ppc$metric[i]]
  ppc <- p$length.ppc[ppc$id[i], ppc$metric[i], ]
  mean(obs - ppc)
})

ppc.smry <- ppc |> 
  group_by(metric) |> 
  summarize(
    median.pct = median(pct.gte.obs),
    lower.pct = unname(quantile(pct.gte.obs, 0.025)),
    upper.pct = unname(quantile(pct.gte.obs, 0.975)),   
    median.diff = median(mean.diff),
    lower.diff = unname(quantile(mean.diff, 0.025)),
    upper.diff = unname(quantile(mean.diff, 0.975)),
    .groups = 'drop'
  )
ppc.smry

ggplot(ppc) +
  geom_histogram(aes(pct.gte.obs), binwidth = 0.05) +
  geom_vline(aes(xintercept = median.pct), data = ppc.smry, color = 'red') +
  geom_vline(aes(xintercept = lower.pct), data = ppc.smry, linetype = 'dashed', color = 'red') +
  geom_vline(aes(xintercept = upper.pct), data = ppc.smry, linetype = 'dashed', color = 'red') +
  facet_wrap(~ metric) +
  labs(x = 'Percent of PPD >= Observed', y = 'Count')

ggplot(ppc) +
  geom_histogram(aes(mean.diff), bins = 50) +
  geom_vline(aes(xintercept = median.diff), data = ppc.smry, color = 'red') +
  geom_vline(aes(xintercept = lower.diff), data = ppc.smry, linetype = 'dashed', color = 'red') +
  geom_vline(aes(xintercept = upper.diff), data = ppc.smry, linetype = 'dashed', color = 'red') +
  facet_wrap(~ metric, scales = 'free_x') +
  labs(x = 'Length Difference (Observed - PPD)', y = 'Count')
