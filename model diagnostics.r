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
dimnames(p$length.ppc)[[2]] <- rownames(length.obs)

ppc <- expand_grid(
  metric = rownames(length.obs),
  id = 1:nrow(length.obs)
) |> 
  mutate(metric = factor(metric, rownames(length.obs)))
  
ppc$pct.gte.obs <- sapply(1:nrow(ppc), function(i) {
  obs <- length.obs[ppc$id[i], ppc$metric[i]]
  ppc <- p$length.ppc[ppc$id[i], ppc$metric[i], ]
  mean(obs >= ppc)
})

ppc.smry <- ppc |> 
  group_by(metric) |> 
  summarize(
    median = median(pct.gte.obs),
    lower = unname(quantile(pct.gte.obs, 0.025)),
    upper = unname(quantile(pct.gte.obs, 0.975)),
    .groups = 'drop'
  )
ppc.smry

ggplot(ppc) +
  geom_histogram(aes(pct.gte.obs), binwidth = 0.05) +
  geom_vline(aes(xintercept = median), data = ppc.smry, color = 'red') +
  geom_vline(aes(xintercept = lower), data = ppc.smry, linetype = 'dashed', color = 'red') +
  geom_vline(aes(xintercept = upper), data = ppc.smry, linetype = 'dashed', color = 'red') +
  facet_wrap(~ metric) +
  labs(x = 'Percent of PPD >= Observed', y = 'Count')

