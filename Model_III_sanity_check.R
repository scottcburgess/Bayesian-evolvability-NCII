rm(list = ls())
library(tidyverse)
source('0_misc_funcs.R')
load('Model_outputs/Model_III_posterior_20260417_2041.rdata')


pr.settle.post <- p$pr.settle |>
  as.data.frame() 
pr.settle.post$sire <- 1:22
pr.settle.post <- pr.settle.post |> 
  pivot_longer(-sire, names_to = 'sample', values_to = 'pr.settle') |> 
  select(sire, pr.settle)
ggplot(pr.settle.post, aes(pr.settle)) +
  geom_histogram() +
  geom_vline(aes(xintercept = pr.settle), data = hatch_settle.df) +
  facet_wrap(~ sire)


# Calculate the mean sire effect for settle probability
sire_effect_pr_settle <- hatch_settle.df |>
  mutate(
    block_pr_settle = median(pr.settle),
    .by = block
  ) |> 
  summarize(
    sire_effect_pr_settle = pr.settle - block_pr_settle,
    .by = sire
  )


# Calculate the mean trunk and tail per sire and block
sire_effect_trunk_tail <- trunk_tail.df |>
  mutate(
    block_trunk = mean(trunk),
    block_tail = mean(tail),
    .by = block
  ) |> 
  summarize(
    Trunk = mean(trunk) - mean(block_trunk),
    Tail = mean(tail) - mean(block_tail),
    .by = sire
  )

# Combine data sets for the 'raw' data
d <- left_join(
  sire_effect_pr_settle, 
  sire_effect_trunk_tail, 
  by = 'sire'
)

# Get the averages of the posteriors for the sire effects
# p$sire.eff is [sires,trait,draw]
d_post <- data.frame(
  apply(p$sire.eff, c(1, 2) , median),
  sire_effect_pr_settle = apply(p$pr.settle.sire.eff, 1, median)
)


# Sanity check 1 ----
## Visualize the raw means per sire vs the posterior means per sire ----
### Trunk ----
ggplot(mapping = aes(Trunk, sire_effect_pr_settle)) +
  geom_point(data = d, color = 'red') +
  geom_point(data = d_post, color = 'blue') +
  theme_bw()

### Tail ----
ggplot(mapping = aes(Tail, sire_effect_pr_settle)) +
  geom_point(data = d, color = 'red') +
  geom_point(data = d_post, color = 'blue') +
  theme_bw()


# Sanity check 2 ----
## Visualize posteriors
beta_long <- beta %>% 
  as.data.frame() %>%
  mutate(draw = row_number()) %>%
  pivot_longer(
    cols = -draw,
    names_to = "parameter",
    values_to = "value")

ggplot(data = beta_long,
       aes(x = value)) +
  geom_histogram(bins = 300, alpha = 0.4) +
  scale_x_continuous(limits = c(-1, 1)) +
  facet_wrap( ~ parameter, scales = "free") +
  theme_bw()


# Sanity check 3 ----
## Covariance from raw data ----
d |> summarize(trunk_cov = cov(sire_effect_trunk, sire_effect_pr_settle),
               tail_cov = cov(sire_effect_tail, sire_effect_pr_settle)) 

## Covariance from Model III ----
apply(beta, 2, vecSmry)
