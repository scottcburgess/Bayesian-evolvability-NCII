rm(list = ls())
library(tidyverse)
source('0_misc_funcs.R')
load('Model_outputs/Model_III_posterior_20260417_2041.rdata')


pr.settle.post <- as.data.frame(p$pr.settle) 
colnames(pr.settle.post) <- make.unique(colnames(pr.settle.post))
pr.settle.post |> 
  mutate(sire = 1:22) |> 
  pivot_longer(-sire, names_to = 'sample', values_to = 'pr.settle') |> 
  select(sire, pr.settle) |> 
  ggplot(aes(pr.settle)) +
  geom_histogram(alpha = 0.5) +
  geom_vline(aes(xintercept = pr.settle), data = hatch_settle.df) +
  facet_wrap(~ sire) +
  theme_minimal()


pr.settle.sire.eff <- as.data.frame(p$pr.settle.sire.eff)
colnames(pr.settle.sire.eff) <- make.unique(colnames(pr.settle.sire.eff))
pr.settle.sire.eff |> 
  mutate(sire = 1:22) |> 
  pivot_longer(-sire, names_to = 'sample', values_to = 'pr.settle') |> 
  select(sire, pr.settle) |> 
  ggplot(aes(pr.settle)) +
  geom_histogram(alpha = 0.5) +
  geom_vline(aes(xintercept = 0)) +
  facet_wrap(~ sire) +
  theme_minimal()


# plot pr.settle.sire.eff vs trunk/tail sire.eff (are they random?)
sire.eff <- p$sire.eff
x <- dimnames(sire.eff)
names(x) <- c('sire', 'metric', 'sample')
x$sire <- as.character(1:22)
x$sample <- make.unique(x$sample)
dimnames(sire.eff) <- x
sire.eff <- sire.eff |> 
  as.data.frame.table(responseName = 'length') |> 
  mutate(
    sire = as.character(sire),
    sample = as.character(sample),
    metric = as.character(metric)
  ) |> 
  pivot_wider(
    id_cols = c(sire, sample),
    names_from = metric,
    values_from = length
  )

df <- pr.settle.sire.eff |> 
  mutate(sire = as.character(1:22)) |> 
  pivot_longer(-sire, names_to = 'sample', values_to = 'pr.settle.sire.eff') |> 
  left_join(sire.eff, by = c('sire', 'sample')) 
    
p1 <- ggplot(df, aes(Trunk, pr.settle.sire.eff)) +
  geom_point(color = 'grey', alpha = 0.5, size = 0.5) +
  geom_density2d(color = 'red') +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

p2 <- ggplot(df, aes(Tail, pr.settle.sire.eff)) +
  geom_point(color = 'grey', alpha = 0.5, size = 0.5) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_density2d(color = 'red') +
  theme_bw()

gridExtra::grid.arrange(p1, p2)



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
beta %>% 
  as.data.frame() %>%
  mutate(draw = row_number()) %>%
  pivot_longer(
    cols = -draw,
    names_to = "parameter",
    values_to = "value") |> 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 300, alpha = 0.4) +
  geom_vline(xintercept = 0) +
  # scale_x_continuous(limits = c(-1, 1)) +
  facet_wrap(~ parameter, scales = "free") +
  theme_bw()


# Sanity check 3 ----
## Covariance from raw data ----
# d |> summarize(trunk_cov = cov(sire_effect_trunk, sire_effect_pr_settle),
#                tail_cov = cov(sire_effect_tail, sire_effect_pr_settle)) 
# 
# ## Covariance from Model III ----
# apply(beta, 2, vecSmry)
