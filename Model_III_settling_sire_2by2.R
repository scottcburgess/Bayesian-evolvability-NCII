rm(list = ls())
library(tidyverse)
library(runjags)
source('0_misc_funcs.R')

# MCMC parameters
chains <- 6 #50
adapt <- 100
burnin <- 1000 #100000
total.sample <- 1000 #8000 
thin <- 1 #5000

# Load data
trunk_tail.df <- readRDS('Data/trunk_tail_data.rds') 
hatch_settle.df <- readRDS('Data/hatch_settle_data.rds') 

# Select blocks and interactions to shared between trunk/tail and hatching/settling data
blocks <- hatch_settle.df |> 
  group_by(block) |> 
  summarize(
    n.hatch = sum(metric == 'hatching'),
    n.settle = sum(metric == 'settling'),
    .groups = 'drop'
  ) |> 
  filter(n.hatch > 0 & n.settle > 0) |> 
  pull(block) |> 
  unique()

interactions <- intersect(trunk_tail.df$interaction, hatch_settle.df$interaction)

trunk_tail.df <- trunk_tail.df |> 
  # filter trunk/tail data for blocks and interactions to use
  filter(interaction %in% interactions & block %in% blocks) |> 
  mutate(
    block = as.numeric(factor(block)),
    sire = as.numeric(factor(sire)),
    dam = as.numeric(factor(dam)),
    interaction = as.numeric(factor(interaction))
  )


hatch_settle.df <- hatch_settle.df |> 
  # filter hatching/settling data for blocks and interactions to use 
  filter(interaction %in% interactions & block %in% blocks) |> 
  mutate(
    block = as.numeric(factor(block)),
    sire = as.numeric(factor(sire))
  ) |> 
  # compress to number of successes and trials for hatching and settling by interaction
  group_by(block, sire) |> 
  summarize(
    n.hatched.trial = sum(metric == 'hatching'),
    n.hatched = sum(metric == 'hatching' & outcome == 1),
    pct.hatched = n.hatched / n.hatched.trial,
    n.settled.hatched.trial = sum(metric == 'settling'),
    n.settled.hatched = sum(metric == 'settling' & outcome == 1),
    pct.settled.hatched= n.settled.hatched / n.settled.hatched.trial,
    pr.settle = pct.hatched * pct.settled.hatched,
    .groups = 'drop'
  )


# Run model
post <- run.jags(
  data = list(
    n.blocks = n_distinct(trunk_tail.df$block),
    n.sires = n_distinct(trunk_tail.df$sire),
    n.dams = n_distinct(trunk_tail.df$dam),
    n.interactions = n_distinct(trunk_tail.df$interaction),
    n.larvae = nrow(trunk_tail.df),
    tt.block = trunk_tail.df$block,
    tt.sire = trunk_tail.df$sire,
    tt.dam = trunk_tail.df$dam,
    tt.interaction = trunk_tail.df$interaction,
    block.mean.range = cbind(
      round(range(trunk_tail.df$trunk)), 
      round(range(trunk_tail.df$tail)),
      qlogis(c(0.2, 0.95))
    ),
    length1 = cbind(trunk_tail.df$trunk, trunk_tail.df$tail),
    length2 = cbind(trunk_tail.df$trunk, trunk_tail.df$tail),
    length3 = cbind(trunk_tail.df$trunk, trunk_tail.df$tail),
    hs.block = hatch_settle.df$block,
    k.hatch = hatch_settle.df$n.hatched.trial,
    n.hatch = hatch_settle.df$n.hatched,
    k.settle.hatch = hatch_settle.df$n.settled.hatched.trial,
    n.settle.hatch = hatch_settle.df$n.settled.hatched
  ),
  model = 'model {
    # for each t-trait...
    for(t in 1:2) {  
      # ---- prior for overall mean and variance ----
      mean.overall[t] ~ dunif(block.mean.range[1, t], block.mean.range[2, t])
      var.overall[t] ~ dunif(0, 1e5)
      
      # ---- priors for block and effect means ----
      for(b in 1:n.blocks) {
        overall.block.mean[b, t] ~ dunif(block.mean.range[1, t], block.mean.range[2, t])
        overall.block.var[b, t] ~ dunif(0, 1e5)
        block.mean[b, t] ~ dunif(block.mean.range[1, t], block.mean.range[2, t])
      }
      sire.mean[t] ~ dnorm(0, 1e-5)
      dam.mean[t] ~ dnorm(0, 1e-5)
      interaction.mean[t] ~ dnorm(0, 1e-5)
      
      # ---- prior for variances ----
      sire.vcov[t, t] ~ dunif(0, 1e3)
      dam.vcov[t, t] ~ dunif(0, 1e3)
      interaction.vcov[t, t] ~ dunif(0, 1e3)
      resid.vcov[t, t] ~ dunif(0, 1e3)
    }
    
    # ---- correlation priors for covariances ----
    sire.corr ~ dunif(-1, 1)
    dam.corr ~ dunif(-1, 1)
    interaction.corr ~ dunif(-1, 1)
    resid.corr ~ dunif(-1, 1)
    
    # ---- construct variance/covariance matrices ----
    sire.vcov[1, 2] <- sire.corr * sqrt(sire.vcov[1, 1] * sire.vcov[2, 2])
    sire.vcov[2, 1] <- sire.vcov[1, 2]
    dam.vcov[1, 2]  <- dam.corr * sqrt(dam.vcov[1, 1] * dam.vcov[2, 2])
    dam.vcov[2, 1] <- dam.vcov[1, 2]
    interaction.vcov[1, 2] <- interaction.corr * sqrt(interaction.vcov[1, 1] * interaction.vcov[2, 2])
    interaction.vcov[2, 1] <- interaction.vcov[1, 2]
    resid.vcov[1, 2] <- resid.corr * sqrt(resid.vcov[1, 1] * resid.vcov[2, 2])
    resid.vcov[2, 1] <- resid.vcov[1, 2]
    
    # ---- priors on block mean of hatching and settling given hatching
    for(b in 1:n.blocks) {
      hatch.block.mean[b] ~ dunif(block.mean.range[1, 3], block.mean.range[2, 3])
      settle.hatch.block.mean[b] ~ dunif(block.mean.range[1, 3], block.mean.range[2, 3])
    }
    
    for(s in 1:n.sires) {
      # ---- prior for additive sire effect for trunk/tail (for each sire) ----
      sire.eff[s, 1:2] ~ dmnorm.vcov(sire.mean, sire.vcov)
      
      # ---- probability of hatching ----
      # prior on hatching sire effect
      hatch.sire.eff[s] ~ dnorm(0, 1e-3)
      
      # likelihood of hatching
      pr.hatch[s] <- ilogit(hatch.block.mean[hs.block[s]] + hatch.sire.eff[s])
      n.hatch[s] ~ dbinom(pr.hatch[s], k.hatch[s])
      
      # draw for posterior predictive check
      hatch.ppd[s] ~ dbinom(pr.hatch[s], k.hatch[s])
    
      # ---- probability of settling given hatching ----
      # prior on settling given hatching sire effect
      settle.hatch.sire.eff[s] ~ dnorm(0, 1e-3)
      
      # likelihood of settling given hatching
      pr.settle.hatch[s] <- ilogit(settle.hatch.block.mean[hs.block[s]] + settle.hatch.sire.eff[s])
      n.settle.hatch[s] ~ dbinom(pr.settle.hatch[s], k.settle.hatch[s])
      
      # draw for posterior predictive check
      settle.hatch.ppd[s] ~ dbinom(pr.settle.hatch[s], k.settle.hatch[s])
      
      # ---- overall probability of settling ----
      pr.settle[s] <- pr.hatch[s] * pr.settle.hatch[s]
      
      # ---- effect of sire on overall probability of settling ----
      pr.settle.wo.sire[s] <- ilogit(hatch.block.mean[hs.block[s]]) * ilogit(settle.hatch.block.mean[hs.block[s]])
      pr.settle.sire.eff[s] <- pr.settle[s] - pr.settle.wo.sire[s]
    }
    
    # ---- prior for maternal effect (for each dam) ----
    for(d in 1:n.dams) {
      dam.eff[d, 1:2] ~ dmnorm.vcov(dam.mean, dam.vcov)
    }
    
    # ---- prior for interaction effect (for each sire x dam interaction) ----
    for(int in 1:n.interactions) {
      interaction.eff[int, 1:2] ~ dmnorm.vcov(interaction.mean, interaction.vcov)
    }
    
    # ---- trunk/tail likelihood ----
    for(l in 1:n.larvae) {
      for(t in 1:2) {        
        # likelihood of overall mean for computing evolvability
        length1[l, t] ~ dnorm(mean.overall[t], 1 / var.overall[t])
        length2[l, t] ~ dnorm(overall.block.mean[tt.block[l], t], 1 / overall.block.var[tt.block[l], t])
        
        # expected mean for the l-th larvae and t-th trait
        length.mu[l, t] <- block.mean[tt.block[l], t] + 
          sire.eff[tt.sire[l], t] + 
          dam.eff[tt.dam[l], t] +
          interaction.eff[tt.interaction[l], t]
      }
      # likelihood of l-th larvae for both traits from multivariate normal
      length3[l, ] ~ dmnorm.vcov(length.mu[l, ], resid.vcov[1:2, 1:2])
      
      # draw for posterior predictive check
      length.ppd[l, 1:2] ~ dmnorm.vcov(length.mu[l, ], resid.vcov[1:2, 1:2])
    }
  }',
  monitor = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov', 
    'resid.vcov', 'overall.block.mean', 'mean.overall', 
    'length.ppd',  'sire.eff',
    'pr.hatch', 'pr.settle.hatch', 'pr.settle', 'pr.settle.sire.eff',
    'hatch.ppd', 'settle.hatch.ppd'
  ), 
  inits = function() list(
    .RNG.name = 'lecuyer::RngStream',
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c('glm', 'lecuyer'),
  n.chains = chains,
  adapt = adapt,
  burnin = burnin,
  sample = ceiling(total.sample / chains),
  thin = thin,
  method = 'parallel',
  summarise = FALSE
)
end <- Sys.time()
elapsed <- swfscMisc::autoUnits(post$timetaken)



# Extract posterior to list of arrays - p
p <- swfscMisc::runjags2list(post)
dimnames(p$overall.block.mean)[[2]] <- 
  dimnames(p$sire.eff)[[2]] <-
  dimnames(p$length.ppd)[[2]] <- c('Trunk', 'Tail')
dimnames(p$sire.vcov)[1:2] <-
  dimnames(p$dam.vcov)[1:2] <-
  dimnames(p$interaction.vcov)[1:2] <-
  dimnames(p$resid.vcov)[1:2] <- 
  list(dimnames(p$overall.block.mean)[[2]], dimnames(p$overall.block.mean)[[2]])

# Add QG metrics to list
p <- addQGmetrics(p)

beta <- sapply(1:dim(p$pr.settle.sire)[2], function(i) {
  cov.i <- c(
    cov.trunk.settle = cov(p$sire.eff[, 'Trunk', i], p$pr.settle.sire.eff[, i]),
    cov.tail.settle = cov(p$sire.eff[, 'Tail', i], p$pr.settle.sire.eff[, i])
  )
  inv.vcov <- solve(p$sire.vcov[c('Trunk', 'Tail'), c('Trunk', 'Tail'), i])
  beta <- cov.i %*% inv.vcov
  c(cov.i, beta = beta[1, ])
}) |> 
  t()


# CODA summary ------------------------------------------------------------

post.smry <- smrzPost(post, c(
  'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov', 'resid.vcov',
  'pr.hatch', 'pr.settle.hatch', 'pr.settle', 'pr.settle.sire.eff' 
))


# Posterior Predictive Check ----------------------------------------------

length.obs <- cbind(Trunk = trunk_tail.df$trunk, Tail = trunk_tail.df$tail)

ppc <- expand_grid(
  metric = colnames(length.obs),
  id = 1:nrow(length.obs)
) |> 
  mutate(metric = factor(metric, colnames(length.obs))) |> 
  bind_rows(
    data.frame(metric = 'Settling', id = 1:nrow(hatch_settle.df))
  )

ppc <- ppc |> 
  cbind(sapply(1:nrow(ppc), function(i) {
    id <- ppc$id[i]
    
    obs <- if(ppc$metric[i] == 'Settling') {
      hatch_settle.df$n.settled.hatched[id]
    } else {
      length.obs[id, ppc$metric[id]]
    }
    
    ppd <- if(ppc$metric[i] == 'Settling') {
      p$settle.hatch.ppd[id, ]
    } else {
      p$length.ppd[id, ppc$metric[id], ]
    }
    
    c(
      pct.gte.obs = mean(obs >= ppd), 
      mean.diff = mean(obs - ppd)
    )
  }) |> 
    t()
  )

ppc.smry <- smrzPPC(ppc)



# Save all objects
save.image(format(end, 'Model_outputs/Model_III_posterior_%Y%m%d_%H%M.rdata'))


# Plot posterior distributions
plot(
  post,
  vars = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov', 'resid.vcov', 
    'pr.hatch', 'pr.settle.hatch', 'pr.settle', 'pr.settle.sire.eff' 
  ),
  file = format(end, 'Model_outputs/Model_III_plots_%Y%m%d_%H%M.pdf')
)


# Plot diagnostics
pdf(format(end, "Model_outputs/Model_III_diagnostics_%Y%m%d_%H%M.pdf"))

ggplot(post.smry$post) +
  geom_histogram(aes(values), bins = 20) +
  facet_wrap(~diag, scales = 'free_x')

ggplot(ppc) +
  geom_histogram(aes(pct.gte.obs), binwidth = 0.05) +
  geom_vline(aes(xintercept = median.pct), data = ppc.smry, color = 'red') +
  geom_vline(aes(xintercept = lower.pct), data = ppc.smry, linetype = 'dashed', color = 'red') +
  geom_vline(aes(xintercept = upper.pct), data = ppc.smry, linetype = 'dashed', color = 'red') +
  facet_wrap(~ metric, scales = 'free') +
  labs(x = 'Percent of PPD >= Observed', y = 'Count')

ggplot(ppc) +
  geom_histogram(aes(mean.diff), bins = 50) +
  geom_vline(aes(xintercept = median.diff), data = ppc.smry, color = 'red') +
  geom_vline(aes(xintercept = lower.diff), data = ppc.smry, linetype = 'dashed', color = 'red') +
  geom_vline(aes(xintercept = upper.diff), data = ppc.smry, linetype = 'dashed', color = 'red') +
  facet_wrap(~ metric, scales = 'free') +
  labs(x = 'Metric Difference (Observed - PPD)', y = 'Count')

dev.off()


print(elapsed)