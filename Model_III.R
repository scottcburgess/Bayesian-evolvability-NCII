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
  filter(interaction %in% interactions & block %in% blocks)
hatch_settle.df <- hatch_settle.df |> 
  filter(interaction %in% interactions & block %in% blocks)

hs.df <- hatch_settle.df |> 
  mutate(
    block = as.numeric(factor(block)),
    sire = as.numeric(factor(sire)),
    dam = as.numeric(factor(dam)),
    interaction = as.numeric(factor(interaction))
  )
hatch.df <- filter(hs.df, metric == 'hatching')
settle.df <- filter(hs.df, metric == 'settling')


# Summarize hatching and settling rate across blocks
hatch_settle.df |> 
  group_by(interaction, metric) |> 
  summarize(pct = mean(outcome, na.rm = TRUE), .groups = 'drop') |> 
  pivot_wider(names_from = 'metric', values_from = 'pct') |> 
  as.data.frame()

hatch_settle.df |> 
  group_by(block, metric) |> 
  summarize(
    n = n(),
    n.settle = sum(outcome), 
    pr.settle = n.settle / n,
    .groups = 'drop'
  ) |> 
  mutate(
    a = n.settle + 1,
    b = n - n.settle + 1,
    lower = qbeta(0.0001, a, b),
    upper = qbeta(0.9999, a, b)
  ) |> 
  as.data.frame()


# Run model
post <- run.jags(
  data = list(
    n.blocks = length(unique(trunk_tail.df$block)),
    n.sires = length(unique(trunk_tail.df$sire)),
    n.dams = length(unique(trunk_tail.df$dam)),
    n.interactions = length(unique(trunk_tail.df$interaction)),
    n.larvae = nrow(trunk_tail.df),
    block = as.numeric(factor(trunk_tail.df$block)),
    sire = as.numeric(factor(trunk_tail.df$sire)),
    dam = as.numeric(factor(trunk_tail.df$dam)),
    interaction = as.numeric(factor(trunk_tail.df$interaction)),
    block.mean.range = cbind(
      round(range(trunk_tail.df$trunk)), 
      round(range(trunk_tail.df$tail)),
      qlogis(c(0.2, 0.95))
    ),
    length1 = cbind(trunk_tail.df$trunk, trunk_tail.df$tail),
    length2 = cbind(trunk_tail.df$trunk, trunk_tail.df$tail),
    length3 = cbind(trunk_tail.df$trunk, trunk_tail.df$tail),
    n.hatch = nrow(hatch.df),
    hatch.sire = hatch.df$sire,
    hatch = hatch.df$outcome,
    n.settle = nrow(settle.df),
    settle.block = settle.df$block,
    settle.sire = settle.df$sire,
    settle.dam = settle.df$dam,
    settle.interaction = settle.df$interaction,
    settle1 = settle.df$outcome,    
    settle2 = settle.df$outcome, 
    settle3 = settle.df$outcome,
    i = c(1, 1, 2),
    j = c(2, 3, 3)
  ),
  model = 'model {
    # for each t-trait...
    for(t in 1:3) {  
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
    for(k in 1:3) {
      sire.vcov[i[k], j[k]] <- sire.corr * sqrt(sire.vcov[i[k], i[k]] * sire.vcov[j[k], j[k]])
      sire.vcov[j[k], i[k]] <- sire.vcov[i[k], j[k]]
      dam.vcov[i[k], j[k]]  <- dam.corr * sqrt(dam.vcov[i[k], i[k]] * dam.vcov[j[k], j[k]])
      dam.vcov[j[k], i[k]] <- dam.vcov[i[k], j[k]]
      interaction.vcov[i[k], j[k]] <- interaction.corr * sqrt(interaction.vcov[i[k], i[k]] * interaction.vcov[j[k], j[k]])
      interaction.vcov[j[k], i[k]] <- interaction.vcov[i[k], j[k]]
      resid.vcov[i[k], j[k]] <- resid.corr * sqrt(resid.vcov[i[k], i[k]] * resid.vcov[j[k], j[k]])
      resid.vcov[j[k], i[k]] <- resid.vcov[i[k], j[k]]
    }
    
    # ---- prior for additive sire effect (for each sire) ----
    for(s in 1:n.sires) {
      sire.eff[s, 1:3] ~ dmnorm.vcov(sire.mean, sire.vcov)
      sire.eff.hatch[s] ~ dnorm(0, 1e-5)
      
      # effect on probability of settling for each sire
      logit(p.settle.sire[s]) <- sire.eff.hatch[s] + sire.eff[s, 3] 
    }
    
    # ---- prior for maternal effect (for each dam) ----
    for(d in 1:n.dams) {
      dam.eff[d, 1:3] ~ dmnorm.vcov(dam.mean, dam.vcov)
    }
    
    # ---- prior for interaction effect (for each sire x dam interaction) ----
    for(int in 1:n.interactions) {
      interaction.eff[int, 1:3] ~ dmnorm.vcov(interaction.mean, interaction.vcov)
    }
    
    # ---- trunk/tail likelihood ----
    for(l in 1:n.larvae) {
      for(t in 1:2) {        
        # likelihood of overall mean for computing evolvability
        length1[l, t] ~ dnorm(mean.overall[t], 1 / var.overall[t])
        length2[l, t] ~ dnorm(overall.block.mean[block[l], t], 1 / overall.block.var[block[l], t])
        
        # expected mean for the l-th larvae and t-th trait
        length.mu[l, t] <- block.mean[block[l], t] + 
          sire.eff[sire[l], t] + 
          dam.eff[dam[l], t] +
          interaction.eff[interaction[l], t]
      }
      # likelihood of l-th larvae for both traits from multivariate normal
      length3[l, ] ~ dmnorm.vcov(length.mu[l, ], resid.vcov[1:2, 1:2])
      
      # draw for posterior predictive check
      length.ppd[l, 1:2] ~ dmnorm.vcov(length.mu[l, ], resid.vcov[1:2, 1:2])
    }
    
    # ---- likelihood of hatching ----
    for(h in 1:n.hatch) {    
      # prior on non-sire effects 
      non.sire.eff.hatch[h] ~ dnorm(0, 1e-5)
      
      # linear model to compute probability of hatching 
      logit(pr.hatch[h]) <- non.sire.eff.hatch[h] +
        sire.eff.hatch[hatch.sire[h]] 
        
      hatch[h] ~ dbern(pr.hatch[h])
    }
    
    # ---- likelihood of settling given hatching ----
    for(s in 1:n.settle) {
      # settling block mean is on logit scale, so must take inverse-logit for bernoulli likelihood
      settle1[s] ~ dbern(ilogit(mean.overall[3]))
      settle2[s] ~ dbern(ilogit(overall.block.mean[settle.block[s], 3]))
      
      # linear model to compute probability of settling given hatching
      logit(pr.settle.hatch[s]) <- block.mean[settle.block[s], 3] +
        sire.eff[settle.sire[s], 3] +
        dam.eff[settle.dam[s], 3] +
        interaction.eff[settle.interaction[s], 3]
        
      # probability of settling given hatching likelihood
      settle3[s] ~ dbern(pr.settle.hatch[s])
      
      # draw for posterior predictive check
      settle.ppd[s] ~ dbern(pr.settle.hatch[s])
    }
  }',
  monitor = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov', 
    'resid.vcov', 'overall.block.mean', 'mean.overall', 
    'length.ppd', 'settle.ppd', 'sire.eff', 'p.settle.sire'
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
  dimnames(p$sire.eff)[[2]] <- c('Trunk', 'Tail', 'Settling')
dimnames(p$length.ppd)[[2]] <- c('Trunk', 'Tail')
dimnames(p$sire.vcov)[1:2] <-
  dimnames(p$dam.vcov)[1:2] <-
  dimnames(p$interaction.vcov)[1:2] <-
  dimnames(p$resid.vcov)[1:2] <- 
  list(dimnames(p$overall.block.mean)[[2]], dimnames(p$overall.block.mean)[[2]])


# Hardcode settling resid.vcov off-diag to 0 and diag to 1
p$resid.vcov['Settling', , ] <- p$resid.vcov[, 'Settling', ] <- 0
p$resid.vcov['Settling', 'Settling', ] <- 1

# Add QG metrics to list
p <- addQGmetrics(p)

vcv.obs <- convertVCVscale.III(p)

beta <- sapply(1:dim(p$p.settle.sire)[2], function(i) {
  cov.i <- c(
    cov.trunk.settle = cov(p$sire.eff[, 'Trunk', i], p$p.settle.sire[, i]),
    cov.tail.settle = cov(p$sire.eff[, 'Tail', i], p$p.settle.sire[, i])
  )
  inv.vcov <- solve(p$sire.vcov[c('Trunk', 'Tail'), c('Trunk', 'Tail'), i])
  beta <- cov.i %*% inv.vcov
  c(cov.i, beta = beta[1, ])
}) |> 
  t()


# CODA summary ------------------------------------------------------------

post.smry <- smrzPost(post, c(
  'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov', 'resid.vcov'
))


# Posterior Predictive Check ----------------------------------------------

length.obs <- cbind(Trunk = trunk_tail.df$trunk, Tail = trunk_tail.df$tail)

ppc <- expand_grid(
  metric = colnames(length.obs),
  id = 1:nrow(length.obs)
) |> 
  mutate(metric = factor(metric, colnames(length.obs))) |> 
  bind_rows(
    data.frame(metric = 'Settling', id = 1:nrow(settle.df))
  )

ppc$pct.gte.obs <- sapply(1:nrow(ppc), function(i) {
  id <- ppc$id[i]
  
  obs <- if(ppc$metric[i] == 'Settling') {
    settle.df$outcome[id]
  } else {
    length.obs[id, ppc$metric[id]]
  }
  
  ppd <- if(ppc$metric[i] == 'Settling') {
    p$settle.ppd[id, ]
  } else {
    p$length.ppd[id, ppc$metric[id], ]
  }
  
  mean(obs >= ppd)
})

ppc$mean.diff <- sapply(1:nrow(ppc), function(i) {
  id <- ppc$id[i]
  
  obs <- if(ppc$metric[i] == 'Settling') {
    settle.df$outcome[id]
  } else {
    length.obs[id, ppc$metric[id]]
  }
  
  ppd <- if(ppc$metric[i] == 'Settling') {
    p$settle.ppd[id, ]
  } else {
    p$length.ppd[id, ppc$metric[id], ]
  }
  
  mean(obs - ppd)
})

ppc.smry <- smrzPPC(ppc)



# Save all objects
save.image(format(end, 'Model_outputs/Model_III_posterior_%Y%m%d_%H%M.rdata'))


# Plot posterior distributions
plot(
  post,
  vars = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov', 'resid.vcov'
  ),
  file = format(end, 'Model_outputs/Model_III_plots_%Y%m%d_%H%M.pdf')
)


# Plot diagnostics
pdf(format(end, "Model_outputs/Model_III_diagnostics_%Y%m%d_%H%M.pdf"))

ggplot(post.smry) +
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