rm(list = ls())
library(tidyverse)
library(runjags)

# MCMC parameters
chains <- 10
adapt <- 100
burnin <- 50000
total.sample <- 50000 
thin <- 100

# Load data
trunk_tail.df <- readRDS('../1_Data/head_tail_data.rds') 
hatch_settle.df <- readRDS('../1_Data/hatch_settle_data.rds') 

# Filter for blocks that occur in both
blocks.to.keep <- table(
  block = hatch_settle.df$block,
  metric = hatch_settle.df$metric
) |> 
  as.data.frame() |> 
  filter(Freq > 0 & metric == 'settling') |>
  pull(block) |> 
  as.character() |> 
  as.integer()

hatch_settle.df <- filter(hatch_settle.df, block %in% blocks.to.keep) 
interactions <- intersect(trunk_tail.df$interaction, hatch_settle.df$interaction)
trunk_tail.df <- filter(trunk_tail.df, interaction %in% interactions)
hatch_settle.df <- filter(hatch_settle.df, interaction %in% interactions)

# Summarize settling rate across blocks
block.settle <- hatch_settle.df |> 
  group_by(block) |> 
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
  )


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
      round(range(trunk_tail.df$head)), 
      round(range(trunk_tail.df$tail)),
      qlogis(c(0.4, 0.95))
    ),
    length1 = cbind(trunk_tail.df$head, trunk_tail.df$tail),
    length2 = cbind(trunk_tail.df$head, trunk_tail.df$tail),
    length3 = cbind(trunk_tail.df$head, trunk_tail.df$tail),
    n.settle = nrow(hatch_settle.df),
    settle.block = as.numeric(factor(hatch_settle.df$block)),
    settle.sire = as.numeric(factor(hatch_settle.df$sire)),
    settle.dam = as.numeric(factor(hatch_settle.df$dam)),
    settle.interaction = as.numeric(factor(hatch_settle.df$interaction)),
    settle1 = hatch_settle.df$outcome,    
    settle2 = hatch_settle.df$outcome, 
    settle3 = hatch_settle.df$outcome,
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
    }
    
    # ---- prior for maternal effect (for each dam) ----
    for(d in 1:n.dams) {
      dam.eff[d, 1:3] ~ dmnorm.vcov(dam.mean, dam.vcov)
    }
    
    # ---- prior for interaction effect (for each sire x dam interaction) ----
    for(int in 1:n.interactions) {
      interaction.eff[int, 1:3] ~ dmnorm.vcov(interaction.mean, interaction.vcov)
    }
    
    
    # ---- head/tail likelihood ----
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
      length.ppc[l, 1:2] ~ dmnorm.vcov(length.mu[l, ], resid.vcov[1:2, 1:2])
    }
    
    
    # ---- likelihood of settling ----
    for(s in 1:n.settle) {
      # settling block mean is on logit scale, so must take inverse-logit for bernoulli likelihood
      settle1[s] ~ dbern(ilogit(mean.overall[3]))
      settle2[s] ~ dbern(ilogit(overall.block.mean[settle.block[s], 3]))
      
      logit(pr.settle[s]) <- block.mean[settle.block[s], 3] +
        sire.eff[settle.sire[s], 3] +
        dam.eff[settle.dam[s], 3] +
        interaction.eff[settle.interaction[s], 3]
      settle3[s] ~ dbern(pr.settle[s])
      
      # draw for psoterior predictive check
      settle.ppc[s] ~ dbern(pr.settle[s])
    }
  }',
  monitor = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov', 
    'resid.vcov', 'block.mean', 'sire.eff', 'dam.eff', 
    'interaction.eff', 'mean.overall', 'overall.block.mean',
    'length.ppc', 'settle.ppc'
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
dimnames(p$block.mean)[[2]] <- 
  dimnames(p$sire.eff)[[2]] <- 
  dimnames(p$dam.eff)[[2]] <- 
  dimnames(p$interaction.eff)[[2]] <- 
  dimnames(p$mean.overall)[[1]] <- 
  dimnames(p$overall.block.mean)[[2]] <- c('Trunk', 'Tail', 'Settling')
dimnames(p$sire.vcov)[1:2] <-
  dimnames(p$dam.vcov)[1:2] <-
  dimnames(p$interaction.vcov)[1:2] <-
  dimnames(p$resid.vcov)[1:2] <- 
  list(dimnames(p$sire.eff)[[2]], dimnames(p$sire.eff)[[2]])

# Hardcode settling resid.vcov off-diag to 0 and diag to 1
p$resid.vcov['Settling', , ] <- p$resid.vcov[, 'Settling', ] <- 0
p$resid.vcov['Settling', 'Settling', ] <- 1

# Add QG metrics to list
p$VA <- 4 * p$sire.vcov
p$VM <- p$dam.vcov - p$sire.vcov
p$VD <- 4 * p$interaction.vcov
p$VP <- p$VA + p$VM + p$VD + p$resid.vcov
p$H <- p$VA / p$VP


# Compute heritability and evolvability based on deVillemereuil et al 2016
qgparams.post <- sapply(dimnames(p$overall.block.mean)[[2]], function(m) {
  parallel::mclapply(1:dim(p$VA)[3], function(i) {
    QGglmm::QGparams(
      var.a = p$VA[m, m, i],
      var.p = p$VP[m, m, i],
      predict = p$overall.block.mean[, m, i],
      model = if(m == 'Settling') 'binom1.logit' else 'Gaussian',
      verbose = FALSE
    )
  }, mc.cores = 14) |>
    bind_rows() |>
    mutate(E = var.a.obs / (p$mean.overall[m, ] ^ 2))
}, simplify = FALSE)


# Save all objects and plot posterior summaries
save.image(format(end, '../3_Model_outputs/Model_I_posterior_%Y%m%d_%H%M.rdata'))

plot(
  post,
  file = format(end, '../3_Model_outputs/Model_I_plots_%Y%m%d_%H%M.pdf')
)

print(elapsed)
