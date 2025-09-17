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
df <- readRDS("1_Data/head_tail_data.rds")
df$log.ratio <- log(df$head / df$tail)

# Run model
post <- run.jags(
  data = list(
    n.blocks = length(unique(df$block)),
    n.sires = length(unique(df$sire)),
    n.dams = length(unique(df$dam)),
    n.interactions = length(unique(df$interaction)),
    n.larvae = nrow(df),
    block = as.numeric(factor(df$block)),
    sire = as.numeric(factor(df$sire)),
    dam = as.numeric(factor(df$dam)),
    interaction = as.numeric(factor(df$interaction)),
    log.ratio.range = range(df$log.ratio),
    log.ratio1 = df$log.ratio,
    log.ratio2 = df$log.ratio,
    log.ratio3 = df$log.ratio
  ),
  model = "model {
    # ---- overall mean ----
    mean.overall ~ dunif(log.ratio.range[1], log.ratio.range[2])
    var.overall ~ dunif(0, 1e5)
    
    # ---- priors for effect means ----
    for(b in 1:n.blocks) {
      # prior for mean of each block
      block.mean[b] ~ dunif(log.ratio.range[1], log.ratio.range[2])
      
      # prior of mean and variance of log ratio in each block for computing heritability
      mean.log.ratio.block[b] ~ dunif(log.ratio.range[1], log.ratio.range[2])
      var.log.ratio.block[b] ~ dunif(0, 1e3)
    }
    sire.mean ~ dnorm(0, 1e-3)
    dam.mean ~ dnorm(0, 1e-3)
    interaction.mean ~ dnorm(0, 1e-3)
      
    # ---- variance priors ----
    sire.var ~ dunif(0, 1e2)
    dam.var ~ dunif(0, 1e2)
    interaction.var ~ dunif(0, 1e2)
    resid.var ~ dunif(0, 1e2)
  
    # ---- sire priors ----
    for(s in 1:n.sires) {
      sire.eff[s] ~ dnorm(sire.mean, 1 / sire.var)
    }
    
    # ---- dam priors ----
    for(d in 1:n.dams) {
      dam.eff[d] ~ dmnorm.vcov(dam.mean, 1 / dam.var)
    }
    
    # ---- interaction priors ----
    for(i in 1:n.interactions) {
      interaction.eff[i] ~ dnorm(interaction.mean, 1 / interaction.var)
    }
    
    # ---- likelihood ----
    for(l in 1:n.larvae) {
      log.ratio1[l] ~ dnorm(mean.overall, 1 / var.overall)
      log.ratio2[l] ~ dnorm(mean.log.ratio.block[block[l]], 1 / var.log.ratio.block[block[l]]) 
      
      mu[l] <- block.mean[block[l]] + 
        sire.eff[sire[l]] + 
        dam.eff[dam[l]] +
        interaction.eff[interaction[l]]
      log.ratio3[l] ~ dnorm(mu[l], 1 / resid.var)
      
      # draw for posterior predictive check
      log.ratio.ppc[l] ~ dnorm(mu[l], 1 / resid.var)
    }
  }",
  monitor = c(
    "deviance", "sire.var", "dam.var", "interaction.var", "resid.var",
    "mean.overall", "mean.log.ratio.block", 'log.ratio.ppc[l]'
  ), 
  inits = function() list(
    .RNG.name = "lecuyer::RngStream",
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c("glm", "lecuyer"),
  n.chains = chains,
  adapt = adapt,
  burnin = burnin,
  sample = ceiling(total.sample / chains),
  thin = thin,
  method = "parallel"
)
end <- Sys.time()
elapsed <- swfscMisc::autoUnits(post$timetaken)

# Extract list of posteriors
p <- swfscMisc::runjags2list(post)

# Add QG metrics to list
p$VA <- 4 * p$sire.var
p$VM <- p$dam.var - p$sire.var
p$VD <- 4 * p$interaction.var
p$VP <- p$VA + p$VM + p$VD + p$resid.var

# Compute heritability and evolvability based on deVillemereuil et al 2016
convertVarScale <- function(metric, p) {
  lapply(1:length(p[[metric]]), function(i) {
    QGglmm::QGparams(
      mu = p$mean.overall[i],
      var.a = p[[metric]][i],
      var.p = p$VP[i],
      predict = p$mean.log.ratio.block[, i],
      custom.model = list(
        inv.link = function(x) {exp(x)},
        var.func = function(x) {0},
        d.inv.link = function(x) {exp(x)}
      ),
      verbose = FALSE
    )
  }) |> 
    bind_rows()
}
var.obs <- sapply(c('VA', 'VM', 'VD'), convertVarScale, p = p, simplify = FALSE)

p$H <- var.obs$VA$h2.obs
p$E <- var.obs$VA |> 
  mutate(E = var.a.obs / (p$mean.overall ^ 2)) |> 
  pull('E')

save.image(format(end, "3_Model_outputs/Model_II_posterior_%Y%m%d_%H%M.rdata"))

plot(post, file = format(end, "3_Model_outputs/Model_II_plots_%Y%m%d_%H%M.pdf"))

print(elapsed)
