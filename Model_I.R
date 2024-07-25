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
df <- readRDS("head_tail_data.rds") 

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
    length.range = cbind(round(range(df$head)), round(range(df$tail))),
    length1 = cbind(df$head, df$tail),
    length2 = cbind(df$head, df$tail)
  ),
  model = "model {
    # for each t-trait...
    for(t in 1:2) {  
      # ---- prior for overall mean and variance ----
      mean.overall[t] ~ dunif(length.range[1, t], length.range[2, t])
      var.overall[t] ~ dunif(0, 1e5)
      
      # ---- prior for block and effect means ----
      for(b in 1:n.blocks) {
        block.mean[b, t] ~ dunif(length.range[1, t], length.range[2, t])
      }
      additive.mean[t] ~ dnorm(0, 1e-5)
      maternal.mean[t] ~ dnorm(0, 1e-5)
      interaction.mean[t] ~ dnorm(0, 1e-5)
      
      # ---- prior for variances ----
      additive.vcov[t, t] ~ dunif(0, 1e3)
      maternal.vcov[t, t] ~ dunif(0, 1e3)
      interaction.vcov[t, t] ~ dunif(0, 1e3)
      resid.vcov[t, t] ~ dunif(0, 1e3)
    }
  
    # ---- prior for additive covariance ----
    additive.corr ~ dunif(-1, 1)
    additive.vcov[1, 2] <- additive.corr * sqrt(additive.vcov[1, 1] * additive.vcov[2, 2])
    additive.vcov[2, 1] <- additive.vcov[1, 2]
    
    # ---- prior maternal covariance ----
    maternal.corr ~ dunif(-1, 1)
    maternal.vcov[1, 2] <- maternal.corr * sqrt(maternal.vcov[1, 1] * maternal.vcov[2, 2])
    maternal.vcov[2, 1] <- maternal.vcov[1, 2]
    
    # ---- prior for interaction covariance ----
    interaction.corr ~ dunif(-1, 1)
    interaction.vcov[1, 2] <- interaction.corr * sqrt(interaction.vcov[1, 1] * interaction.vcov[2, 2])
    interaction.vcov[2, 1] <- interaction.vcov[1, 2]
    
    # ---- prior for residual covariance ----
    resid.corr ~ dunif(-1, 1)
    resid.vcov[1, 2] <- resid.corr * sqrt(resid.vcov[1, 1] * resid.vcov[2, 2])
    resid.vcov[2, 1] <- resid.vcov[1, 2]
    
    # ---- prior for additive sire effect (for each sire) ----
    for(s in 1:n.sires) {
      additive.sire.eff[s, 1:2] ~ dmnorm.vcov(additive.mean, additive.vcov)
    }
    
    # ---- priors for additive dam & maternal effect (for each dam) ----
    for(d in 1:n.dams) {
      additive.dam.eff[d, 1:2] ~ dmnorm.vcov(additive.mean, additive.vcov)
      maternal.eff[d, 1:2] ~ dmnorm.vcov(maternal.mean, maternal.vcov)
    }
    
    # ---- prior for interaction effect (for each sire x dam interaction) ----
    for(i in 1:n.interactions) {
      interaction.eff[i, 1:2] ~ dmnorm.vcov(interaction.mean, interaction.vcov)
    }
    
    # ---- likelihood ----
    for(l in 1:n.larvae) {
      for(t in 1:2) {
        # likelihood of overall mean for computing evolvability
        length1[l, t] ~ dnorm(mean.overall[t], 1 / var.overall[t])
        
        # expected mean for the l-th larvae and t-th trait
        mu[l, t] <- block.mean[block[l], t] + 
          additive.sire.eff[sire[l], t] + 
          additive.dam.eff[dam[l], t] + 
          maternal.eff[dam[l], t] +
          interaction.eff[interaction[l], t]
      }
      # likelihood of l-th larvae for both traits from multivariate normal
      length2[l, ] ~ dmnorm.vcov(mu[l, ], resid.vcov)
    }
  }",
  monitor = c(
    "deviance", "additive.vcov", "maternal.vcov", "interaction.vcov", 
    "resid.vcov", "mean.overall"
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

# Extract posterior to list of arrays - p
p <- swfscMisc::runjags2list(post)
rownames(p$mean.overall) <- c("head", "tail")
dimnames(p$additive.vcov)[1:2] <- 
  dimnames(p$maternal.vcov)[1:2] <- 
  dimnames(p$interaction.vcov)[1:2] <- 
  dimnames(p$resid.vcov)[1:2] <- 
  list(rownames(p$mean.overall), rownames(p$mean.overall))

# Add QG metrics to list
p$VA <- 4 * p$additive.vcov
p$VM <- p$maternal.vcov
p$VD <- 4 * p$interaction.vcov
p$VP <- 2 * p$additive.vcov + p$VM + p$interaction.vcov + p$resid.vcov
p$H <- p$VA / p$VP
p$E <- rbind(
  head = p$VA[1, 1, ] / (p$mean.overall[1, ] ^ 2),
  tail = p$VA[2, 2, ] / (p$mean.overall[2, ] ^ 2)
)

# Compute evolvability
e.params_means <- do.call(
  rbind,
  parallel::mclapply(1:dim(p$VA)[3], function(i) {
    evolvability::evolvabilityMeans(
      G = as.vector(p$VA[, , i]),
      means = p$mean.overall[, i]
    )
  }, mc.cores = 14) 
)
  
e.params_BetaMCMC <- evolvability::evolvabilityBetaMCMC(
  G_mcmc = evolvability::meanStdGMCMC(
    t(apply(p$VA, 3, as.vector)),
    t(p$mean.overall)
  ),
  Beta = evolvability::randomBeta(1000, 2),
  post.dist = TRUE
)

# Save all objects and plot posterior summaries
save.image(format(end, "Model_I_posterior_%Y%m%d_%H%M.rdata"))

plot(
  post, 
  file = format(end, "Model_I_plots_%Y%m%d_%H%M.pdf")
)

print(elapsed)