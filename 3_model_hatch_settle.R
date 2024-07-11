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
df <- readRDS("hatch_settle_data.rds") 

# Only keep blocks with both hatching and settling data
blocks.to.keep <- table(block = df$block, metric = df$metric) %>% 
  as.data.frame() %>% 
  filter(Freq > 0) %>% 
  group_by(block) %>% 
  summarize(n = n(), .groups = "drop") %>% 
  filter(n == 2) %>% 
  pull(block) %>% 
  as.character() %>% 
  as.integer()
df <- filter(df, block %in% blocks.to.keep)

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
    metric = as.numeric(factor(df$metric)),
    outcome1 = df$outcome,
    outcome2 = df$outcome,
    outcome3 = df$outcome
  ),
  model = "model {
    for(m in 1:2) { 
      # ---- prior for overall probability of outcome for computing evolvability ----
      pr.overall[m] ~ dunif(0, 1)    
    
      # ---- prior for block and effect means ----
      for(b in 1:n.blocks) {
        block.mean.pr[b, m] ~ dunif(0, 1)
        block.mean[b, m] <- logit(block.mean.pr[b, m])
        
        # pior for block probability of outcome for computing heritability
        pr.block[b, m] ~ dunif(0, 1)
      }
      additive.mean[m] ~ dnorm(0, 1e-6)
      maternal.mean[m] ~ dnorm(0, 1e-6)
      interaction.mean[m] ~ dnorm(0, 1e-6)
      
      # ---- prior for variances ----
      additive.vcov[m, m] ~ dunif(0, 1000)
      maternal.vcov[m, m] ~ dunif(0, 1000)
      interaction.vcov[m, m] ~ dunif(0, 1000)
    }
    
    # ---- prior for additive covariance ---
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
      outcome1[l] ~ dbern(pr.overall[metric[l]])
      outcome2[l] ~ dbern(pr.block[block[l], metric[l]])
      
      logit(pr[l]) <- block.mean[block[l], metric[l]] + 
         additive.sire.eff[sire[l], metric[l]] + 
         additive.dam.eff[dam[l], metric[l]] +
         maternal.eff[dam[l], metric[l]] + 
         interaction.eff[interaction[l], metric[l]] 
      outcome3[l] ~ dbern(pr[l])
    }
  }",
  monitor = c(
    "deviance", "additive.vcov", "maternal.vcov", "interaction.vcov",
    "pr.overall", "pr.block"
  ), 
  inits = function() list(
    .RNG.name = "lecuyer::RngStream",
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c("glm", "lecuyer"),
  summarise = FALSE,
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
dimnames(p$pr.overall)[[1]] <- dimnames(p$pr.block)[[2]] <- sort(unique(df$metric))
dimnames(p$additive.vcov)[1:2] <-
  dimnames(p$maternal.vcov)[1:2] <-
  dimnames(p$interaction.vcov)[1:2] <-
  list(dimnames(p$pr.overall)[[1]], dimnames(p$pr.overall)[[1]])

# Add QG metrics to list
p$VA <- 4 * p$additive.vcov
p$VM <- p$maternal.vcov
p$VD <- 4 * p$interaction.vcov
# p$VP <- p$VA + p$VM + p$VD
p$VP <- 2 * p$additive.vcov + p$VM + p$interaction.vcov

# Compute heritability and evolvability based on deVillemereuil et al 2016
qgparams.post <- sapply(dimnames(p$pr.overall)[[1]], function(m) {
  parallel::mclapply(1:dim(p$VA)[3], function(i) {
    QGglmm::QGparams(
      var.a = p$VA[m, m, i],
      var.p = p$VP[m, m, i],
      predict = swfscMisc::logOdds(p$pr.block[, m, i]),
      model = "binom1.logit",
      verbose = FALSE
    )
  }, mc.cores = 14) |> 
    bind_rows() |> 
    mutate(E = var.a.obs / (p$pr.overall[m, ] ^ 2))
}, simplify = FALSE)

p$H <- t(sapply(qgparams.post, function(x) x$h2.obs))
p$E <- t(sapply(qgparams.post, function(x) x$E))


e.params_BetaMCMC <- evolvability::evolvabilityBetaMCMC(
  G_mcmc = t(apply(p$VA, 3, as.vector)),
  Beta = evolvability::randomBeta(1000, 2),
  post.dist = TRUE
)

e.params_BetaMCMC <- evolvability::evolvabilityBetaMCMC(
  G_mcmc = evolvability::meanStdGMCMC(
    t(apply(p$VA, 3, as.vector)),
    t(p$pr.overall)
  ),
  Beta = evolvability::randomBeta(1000, 2),
  post.dist = TRUE
)

# Save all objects and plot posterior summaries
save.image(format(end, "hatch_settle_posterior_%Y%m%d_%H%M.rdata"))

plot(
  post, 
  file = format(end, "hatch_settle_plots_%Y%m%d_%H%M.pdf")
)

print(elapsed)