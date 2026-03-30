rm(list = ls())
library(tidyverse)
library(runjags)
source('0_misc_funcs.R')

# MCMC parameters
chains <- 50
adapt <- 100
burnin <- 80000
total.sample <- 8000
thin <- 20000

# Load data
df <- readRDS('Data/hatch_settle_data.rds') 

# Only keep blocks with both hatching and settling data
blocks.to.keep <- table(block = df$block, metric = df$metric) |> 
  as.data.frame() |> 
  filter(Freq > 0) |> 
  group_by(block) |> 
  summarize(n = n(), .groups = 'drop') |> 
  filter(n == 2) |> 
  pull(block) |> 
  as.character() |> 
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
  model = 'model {
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
      sire.mean[m] ~ dnorm(0, 1e-6)
      dam.mean[m] ~ dnorm(0, 1e-6)
      interaction.mean[m] ~ dnorm(0, 1e-6)
      
      # ---- prior for variances ----
      sire.vcov[m, m] ~ dunif(0, 1000)
      dam.vcov[m, m] ~ dunif(0, 1000)
      interaction.vcov[m, m] ~ dunif(0, 1000)
    }
    
    # ---- prior for sire covariance ---
    sire.corr ~ dunif(-1, 1)
    sire.vcov[1, 2] <- sire.corr * sqrt(sire.vcov[1, 1] * sire.vcov[2, 2])
    sire.vcov[2, 1] <- sire.vcov[1, 2]
    
    # ---- prior for dam covariance ----
    dam.corr ~ dunif(-1, 1)
    dam.vcov[1, 2] <- dam.corr * sqrt(dam.vcov[1, 1] * dam.vcov[2, 2])
    dam.vcov[2, 1] <- dam.vcov[1, 2]
    
    # ---- prior for interaction covariance ----
    interaction.corr ~ dunif(-1, 1)
    interaction.vcov[1, 2] <- interaction.corr * sqrt(interaction.vcov[1, 1] * interaction.vcov[2, 2])
    interaction.vcov[2, 1] <- interaction.vcov[1, 2]
    
    # ---- prior for sire effect (for each sire) ----
    for(s in 1:n.sires) {
      sire.eff[s, 1:2] ~ dmnorm.vcov(sire.mean, sire.vcov)
    }
    
    # ---- priors for dam effect (for each dam) ----
    for(d in 1:n.dams) {
      dam.eff[d, 1:2] ~ dmnorm.vcov(dam.mean, dam.vcov)
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
         sire.eff[sire[l], metric[l]] + 
         dam.eff[dam[l], metric[l]] + 
         interaction.eff[interaction[l], metric[l]] 
      outcome3[l] ~ dbern(pr[l])
      
      # draw for posterior predictive check
      outcome.ppd[l] ~ dbern(pr[l])
    }
  }',
  monitor = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov',
    'pr.overall', 'pr.block', 'outcome.ppd'
  ), 
  inits = function() list(
    .RNG.name = 'lecuyer::RngStream',
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c('glm', 'lecuyer'),
  summarise = FALSE,
  n.chains = chains,
  adapt = adapt,
  burnin = burnin,
  sample = ceiling(total.sample / chains),
  thin = thin,
  method = 'parallel'
)
end <- Sys.time()
elapsed <- swfscMisc::autoUnits(post$timetaken)

# Extract posterior to list of arrays - p
p <- swfscMisc::runjags2list(post)
dimnames(p$pr.overall)[[1]] <- 
  dimnames(p$pr.block)[[2]] <- c('Hatching', 'Settling')
dimnames(p$sire.vcov)[1:2] <-
  dimnames(p$dam.vcov)[1:2] <-
  dimnames(p$interaction.vcov)[1:2] <-
  list(dimnames(p$pr.overall)[[1]], dimnames(p$pr.overall)[[1]])

# Add QG metrics to list
p <- addQGmetrics(p)

# Compute heritability and evolvability based on deVillemereuil et al 2016
qgparams.post <- sapply(dimnames(p$pr.overall)[[1]], function(m) {
  parallel::mclapply(1:dim(p$VA)[3], function(i) {
    QGglmm::QGparams(
      var.a = p$VA[m, m, i],
      var.p = p$VP[m, m, i],
      predict = qlogis(p$pr.block[, m, i]),
      model = 'binom1.logit',
      verbose = FALSE
    )
  }, mc.cores = 14) |> 
    bind_rows() |> 
    mutate(E = var.a.obs / (p$pr.overall[m, ] ^ 2))
}, simplify = FALSE)

p$H <- t(sapply(qgparams.post, function(x) x$h2.obs))
p$E <- t(sapply(qgparams.post, function(x) x$E))


vcv.obs <- sapply(
  c('VA', 'VM', 'VD'), 
  convertVCVscale.II, 
  p = p, 
  simplify = FALSE
)


# Compute evolvability 
e.params_BetaMCMC <- evolvability::evolvabilityBetaMCMC(
  G_mcmc = evolvability::meanStdGMCMC(
    t(apply(vcv.obs$VA$vcv.G.obs, 3, as.vector)),
    t(p$pr.overall)
  ),
  Beta = evolvability::randomBeta(1000, 2),
  post.dist = TRUE
)


# CODA summary ------------------------------------------------------------

post.smry <- smrzPost(post, c(
  'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov',
  'pr.overall', 'pr.block'
))


# Posterior Predictive Check ----------------------------------------------

ppc <- data.frame(id = 1:nrow(df))

ppc$pct.gte.obs <- sapply(1:nrow(ppc), function(i) {
  obs <- df$outcome[i]
  ppd <- p$outcome.ppd[i, ]
  mean(obs >= ppd)
})

ppc$mean.diff <- sapply(1:nrow(ppc), function(i) {
  obs <- df$outcome[i]
  ppd <- p$outcome.ppd[i, ]
  mean(obs - ppd)
})

ppc.smry <- smrzPPC(ppc)


# Save all objects
save.image(format(end, 'Model_outputs/Model_II_posterior_%Y%m%d_%H%M.rdata'))


# Plot posterior distributions
plot(
  post, 
  vars = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov',
    'pr.overall', 'pr.block'
  ),
  file = format(end, 'Model_outputs/Model_II_plots_%Y%m%d_%H%M.pdf')
)


# Plot diagnostics
pdf(format(end, "Model_outputs/Model_II_diagnostics_%Y%m%d_%H%M.pdf"))

ggplot(post.smry) +
  geom_histogram(aes(values), bins = 20) +
  facet_wrap(~diag, scales = 'free_x')

ggplot(ppc) +
  geom_histogram(aes(pct.gte.obs), binwidth = 0.05) +
  geom_vline(aes(xintercept = median.pct), data = ppc.smry, color = 'red') +
  geom_vline(aes(xintercept = lower.pct), data = ppc.smry, linetype = 'dashed', color = 'red') +
  geom_vline(aes(xintercept = upper.pct), data = ppc.smry, linetype = 'dashed', color = 'red') +
  labs(x = 'Percent of PPD >= Observed', y = 'Count')

ggplot(ppc) +
  geom_histogram(aes(mean.diff), bins = 50) +
  geom_vline(aes(xintercept = median.diff), data = ppc.smry, color = 'red') +
  geom_vline(aes(xintercept = lower.diff), data = ppc.smry, linetype = 'dashed', color = 'red') +
  geom_vline(aes(xintercept = upper.diff), data = ppc.smry, linetype = 'dashed', color = 'red') +
  labs(x = 'Pr(Settling) Difference (Observed - PPD)', y = 'Count')

dev.off()


print(elapsed)