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


# hatch_settle.df <- hatch_settle.df |> 
#   # filter hatching/settling data for blocks and interactions to use 
#   filter(interaction %in% interactions & block %in% blocks) |> 
#   mutate(
#     block = as.numeric(factor(block)),
#     sire = as.numeric(factor(sire))
#   ) |> 
#   # compress to number of successes and trials for hatching and settling by interaction
#   group_by(block, sire) |> 
#   summarize(
#     n.hatched.trial = sum(metric == 'hatching'),
#     n.hatched = sum(metric == 'hatching' & outcome == 1),
#     pct.hatched = n.hatched / n.hatched.trial,
#     n.settled.hatched.trial = sum(metric == 'settling'),
#     n.settled.hatched = sum(metric == 'settling' & outcome == 1),
#     pct.settled.hatched= n.settled.hatched / n.settled.hatched.trial,
#     pr.settle = pct.hatched * pct.settled.hatched,
#     .groups = 'drop'
#   )

hatch_settle.df <- hatch_settle.df |> 
  filter(interaction %in% interactions & block %in% blocks) |> 
  mutate(
    block = as.numeric(factor(block)),
    sire = as.numeric(factor(sire)),
    dam = as.numeric(factor(dam)),
    interaction = as.numeric(factor(interaction))
  )
    
h.df <- filter(hatch_settle.df, metric == 'hatching') 
sh.df <- filter(hatch_settle.df, metric == 'settling')

# Run model
post <- run.jags(
  data = list(
    n.blocks = n_distinct(trunk_tail.df$block),
    n.sires = n_distinct(trunk_tail.df$sire),
    n.dams = n_distinct(trunk_tail.df$dam),
    n.int = n_distinct(trunk_tail.df$interaction),
    n.larvae = nrow(trunk_tail.df),
    tt.block = trunk_tail.df$block,
    tt.sire = trunk_tail.df$sire,
    tt.dam = trunk_tail.df$dam,
    tt.int = trunk_tail.df$interaction,
    block.mean.range = cbind(
      round(range(trunk_tail.df$trunk)), 
      round(range(trunk_tail.df$tail)),
      qlogis(c(0.2, 0.95))
    ),
    length1 = cbind(trunk_tail.df$trunk, trunk_tail.df$tail),
    length2 = cbind(trunk_tail.df$trunk, trunk_tail.df$tail),
    length3 = cbind(trunk_tail.df$trunk, trunk_tail.df$tail),
    n.hatch = nrow(h.df),
    h.block = h.df$block,
    h.sire = h.df$sire,
    h.dam = h.df$dam,
    h.int = h.df$interaction,
    hatch = h.df$outcome,
    n.settle_hatch = nrow(sh.df),
    sh.block = sh.df$block,
    sh.sire = sh.df$sire,
    sh.dam = sh.df$dam,
    sh.int = sh.df$interaction,
    settle_hatch = sh.df$outcome
  ),
  model = 'model {
    # ---- trunk and tail specific priors ----
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
      
      resid.vcov[t, t] ~ dunif(0, 1e3)
    }
    resid.corr ~ dunif(-1, 1)
    
    # ---- hatch and settle|hatch specific priors ----
    for(t in 3:4) {
      for(b in 1:n.blocks) {
        block.mean[b, t] ~ dunif(block.mean.range[1, 3], block.mean.range[2, 3])
      }
    }
    
    # ---- priors for effects ----
    for(t in 1:4) {
      # ---- effect magnitudes ----
      sire.mean[t] ~ dnorm(0, 1e-5)
      dam.mean[t] ~ dnorm(0, 1e-5)
      int.mean[t] ~ dnorm(0, 1e-5)
      
      # ---- variances ----
      sire.vcov[t, t] ~ dunif(0, 1e3)
      dam.vcov[t, t] ~ dunif(0, 1e3)
      int.vcov[t, t] ~ dunif(0, 1e3)
    }
    
    # ---- construct variance/covariance matrices ----
    for(t1 in 1:3) {
      for(t2 in (t1+1):4) {
        # ---- correlations  ----
        sire.corr[t1, t2] ~ dunif(-1, 1)
        sire.corr[t2, t1] <- sire.corr[t1, t2]
        dam.corr[t1, t2] ~ dunif(-1, 1)
        dam.corr[t2, t1] <- dam.corr[t1, t2]
        int.corr[t1, t2] ~ dunif(-1, 1)
        int.corr[t2, t1] <- int.corr[t1, t2]
      
        # ---- covariances ----
        sire.vcov[t1, t2] <- sire.corr[t1, t2] * sqrt(sire.vcov[t1, t1] * sire.vcov[t2, t2])
        sire.vcov[t2, t1] <- sire.vcov[t1, t2]
        dam.vcov[t1, t2]  <- dam.corr[t1, t2] * sqrt(dam.vcov[t1, t1] * dam.vcov[t2, t2])
        dam.vcov[t2, t1] <- dam.vcov[t1, t2]
        int.vcov[t1, t2] <- int.corr[t1, t2] * sqrt(int.vcov[t1, t1] * int.vcov[t2, t2])
        int.vcov[t2, t1] <- int.vcov[t1, t2]
      }
    }
    resid.vcov[1, 2] <- resid.corr * sqrt(resid.vcov[1, 1] * resid.vcov[2, 2])
    resid.vcov[2, 1] <- resid.vcov[1, 2]
    
    # ---- prior for additive sire effect (for each sire) ----
    for(s in 1:n.sires) {
      sire.eff[s, 1:4] ~ dmnorm.vcov(sire.mean, sire.vcov)
    }
    
    # ---- prior for maternal effect (for each dam) ----
    for(d in 1:n.dams) {
      dam.eff[d, 1:4] ~ dmnorm.vcov(dam.mean, dam.vcov)
    }
    
    # ---- prior for interaction effect (for each sire x dam interaction) ----
    for(int in 1:n.int) {
      int.eff[int, 1:4] ~ dmnorm.vcov(int.mean, int.vcov)
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
          int.eff[tt.int[l], t]
      }
      # likelihood of l-th larvae for both traits from multivariate normal
      length3[l, ] ~ dmnorm.vcov(length.mu[l, 1:2], resid.vcov[1:2, 1:2])
      
      # draw for posterior predictive check
      length.ppd[l, 1:2] ~ dmnorm.vcov(length.mu[l, 1:2], resid.vcov[1:2, 1:2])
    }
    
    # ---- hatching likelihood ----
    for(h in 1:n.hatch) {
      hatch.mu[h] <- block.mean[h.block[h], 3] + 
        sire.eff[h.sire[h], 3] + 
        dam.eff[h.dam[h], 3] +
        int.eff[h.int[h], 3]
      hatch[h] ~ dbern(ilogit(hatch.mu[h]))
    }
    
    # ---- settling|hatching likelihood ----
    for(sh in 1:n.settle_hatch) {
      settle_hatch.mu[sh] <- block.mean[sh.block[sh], 4] + 
        sire.eff[sh.sire[sh], 4] + 
        dam.eff[sh.dam[sh], 4] +
        int.eff[sh.int[sh], 4]
      settle_hatch[sh] ~ dbern(ilogit(settle_hatch.mu[sh]))
    }
  }',
  monitor = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'int.vcov', 
    'resid.vcov', 'overall.block.mean', 'mean.overall'
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
  dimnames(p$mean.overall)[[1]] <- c('Trunk', 'Tail')
dimnames(p$resid.vcov)[1:2] <- list(c('Trunk', 'Tail'), c('Trunk', 'Tail'))
dimnames(p$sire.vcov) <-
  dimnames(p$dam.vcov) <-
  dimnames(p$int.vcov) <- list(
    c('Trunk', 'Tail', 'Hatch', 'Settle|Hatch'), 
    c('Trunk', 'Tail', 'Hatch', 'Settle|Hatch')
  )

# Add QG metrics to list
p <- addQGmetrics(p)

# 
# # Compute heritability and evolvability based on deVillemereuil et al 2016
# qgparams.post <- sapply(dimnames(p$pr.overall)[[1]], function(m) {
#   parallel::mclapply(1:dim(p$VA)[3], function(i) {
#     QGglmm::QGparams(
#       var.a = p$VA[m, m, i],
#       var.p = p$VP[m, m, i],
#       predict = qlogis(p$pr.block[, m, i]),
#       model = 'binom1.logit',
#       verbose = FALSE
#     )
#   }, mc.cores = 14) |> 
#     bind_rows() |> 
#     mutate(E = var.a.obs / (p$pr.overall[m, ] ^ 2))
# }, simplify = FALSE)
# 
# p$H <- t(sapply(qgparams.post, function(x) x$h2.obs))
# p$E <- t(sapply(qgparams.post, function(x) x$E))
# 
# 
# beta <- sapply(1:dim(p$pr.settle.sire.eff)[2], function(i) {
#   s_g <- 16 * c(
#     cov(p$sire.eff[, 'Trunk', i], p$pr.settle.sire.eff[, i]),
#     cov(p$sire.eff[, 'Tail', i], p$pr.settle.sire.eff[, i])
#   )
#   inv.G <- solve(p$VA[, , i])
#   beta <- inv.G %*% s_g
#   setNames(
#     c(s_g, beta[, 1]),
#     c('s_g.trunk.settle', 's_g.tail.settle', 'beta.trunk', 'beta.tail')
#   )
# }) |> 
#   t()
# 
# 
# CODA summary ------------------------------------------------------------

post.smry <- smrzPost(
  post, 
  c('deviance', 'sire.vcov', 'dam.vcov', 'int.vcov', 'resid.vcov')
)
# 
# 
# # Posterior Predictive Check ----------------------------------------------
# 
# length.obs <- cbind(Trunk = trunk_tail.df$trunk, Tail = trunk_tail.df$tail)
# 
# ppc <- expand_grid(
#   metric = colnames(length.obs),
#   id = 1:nrow(length.obs)
# ) |> 
#   mutate(metric = factor(metric, colnames(length.obs))) |> 
#   bind_rows(
#     data.frame(metric = 'Hatching', id = 1:nrow(hatch_settle.df)),
#     data.frame(metric = 'Settling.Hatching', id = 1:nrow(hatch_settle.df))
#   )
# 
# ppc <- ppc |> 
#   cbind(sapply(1:nrow(ppc), function(i) {
#     id <- ppc$id[i]
#     
#     obs <- switch(
#       ppc$metric[id],
#       Hatching = hatch_settle.df$n.hatched[id],
#       Settling.Hatching = hatch_settle.df$n.settled.hatched[id],
#       length.obs[id, ppc$metric[id]]
#     )
#     
#     ppd <- switch(
#       ppc$metric[id],
#       Hatching = p$hatch.ppd[id, ],
#       Settling.Hatching = p$settle.hatch.ppd[id, ],
#       p$length.ppd[id, ppc$metric[id], ]
#     )
#     
#     c(pct.gte.obs = mean(obs >= ppd), mean.diff = mean(obs - ppd))
#   }) |> 
#     t()
#   )
# 
# ppc.smry <- smrzPPC(ppc)
# 
# 
# 
# # Save all objects
# save.image(format(end, 'Model_outputs/Model_III_posterior_%Y%m%d_%H%M.rdata'))
# 
# 
# Plot posterior distributions
plot(
  post,
  vars = c('deviance', 'sire.vcov', 'dam.vcov', 'int.vcov', 'resid.vcov'),
  file = format(end, 'Model_outputs/Model_III_plots_%Y%m%d_%H%M.pdf')
)
# 
# 
# # Plot diagnostics
# pdf(format(end, "Model_outputs/Model_III_diagnostics_%Y%m%d_%H%M.pdf"))
# 
# ggplot(post.smry$post) +
#   geom_histogram(aes(values), bins = 20) +
#   facet_wrap(~diag, scales = 'free_x')
# 
# ggplot(ppc) +
#   geom_histogram(aes(pct.gte.obs), binwidth = 0.05) +
#   geom_vline(aes(xintercept = median.pct), data = ppc.smry, color = 'red') +
#   geom_vline(aes(xintercept = lower.pct), data = ppc.smry, linetype = 'dashed', color = 'red') +
#   geom_vline(aes(xintercept = upper.pct), data = ppc.smry, linetype = 'dashed', color = 'red') +
#   facet_wrap(~ metric, scales = 'free') +
#   labs(x = 'Percent of PPD >= Observed', y = 'Count')
# 
# ggplot(ppc) +
#   geom_histogram(aes(mean.diff), bins = 50) +
#   geom_vline(aes(xintercept = median.diff), data = ppc.smry, color = 'red') +
#   geom_vline(aes(xintercept = lower.diff), data = ppc.smry, linetype = 'dashed', color = 'red') +
#   geom_vline(aes(xintercept = upper.diff), data = ppc.smry, linetype = 'dashed', color = 'red') +
#   facet_wrap(~ metric, scales = 'free') +
#   labs(x = 'Metric Difference (Observed - PPD)', y = 'Count')
# 
# dev.off()
# 
# 
# print(elapsed)