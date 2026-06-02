rm(list = ls())
library(tidyverse)
library(runjags)
source('0_misc_funcs.R')

start <- Sys.time()

# ---- MCMC parameters
chains <- 6 #50
adapt <- 100
burnin <- 500 #500000
total.sample <- 1000 #5000 
thin <- 1 #5000


# Load data ---------------------------------------------------------------

trunk_tail.df <- readRDS('Data/trunk_tail_data.rds') 
hatch_settle.df <- readRDS('Data/hatch_settle_data.rds') 


# Format model data -------------------------------------------------------

blocks <- c(trunk_tail.df$block, hatch_settle.df$block) |> 
  unique() |> 
  sort()
sires <- c(trunk_tail.df$sire, hatch_settle.df$sire) |> 
  unique() |> 
  sort()
dams <- c(trunk_tail.df$dam, hatch_settle.df$dam) |> 
  unique() |> 
  sort()
interactions <- c(trunk_tail.df$interaction, hatch_settle.df$interaction) |> 
  unique() |> 
  sort()

trunk_tail.df <- trunk_tail.df |>  
  mutate(
    block = as.numeric(factor(block, levels = blocks)),
    sire = as.numeric(factor(sire, levels = sires)),
    dam = as.numeric(factor(dam, levels = dams)),
    interaction = as.numeric(factor(interaction, levels = interactions))
  )

hatch_settle.df <- hatch_settle.df |> 
  mutate(
    block = as.numeric(factor(block, levels = blocks)),
    sire = as.numeric(factor(sire, levels = sires)),
    dam = as.numeric(factor(dam, levels = dams)),
    interaction = as.numeric(factor(interaction, levels = interactions))
  )

h.df <- filter(hatch_settle.df, metric == 'hatching') 
sh.df <- filter(hatch_settle.df, metric == 'settling')

model.data <- list(
  n.blocks = length(blocks),
  n.sires = length(sires),
  n.dams = length(dams),
  n.int = length(interactions),
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
  length = cbind(Trunk = trunk_tail.df$trunk, Tail = trunk_tail.df$tail),
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
)


# ---- Run model

post <- run.jags(
  data = model.data,
  model = 'model {
    # ---- trunk and tail specific priors ----
    for(t in 1:2) {  
      # ---- priors for block and effect means ----
      for(b in 1:n.blocks) {
        block.mean[b, t] ~ dunif(block.mean.range[1, t], block.mean.range[2, t])
      }
    }
    
    # ---- hatch and settle|hatch specific priors ----
    for(t in 3:4) {
      for(b in 1:n.blocks) {
        block.mean[b, t] ~ dunif(block.mean.range[1, 3], block.mean.range[2, 3])
      }
    }
    
    # ---- priors for effects ----
    for(t in 1:4) {
      # ---- effect means ----
      sire.mean[t] ~ dnorm(0, 1e-5)
      dam.mean[t] ~ dnorm(0, 1e-5)
      int.mean[t] ~ dnorm(0, 1e-5)
      
      # ---- variances ----
      sire.vcov[t, t] ~ dunif(0, 1e3)
      dam.vcov[t, t] ~ dunif(0, 1e3)
      int.vcov[t, t] ~ dunif(0, 1e3)
      resid.vcov[t, t] ~ dunif(0, 1e3)
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
        resid.corr[t1, t2] ~ dunif(-1, 1)
        resid.corr[t2, t1] <- resid.corr[t1, t2]
      
        # ---- covariances ----
        sire.vcov[t1, t2] <- sire.corr[t1, t2] * sqrt(sire.vcov[t1, t1] * sire.vcov[t2, t2])
        sire.vcov[t2, t1] <- sire.vcov[t1, t2]
        dam.vcov[t1, t2]  <- dam.corr[t1, t2] * sqrt(dam.vcov[t1, t1] * dam.vcov[t2, t2])
        dam.vcov[t2, t1] <- dam.vcov[t1, t2]
        int.vcov[t1, t2] <- int.corr[t1, t2] * sqrt(int.vcov[t1, t1] * int.vcov[t2, t2])
        int.vcov[t2, t1] <- int.vcov[t1, t2]
        resid.vcov[t1, t2] <- resid.corr[t1, t2] * sqrt(resid.vcov[t1, t1] * resid.vcov[t2, t2])
        resid.vcov[t2, t1] <- resid.vcov[t1, t2]
      }
    }
    
    # ---- draw additive sire effect ----
    for(s in 1:n.sires) {
      sire.eff[s, 1:4] ~ dmnorm.vcov(sire.mean, sire.vcov)
    }
    
    # ---- draw maternal effect (for each dam) ----
    for(d in 1:n.dams) {
      dam.eff[d, 1:4] ~ dmnorm.vcov(dam.mean, dam.vcov)
    }
    
    # ---- draw interaction effect (for each sire x dam interaction) ----
    for(int in 1:n.int) {
      int.eff[int, 1:4] ~ dmnorm.vcov(int.mean, int.vcov)
    }
    
    # ---- trunk/tail ----
    for(l in 1:n.larvae) {
      for(t in 1:2) {        
        # expected mean for the l-th larvae and t-th trait
        length.mu[l, t] <- block.mean[tt.block[l], t] + 
          sire.eff[tt.sire[l], t] + 
          dam.eff[tt.dam[l], t] +
          int.eff[tt.int[l], t]
      }
      # likelihood of l-th larvae for both traits from multivariate normal
      length[l, ] ~ dmnorm.vcov(length.mu[l, ], resid.vcov[1:2, 1:2])
      # draw for posterior predictive check
      length.ppd[l, 1:2] ~ dmnorm.vcov(length.mu[l, ], resid.vcov[1:2, 1:2])
    }
    
    # ---- hatching ----
    for(h in 1:n.hatch) {
      # expected probability of hatching 
      hatch.p[h] <- ilogit(
        block.mean[h.block[h], 3] + 
        sire.eff[h.sire[h], 3] + 
        dam.eff[h.dam[h], 3] +
        int.eff[h.int[h], 3]
      )
      # likelihood of hatching
      hatch[h] ~ dbern(hatch.p[h])
      # draw for posterior predictive check
      hatch.ppd[h] ~ dbern(hatch.p[h])
    }
    
    # ---- settling|hatching ----
    for(sh in 1:n.settle_hatch) {
      # expected probability of settling|hatching
      settle_hatch.p[sh] <- ilogit(
        block.mean[sh.block[sh], 4] + 
        sire.eff[sh.sire[sh], 4] + 
        dam.eff[sh.dam[sh], 4] +
        int.eff[sh.int[sh], 4]
      )
      # likelihod of settling|hatching
      settle_hatch[sh] ~ dbern(settle_hatch.p[sh])
      # draw for posterior predictive check
      settle_hatch.ppd[sh] ~ dbern(settle_hatch.p[sh])
    }
  }',
  monitor = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'int.vcov', 'resid.vcov', 
    'block.mean', 'length.ppd', 'hatch.ppd', 'settle_hatch.ppd'
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


# Extract posterior to named list of arrays: p -----------------------------

p <- swfscMisc::runjags2list(post)

names.2 <- c('Trunk', 'Tail')
names.4 <- c(names.2, 'Hatch', 'Settle|Hatch')
dimnames(p$sire.vcov) <-
  dimnames(p$dam.vcov) <-
  dimnames(p$int.vcov) <- 
  dimnames(p$resid.vcov) <- list(names.4, names.4)
dimnames(p$block.mean)[[2]] <- names.4
dimnames(p$length.ppd)[[2]] <- names.2

# add QG metrics to list
p$resid.vcov[c('Hatch', 'Settle|Hatch'), , ] <- 0
p$resid.vcov[, c('Hatch', 'Settle|Hatch'), ] <- 0
for(m in dimnames(p$block.mean)[[2]]) {
  b <- switch(
    m, 
    Hatch = setdiff(1:model.data$n.blocks, model.data$h.block),
    'Settle|Hatch' = setdiff(1:model.data$n.blocks, model.data$sh.block),
    setdiff(1:model.data$n.blocks, model.data$tt.block)
  )
  p$block.mean[b, m, ] <- NA
}
p <- addQGmetrics(p)


# QGmvparams: compute posterior on observed scale -------------------------

# iterate over every posterior sample
qgparams.post <- parallel::mclapply(1:dim(p$VA)[3], function(i) { 
  # extract VA and VP matrices for this sample
  VA <- p$VA[, , i]
  VP <- p$VP[, , i]
  # iterate over blocks
  lapply(1:dim(p$block.mean)[1], function(b) {
    mu <- p$block.mean[b, , i]
    va <- VA
    vp <- VP
    # identify metrics without this block
    not.missing <- which(!is.na(mu))
    qg <- QGglmm::QGmvparams(
      mu = mu[not.missing],
      vcv.G = VA[not.missing, not.missing],
      vcv.P = VP[not.missing, not.missing],
      models = c('Gaussian', 'Gaussian', 'binom1.logit', 'binom1.logit')[not.missing],
      verbose = FALSE
    ) 
    # reload results to original vectors/matrices to preserve NAs
    mu[not.missing] <- qg$mean.obs
    va[not.missing, not.missing] <- qg$vcv.G.obs
    vp[not.missing, not.missing] <- qg$vcv.P.obs
    list(mu.obs = mu, va.obs = va, vp.obs = vp)
  }) |> 
    list_transpose(simplify = FALSE)
}, mc.cores = chains) |> 
  list_transpose(simplify = FALSE)

# format 3D array of block means [metric, block, sample]
block.mean.obs <- do.call(
  abind::abind,
  c(
    lapply(qgparams.post$mu.obs, abind::abind, along = 2), 
    list(along = 3)
  )
)
dimnames(block.mean.obs)[[1]] <- names.4

# format 4D array of VA matrices [metric, metric, block, sample]
va.obs <- do.call(
  abind::abind,
  c(
    lapply(qgparams.post$va.obs, abind::abind, along = 3),
    list(along = 4)
  )
)
dimnames(va.obs)[1:2] <- list(names.4, names.4)

# format 4D array of VP matrices [metric, metric, block, sample]
vp.obs <- do.call(
  abind::abind,
  c(
    lapply(qgparams.post$vp.obs, abind::abind, along = 3),
    list(along = 4)
  )
)
dimnames(vp.obs)[1:2] <- list(names.4, names.4)


# Summarize QGmvparams posteriors -----------------------------------------

# median and HDI of mean.obs for each block
block.mean.obs.smry <- block.mean.obs |> 
  apply(c(1, 2), function(x) {
    c(median = median(x, na.rm = TRUE), HDInterval::hdi(x))
  }) |> 
  aperm(c(3, 2, 1))

# median and HDI of VA across blocks
va.obs.smry <- va.obs |> 
  apply(c(1, 2), function(x) {
    x <- as.vector(x)
    c(median = median(x, na.rm = TRUE), HDInterval::hdi(x))
  }) |> 
  aperm(c(3, 2, 1))

# median and HDI of VP across blocks
vp.obs.smry <- va.obs |> 
  apply(c(1, 2), function(x) {
    x <- as.vector(x)
    c(median = median(x, na.rm = TRUE), HDInterval::hdi(x))
  }) |> 
  aperm(c(3, 2, 1))


# CODA summary ------------------------------------------------------------

post.smry <- smrzPost(
  post, 
  c('deviance', 'sire.vcov', 'dam.vcov', 'int.vcov', 'resid.vcov')
)


# Posterior Predictive Check ----------------------------------------------

ppc <- bind_rows(
  data.frame(metric = rep('Trunk', nrow(model.data$length))),
  data.frame(metric = rep('Tail', nrow(model.data$length))),
  data.frame(metric = rep('Hatch', length(model.data$hatch))),
  data.frame(metric = rep('Settle|Hatch', length(model.data$settle_hatch)))
) |> 
  mutate(id = 1:n(), .by = metric)

ppc <- cbind(ppc, sapply(1:nrow(ppc), function(i) {
  m <- ppc$metric[i]
  id <- ppc$id[i]
  
  x <- switch(
    m,
    Hatch = list(obs = model.data$hatch[id], ppd = p$hatch.ppd[id, ]),
    'Settle|Hatch' = list(obs = model.data$settle_hatch[id], ppd = p$settle_hatch.ppd[id, ]),
    list(
      obs = model.data$length[id, m],
      ppd = p$length.ppd[id, m, ]
    )
  )
  
  c(pct.gte.obs = mean(x$obs >= x$ppd), mean.diff = mean(x$obs - x$ppd))
}) |>
  t()
)

ppc.smry <- smrzPPC(ppc)


# Save all objects --------------------------------------------------------

end <- Sys.time()
save.image(format(end, 'Model_outputs/Model_III_posterior_%Y%m%d_%H%M.rdata'))


# Plot posterior distributions --------------------------------------------

plot(
  post,
  vars = c('deviance', 'sire.vcov', 'dam.vcov', 'int.vcov', 'resid.vcov'),
  file = format(end, 'Model_outputs/Model_II_plots_%Y%m%d_%H%M.pdf')
)


# Plot diagnostics --------------------------------------------------------

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


cat('Run start:', format(start))
cat('Run end:', format(end))
cat('Model elapsed:', format(swfscMisc::autoUnits(post$timetaken)))
cat('Run elapsed: ', format(difftime(end, start)))
