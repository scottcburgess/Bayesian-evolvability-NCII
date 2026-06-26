rm(list = ls())
library(tidyverse)
library(runjags)
source('0_misc_funcs.R')

start.time <- Sys.time()

# ---- MCMC parameters --------------------------------------------------------
chains       <- 3  #50
adapt        <- 1000 #5000 
burnin       <- 1000  #50000
total.sample <- 1000
thin         <- 1 #1000


# Load data -------------------------------------------------------------------

trunk_tail.df  <- readRDS('Data/trunk_tail_data.rds')
hatch_settle.df <- readRDS('Data/hatch_settle_data.rds')


# Format model data -----------------------------------------------------------

blocks       <- c(trunk_tail.df$block,       hatch_settle.df$block)       |> unique() |> sort()
sires        <- c(trunk_tail.df$sire,        hatch_settle.df$sire)        |> unique() |> sort()
dams         <- c(trunk_tail.df$dam,         hatch_settle.df$dam)         |> unique() |> sort()
interactions <- c(trunk_tail.df$interaction, hatch_settle.df$interaction) |> unique() |> sort()

trunk_tail.df <- trunk_tail.df |>
  mutate(
    block       = as.numeric(factor(block,       levels = blocks)),
    sire        = as.numeric(factor(sire,        levels = sires)),
    dam         = as.numeric(factor(dam,         levels = dams)),
    interaction = as.numeric(factor(interaction, levels = interactions))
  )

hatch_settle.df <- hatch_settle.df |>
  mutate(
    block       = as.numeric(factor(block,       levels = blocks)),
    sire        = as.numeric(factor(sire,        levels = sires)),
    dam         = as.numeric(factor(dam,         levels = dams)),
    interaction = as.numeric(factor(interaction, levels = interactions))
  )

h.df  <- filter(hatch_settle.df, metric == 'hatching')
sh.df <- filter(hatch_settle.df, metric == 'settling')

# Wishart scale matrix: diagonal = 1/df so marginal variance ~ 1/Gamma(df/2, df/2)
# df = 5 (minimally informative: df must be >= p = 4)
wishart.df <- 5
# Scale matrix calibrated to match the original dunif(0, 1e3) variance prior.
# Under dwish(R, df) with df=5 and p=4, the marginal prior on each variance
# is InvGamma(shape=(df-p+1)/2, scale=R_ii/2) = InvGamma(1, R_ii/2).
# Median of InvGamma(1, b) = b/log(2).
# To match the dunif(0, 1e3) median of 500: R_ii = 2 * 500 * log(2) = 1e3*log(2).
wishart.R  <- diag(1e3 * log(2), 4)  # median variance ~500, matching dunif(0, 1e3)

model.data <- list(
  n.blocks  = length(blocks),
  n.sires   = length(sires),
  n.dams    = length(dams),
  n.int     = length(interactions),
  n.larvae  = nrow(trunk_tail.df),
  tt.block  = trunk_tail.df$block,
  tt.sire   = trunk_tail.df$sire,
  tt.dam    = trunk_tail.df$dam,
  tt.int    = trunk_tail.df$interaction,
  block.mean.range = cbind(
    range(trunk_tail.df$trunk),
    range(trunk_tail.df$tail),
    h.df  |> group_by(interaction) |> summarize(pct = mean(outcome), .groups = 'drop') |>
      filter(pct > 0 & pct < 1) |> pull('pct') |> qlogis() |> range(),
    sh.df |> group_by(interaction) |> summarize(pct = mean(outcome), .groups = 'drop') |>
      filter(pct > 0 & pct < 1) |> pull('pct') |> qlogis() |> range()
  ),
  length1   = cbind(Trunk = trunk_tail.df$trunk, Tail = trunk_tail.df$tail),
  length2   = cbind(Trunk = trunk_tail.df$trunk, Tail = trunk_tail.df$tail),
  n.hatch   = nrow(h.df),
  h.block   = h.df$block,
  h.sire    = h.df$sire,
  h.dam     = h.df$dam,
  h.int     = h.df$interaction,
  hatch     = h.df$outcome,
  n.settle_hatch = nrow(sh.df),
  sh.block  = sh.df$block,
  sh.sire   = sh.df$sire,
  sh.dam    = sh.df$dam,
  sh.int    = sh.df$interaction,
  settle_hatch = sh.df$outcome,
  # Wishart hyperparameters
  wishart.df = wishart.df,
  wishart.R  = wishart.R,
  # Residual Wishart (2x2 for trunk/tail only)
  resid.wishart.df = 3,                   # df >= p = 2
  # Same calibration: InvGamma(1, R_ii/2), median = R_ii/(2*log(2)) = 500
  # => R_ii = 1e3 * log(2)
  resid.wishart.R  = diag(1e3 * log(2), 2)
)


# Run Model -------------------------------------------------------------------

post <- run.jags(
  data = model.data,
  model = 'model {

    # ---- overall trunk/tail mean and variance (for evolvability) ----
    for(t in 1:2) {
      overall.mean[t] ~ dunif(block.mean.range[1, t], block.mean.range[2, t])
      overall.var[t]  ~ dunif(0, 1e5)
    }

    # ---- block means ----
    for(t in 1:4) {
      for(b in 1:n.blocks) {
        block.mean[b, t] ~ dunif(block.mean.range[1, t], block.mean.range[2, t])
      }
    }

    # ---- effect means ----
    for(t in 1:4) {
      sire.mean[t] ~ dnorm(0, 1e-3)
      dam.mean[t]  ~ dnorm(0, 1e-3)
      int.mean[t]  ~ dnorm(0, 1e-3)
    }

    sire.prec[1:4, 1:4]  ~ dwish(wishart.R[,], wishart.df)
    dam.prec[1:4, 1:4]   ~ dwish(wishart.R[,], wishart.df)
    int.prec[1:4, 1:4]   ~ dwish(wishart.R[,], wishart.df)
    resid.prec[1:2, 1:2] ~ dwish(resid.wishart.R[,], resid.wishart.df)

    sire.vcov[1:4, 1:4]  <- inverse(sire.prec[,])
    dam.vcov[1:4, 1:4]   <- inverse(dam.prec[,])
    int.vcov[1:4, 1:4]   <- inverse(int.prec[,])

    # resid.vcov is 4x4 with trunk/tail block in [1:2, 1:2] and zeros elsewhere,
    # matching the structure of sire/dam/int vcov for downstream QG calculations.
    resid.vcov[1:2, 1:2] <- inverse(resid.prec[,])
    for(t in 1:2) {
      resid.vcov[t, 3]     <- 0
      resid.vcov[t, 4]     <- 0
      resid.vcov[3, t]     <- 0
      resid.vcov[4, t]     <- 0
    }
    resid.vcov[3, 3] <- 0
    resid.vcov[3, 4] <- 0
    resid.vcov[4, 3] <- 0
    resid.vcov[4, 4] <- 0

    for(s in 1:n.sires) {
      sire.eff[s, 1:4] ~ dmnorm(sire.mean[], sire.prec[,])
    }

    for(d in 1:n.dams) {
      dam.eff[d, 1:4] ~ dmnorm(dam.mean[], dam.prec[,])
    }

    for(i in 1:n.int) {
      int.eff[i, 1:4] ~ dmnorm(int.mean[], int.prec[,])
    }

    # ---- trunk/tail likelihood ----
    for(l in 1:n.larvae) {
      for(t in 1:2) {
        # likelihood of overall mean for evolvability
        length1[l, t] ~ dnorm(overall.mean[t], 1 / overall.var[t])

        # linear predictor
        length.mu[l, t] <- block.mean[tt.block[l], t] +
          sire.eff[tt.sire[l], t] +
          dam.eff[tt.dam[l], t] +
          int.eff[tt.int[l], t]
      }
      length2[l, ]        ~ dmnorm(length.mu[l, ], resid.prec[,])
      length.ppd[l, 1:2]  ~ dmnorm(length.mu[l, ], resid.prec[,])
    }

    # ---- hatching likelihood ----
    for(h in 1:n.hatch) {
      hatch.p[h] <- ilogit(
        block.mean[h.block[h], 3] +
        sire.eff[h.sire[h], 3] +
        dam.eff[h.dam[h], 3] +
        int.eff[h.int[h], 3]
      )
      hatch[h]     ~ dbern(hatch.p[h])
      hatch.ppd[h] ~ dbern(hatch.p[h])
    }

    # ---- settling|hatching likelihood ----
    for(sh in 1:n.settle_hatch) {
      settle_hatch.p[sh] <- ilogit(
        block.mean[sh.block[sh], 4] +
        sire.eff[sh.sire[sh], 4] +
        dam.eff[sh.dam[sh], 4] +
        int.eff[sh.int[sh], 4]
      )
      settle_hatch[sh]     ~ dbern(settle_hatch.p[sh])
      settle_hatch.ppd[sh] ~ dbern(settle_hatch.p[sh])
    }
  }',
  monitor = c(
    'deviance',
    'sire.vcov', 'dam.vcov', 'int.vcov', 'resid.vcov',
    'overall.mean', 'block.mean',
    'length.ppd', 'hatch.ppd', 'settle_hatch.ppd',
    'sire.mean', 'dam.mean', 'int.mean',
    'sire.eff', 'dam.eff', 'int.eff'
  ),
  inits = function() {
    list(
      sire.prec  = diag(4),
      dam.prec   = diag(4),
      int.prec   = diag(4),
      resid.prec = diag(2),
      .RNG.name  = 'lecuyer::RngStream',
      .RNG.seed  = sample(1:9999, 1)
    )
  },
  modules   = c('glm', 'lecuyer'),
  n.chains  = chains,
  adapt     = adapt,
  burnin    = burnin,
  sample    = ceiling(total.sample / chains),
  thin      = thin,
  method    = 'parallel',
  summarise = FALSE
)
save.image('Model_outputs/last_posterior.rdata')
end.time <- Sys.time()


# Extract posterior to named list of arrays: p --------------------------------

p <- swfscMisc::runjags2list(post)

names.2 <- c('Trunk', 'Tail')
names.4 <- c(names.2, 'Hatch', 'Settle|Hatch')

dimnames(p$sire.vcov)[1:2] <-
  dimnames(p$dam.vcov)[1:2]  <-
  dimnames(p$int.vcov)[1:2]  <-
  dimnames(p$resid.vcov)[1:2] <- list(names.4, names.4)   # now 4x4 from model

dimnames(p$overall.mean)[[1]] <- names.2

dimnames(p$block.mean)[[2]]  <-
  dimnames(p$sire.mean)[[1]] <-
  dimnames(p$dam.mean)[[1]]  <-
  dimnames(p$int.mean)[[1]]  <-
  dimnames(p$sire.eff)[[2]]  <-
  dimnames(p$dam.eff)[[2]]   <-
  dimnames(p$int.eff)[[2]]   <- names.4

dimnames(p$length.ppd)[[2]] <- names.2
dimnames(p$block.mean)[[1]] <- blocks
dimnames(p$sire.eff)[[1]]   <- sires
dimnames(p$dam.eff)[[1]]    <- dams
dimnames(p$int.eff)[[1]]    <- interactions

# Set block means to NA for blocks not observed in each metric
for(m in dimnames(p$block.mean)[[2]]) {
  b <- switch(
    m,
    Hatch          = setdiff(1:model.data$n.blocks, model.data$h.block),
    'Settle|Hatch' = setdiff(1:model.data$n.blocks, model.data$sh.block),
    setdiff(1:model.data$n.blocks, model.data$tt.block)
  )
  p$block.mean[b, m, ] <- NA
}

p <- addQGmetrics(p)


# QGmvparams ------------------------------------------------------------------

qgparams.post <- parallel::mclapply(1:dim(p$VA)[3], function(i) {
  VA <- p$VA[, , i]
  VP <- p$VP[, , i]
  lapply(1:dim(p$block.mean)[1], function(b) {
    mu          <- p$block.mean[b, , i]
    not.missing <- which(!is.na(mu))
    qg <- QGglmm::QGmvparams(
      mu      = mu[not.missing],
      vcv.G   = VA[not.missing, not.missing],
      vcv.P   = VP[not.missing, not.missing],
      models  = c('Gaussian', 'Gaussian', 'binom1.logit', 'binom1.logit')[not.missing],
      verbose = FALSE
    )
    mu[not.missing] <- qg$mean.obs
    va <- VA; va[not.missing, not.missing] <- qg$vcv.G.obs
    vp <- VP; vp[not.missing, not.missing] <- qg$vcv.P.obs
    list(mu.obs = mu, va.obs = va, vp.obs = vp)
  }) |>
    list_transpose(simplify = FALSE)
}, mc.cores = chains) |>
  list_transpose(simplify = FALSE)

save.image('Model_outputs/last_qgparams.rdata')


# Format QGmvparams arrays ----------------------------------------------------

block.mean.obs <- do.call(
  abind::abind,
  c(lapply(qgparams.post$mu.obs, abind::abind, along = 2), list(along = 3))
)
dimnames(block.mean.obs)[[1]] <- names.4
dimnames(block.mean.obs)[[2]] <- blocks

va.obs <- do.call(
  abind::abind,
  c(lapply(qgparams.post$va.obs, abind::abind, along = 3), list(along = 4))
)
dimnames(va.obs)[1:2] <- list(names.4, names.4)
dimnames(va.obs)[[3]] <- blocks

vp.obs <- do.call(
  abind::abind,
  c(lapply(qgparams.post$vp.obs, abind::abind, along = 3), list(along = 4))
)
dimnames(vp.obs)[1:2] <- list(names.4, names.4)
dimnames(vp.obs)[[3]] <- blocks


# Summarize QGmvparams posteriors ---------------------------------------------

block.mean.obs.smry <- block.mean.obs |> apply(c(1, 2), vecSmry) |> aperm(c(2, 1, 3))
va.obs.smry <- va.obs |> apply(c(1, 2), vecSmry) |> aperm(c(3, 2, 1))
vp.obs.smry <- va.obs |> apply(c(1, 2), vecSmry) |> aperm(c(3, 2, 1))


# Evolvability ----------------------------------------------------------------

e.params_means <- do.call(
  rbind,
  parallel::mclapply(1:dim(p$VA)[3], function(i) {
    evolvability::evolvabilityMeans(
      G     = as.vector(p$VA[names.2, names.2, i]),
      means = p$overall.mean[, i]
    )
  }, mc.cores = 14)
)

e.params_BetaMCMC <- evolvability::evolvabilityBetaMCMC(
  G_mcmc = evolvability::meanStdGMCMC(
    t(apply(p$VA[names.2, names.2, ], 3, as.vector)),
    t(p$overall.mean)
  ),
  Beta      = evolvability::randomBeta(1000, 2),
  post.dist = TRUE
)

B <- matrix(
  c(
    c(0, 1),
    c(-(1/sqrt(2)), -(1/sqrt(2))),
    c( (1/sqrt(2)), -(1/sqrt(2)))
  ),
  nrow = 2, ncol = 3
)

e.params_beta <- do.call(
  rbind,
  parallel::mclapply(1:dim(p$VA)[3], function(i) {
    do.call(rbind, lapply(1:ncol(B), function(j) {
      tmp <- evolvability::evolvabilityBeta(
        G     = p$VA[names.2, names.2, i],
        Beta  = B[, j],
        means = p$overall.mean[, i]
      )
      data.frame(sample = i, Beta_index = j,
                 e = tmp$e, r = tmp$r, c = tmp$c, a = tmp$a, i = tmp$i)
    }))
  }, mc.cores = 14)
)
rownames(e.params_beta) <- NULL


# Compute the selection differentials from the covariance -----------------------

sg.df <- expand.grid(
  z = c('Trunk', 'Tail'), 
  W = c('Hatch', 'Settle|Hatch'),
  stringsAsFactors = FALSE
) |> 
  mutate(z.W = paste0(z, ' : ', W))

sg <- lapply(1:nrow(sg.df), function(i) {
  sg.i <- sapply(1:dim(block.mean.obs)[2], function(b) {
    # delta z = cov(z,W)/mean(W) = R = sg (in units of microns)
    va.obs[sg.df$z[i], sg.df$W[i], b, ] / block.mean.obs[sg.df$W[i], b, ]
  }) |> 
    t() 
  
  sg.i |> 
    as.data.frame() |> 
    setNames(1:ncol(sg.i)) |> 
    mutate(block = blocks) |> 
    pivot_longer(-block, names_to = 'sample', values_to = 'sg') |> 
    mutate(
      z = sg.df$z[i],
      W = sg.df$W[i],
      z.W = sg.df$z.W[i]
    ) |>
    left_join(
      p$block.mean[, sg.df$z[i], ] |> 
        as.data.frame() |> 
        setNames(1:dim(p$block.mean)[3]) |> 
        mutate(block = blocks) |> 
        pivot_longer(-block, names_to = 'sample', values_to = 'block.mean'),
      by = c('block', 'sample')
    ) |> 
    # delta z / mean(z) = sg / mean(z) (in units of percent)
    mutate(sg.pct = 100 * sg / block.mean)
}) |> 
  bind_rows()

# Add total response to selection
sg <- sg |> 
  bind_rows(
    sg |> 
      group_by(block, sample, z) |> 
      summarize(
        sg = sum(sg),
        block.mean = mean(block.mean),
        sg.pct = sum(sg.pct), 
        .groups = 'drop'
      ) |> 
      mutate(
        W = 'Total',
        z.W = paste0(z, ' : ', W)
      )
  )

# genetic selection gradient
# beta_g = G^-1 * s_g
traits <- unique(sg$z)
beta_g <- lapply(unique(sg$W), function(w) {
  lapply(dimnames(va.obs)[[3]], function(b) {
    sg.w.b <- sg |> 
      filter(block == b, W == w) |> 
      select(sample, z, sg) |> 
      pivot_wider(names_from = 'z', values_from = 'sg') |> 
      arrange(sample) |> 
      select(-sample) |> 
      as.matrix()
    
    sapply(1:dim(sg.w.b)[1], function(i) {
      if(any(is.na(sg.w.b[i, traits]))) return(setNames(c(NA, NA), traits))
      solve(va.obs[traits, traits, b, i], sg.w.b[i, traits])
    }) |> 
      t() |> 
      as.data.frame() |> 
      mutate(W = w, block = b) |> 
      select(W, block, everything())
  })
}) |> 
  bind_rows()



# CODA summary ------------------------------------------------------------

post.smry <- smrzPost(
  post, 
  c('deviance', 'sire.vcov', 'dam.vcov', 'int.vcov', 'resid.vcov')
)


# Posterior Predictive Check ----------------------------------------------

ppc <- bind_rows(
  data.frame(metric = rep('Trunk', nrow(model.data$length1))),
  data.frame(metric = rep('Tail', nrow(model.data$length1))),
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
    'Settle|Hatch' = list(
      obs = model.data$settle_hatch[id], 
      ppd = p$settle_hatch.ppd[id, ]
    ),
    list(
      obs = model.data$length1[id, m],
      ppd = p$length.ppd[id, m, ]
    )
  )
  
  c(
    pct.gte.obs = mean(x$obs >= x$ppd, na.rm = TRUE), 
    mean.diff = mean(x$obs - x$ppd, na.rm = TRUE)
  )
}) |>
  t()
)

ppc.smry <- smrzPPC(ppc)


# Save all objects --------------------------------------------------------

end.time <- if(exists('end.time')) end.time else Sys.time()
post.file <- format(end.time, 'posterior_%Y%m%d_%H%M.rdata', tz = 'GMT')
save.image(file.path('Model_outputs', post.file))


# Plot posterior distributions --------------------------------------------

plot(
  post,
  vars = c('deviance', 'sire.vcov', 'dam.vcov', 'int.vcov', 'resid.vcov'),
  file = format(end.time, 'Model_outputs/plots_%Y%m%d_%H%M.pdf', tz = 'GMT')
)


# Plot diagnostics --------------------------------------------------------

pdf(format(end.time, "Model_outputs/diagnostics_%Y%m%d_%H%M.pdf", tz = 'GMT'))

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


rmarkdown::render(
  'posterior_summary.Rmd',
  params = list(posterior_file = post.file),
  output_file = format(end.time, 'Figures and Tables/posterior_summary_%Y%m%d_%H%M.html', tz = 'GMT')
)


cat(
  'Run start: ', format(start.time, tz = 'GMT'), '\n',
  'Run end: ', format(end.time, tz = 'GMT'), '\n',
  'Model elapsed: ', format(swfscMisc::autoUnits(post$timetaken)), '\n',
  'Run elapsed: ', format(difftime(end.time, start.time)), '\n',
  sep = ''
)
