rm(list = ls())
library(tidyverse)
library(rstan)
source('0_misc_funcs.R')

# Compile Stan model once; reuse across runs
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

start.time <- Sys.time()


# ---- MCMC parameters --------------------------------------------------------

chains    <- 6
warmup    <- 1000   # equivalent to adapt + burnin in JAGS
iter      <- warmup + ceiling(1000 / chains) * chains  # warmup + kept samples
thin      <- 1
adapt_delta      <- 0.95   # increase if divergences appear
max_treedepth    <- 12


# Load data -------------------------------------------------------------------

trunk_tail.df  <- readRDS('Data/trunk_tail_data.rds')
hatch_settle.df <- readRDS('Data/hatch_settle_data.rds')


# Format model data -----------------------------------------------------------

blocks       <- sort(unique(c(trunk_tail.df$block,       hatch_settle.df$block)))
sires        <- sort(unique(c(trunk_tail.df$sire,        hatch_settle.df$sire)))
dams         <- sort(unique(c(trunk_tail.df$dam,         hatch_settle.df$dam)))
interactions <- sort(unique(c(trunk_tail.df$interaction, hatch_settle.df$interaction)))

trunk_tail.df <- trunk_tail.df |>
  mutate(
    block       = as.integer(factor(block,       levels = blocks)),
    sire        = as.integer(factor(sire,        levels = sires)),
    dam         = as.integer(factor(dam,         levels = dams)),
    interaction = as.integer(factor(interaction, levels = interactions))
  )

hatch_settle.df <- hatch_settle.df |>
  mutate(
    block       = as.integer(factor(block,       levels = blocks)),
    sire        = as.integer(factor(sire,        levels = sires)),
    dam         = as.integer(factor(dam,         levels = dams)),
    interaction = as.integer(factor(interaction, levels = interactions))
  )

h.df  <- filter(hatch_settle.df, metric == 'hatching')
sh.df <- filter(hatch_settle.df, metric == 'settling')

# Block-mean range matrix [2 x 4]: rows = (min, max), cols = Trunk/Tail/Hatch/Settle
block_mean_range <- cbind(
  range(trunk_tail.df$trunk),
  range(trunk_tail.df$tail),
  h.df  |> group_by(interaction) |> summarize(pct = mean(outcome), .groups = 'drop') |>
    filter(pct > 0 & pct < 1) |> pull(pct) |> qlogis() |> range(),
  sh.df |> group_by(interaction) |> summarize(pct = mean(outcome), .groups = 'drop') |>
    filter(pct > 0 & pct < 1) |> pull(pct) |> qlogis() |> range()
)

stan.data <- list(
  # Dimensions
  n_blocks  = length(blocks),
  n_sires   = length(sires),
  n_dams    = length(dams),
  n_int     = length(interactions),

  # Trunk/tail
  n_larvae  = nrow(trunk_tail.df),
  tt_block  = trunk_tail.df$block,
  tt_sire   = trunk_tail.df$sire,
  tt_dam    = trunk_tail.df$dam,
  tt_int    = trunk_tail.df$interaction,
  length_obs = cbind(trunk_tail.df$trunk, trunk_tail.df$tail),

  # Hatching
  n_hatch   = nrow(h.df),
  h_block   = h.df$block,
  h_sire    = h.df$sire,
  h_dam     = h.df$dam,
  h_int     = h.df$interaction,
  hatch     = h.df$outcome,

  # Settling|Hatching
  n_settle_hatch = nrow(sh.df),
  sh_block  = sh.df$block,
  sh_sire   = sh.df$sire,
  sh_dam    = sh.df$dam,
  sh_int    = sh.df$interaction,
  settle_hatch = sh.df$outcome,

  # Priors
  block_mean_range = block_mean_range
)


# Compile Stan model ----------------------------------------------------------

stan.model <- stan_model('larvae_model.stan')


# Run model -------------------------------------------------------------------

post <- sampling(
  object  = stan.model,
  data    = stan.data,
  chains  = chains,
  iter    = iter,
  warmup  = warmup,
  thin    = thin,
  control = list(
    adapt_delta   = adapt_delta,
    max_treedepth = max_treedepth
  ),
  seed    = 12345,
  pars    = c(
    'overall_mean', 'overall_var',
    'block_mean',
    'sire_mean', 'dam_mean', 'int_mean',
    'sire_vcov', 'dam_vcov', 'int_vcov', 'resid_vcov',
    'sire_eff', 'dam_eff', 'int_eff',
    'length_ppd', 'hatch_ppd', 'settle_hatch_ppd'
  )
)

save(post, file = 'Model_outputs/last_posterior_stan.rdata')
end.time <- Sys.time()


# Extract posterior to named list of arrays: p --------------------------------

# Helper: extract a named parameter as an array [dim1, ..., n_samples]
extract_par <- function(fit, par) {
  e <- rstan::extract(fit, pars = par, permuted = TRUE)[[par]]
  e  # returned as [n_samples, ...]
}

names.2 <- c('Trunk', 'Tail')
names.4 <- c(names.2, 'Hatch', 'Settle|Hatch')

p <- list()

# Overall mean/var [2 x n_samples]
p$overall.mean <- t(extract_par(post, 'overall_mean'))
rownames(p$overall.mean) <- names.2

# Block means [n_blocks x 4 x n_samples]
bm_raw <- extract_par(post, 'block_mean')   # [n_samples, n_blocks, 4]
p$block.mean <- aperm(bm_raw, c(2, 3, 1))
dimnames(p$block.mean) <- list(blocks, names.4, NULL)

# Effect means [4 x n_samples]
for (nm in c('sire', 'dam', 'int')) {
  raw <- extract_par(post, paste0(nm, '_mean'))  # [n_samples, 4]
  p[[paste0(nm, '.mean')]] <- t(raw)
  rownames(p[[paste0(nm, '.mean')]]) <- names.4
}

# Variance-covariance matrices [4 x 4 x n_samples]
for (nm in c('sire', 'dam', 'int')) {
  raw <- extract_par(post, paste0(nm, '_vcov'))   # [n_samples, 4, 4]
  p[[paste0(nm, '.vcov')]] <- aperm(raw, c(2, 3, 1))
  dimnames(p[[paste0(nm, '.vcov')]])[1:2] <- list(names.4, names.4)
}

# Residual vcov [2 x 2 x n_samples]  — expand to 4x4 with zeros for downstream QG calcs
rv_raw <- extract_par(post, 'resid_vcov')   # [n_samples, 2, 2]
rv_full <- array(0, dim = c(4, 4, dim(rv_raw)[1]))
rv_full[1:2, 1:2, ] <- aperm(rv_raw, c(2, 3, 1))
dimnames(rv_full)[1:2] <- list(names.4, names.4)
p$resid.vcov <- rv_full

# Random effects [n_levels x 4 x n_samples]
for (nm in c('sire', 'dam', 'int')) {
  raw <- extract_par(post, paste0(nm, '_eff'))   # [n_samples, n_levels, 4]
  p[[paste0(nm, '.eff')]] <- aperm(raw, c(2, 3, 1))
  dimnames(p[[paste0(nm, '.eff')]])[2] <- list(names.4)
}
dimnames(p$sire.eff)[[1]] <- sires
dimnames(p$dam.eff)[[1]]  <- dams
dimnames(p$int.eff)[[1]]  <- interactions

# PPD arrays
lppd_raw  <- extract_par(post, 'length_ppd')        # [n_samples, n_larvae, 2]
p$length.ppd <- aperm(lppd_raw, c(2, 3, 1))
dimnames(p$length.ppd)[[2]] <- names.2

hppd_raw  <- extract_par(post, 'hatch_ppd')         # [n_samples, n_hatch]
p$hatch.ppd <- t(hppd_raw)

shppd_raw <- extract_par(post, 'settle_hatch_ppd')  # [n_samples, n_settle_hatch]
p$settle_hatch.ppd <- t(shppd_raw)

# Set block.mean to NA for blocks not observed in each metric
for (m in dimnames(p$block.mean)[[2]]) {
  b <- switch(
    m,
    Hatch         = setdiff(seq_len(stan.data$n_blocks), stan.data$h_block),
    'Settle|Hatch' = setdiff(seq_len(stan.data$n_blocks), stan.data$sh_block),
    setdiff(seq_len(stan.data$n_blocks), stan.data$tt_block)
  )
  p$block.mean[b, m, ] <- NA
}

# Add QG metrics (VA, VP, etc.) — same function as before
p <- addQGmetrics(p)


# QGmvparams: compute posterior on observed scale ----------------------------

qgparams.post <- parallel::mclapply(1:dim(p$VA)[3], function(i) {
  VA <- p$VA[, , i]
  VP <- p$VP[, , i]
  lapply(seq_len(dim(p$block.mean)[1]), function(b) {
    mu <- p$block.mean[b, , i]
    not.missing <- which(!is.na(mu))
    qg <- QGglmm::QGmvparams(
      mu      = mu[not.missing],
      vcv.G   = VA[not.missing, not.missing],
      vcv.P   = VP[not.missing, not.missing],
      models  = c('Gaussian', 'Gaussian', 'binom1.logit', 'binom1.logit')[not.missing],
      verbose = FALSE
    )
    mu_full <- rep(NA_real_, 4);  names(mu_full) <- names.4
    mu_full[not.missing] <- qg$mean.obs
    va <- VA * 0
    va[not.missing, not.missing] <- qg$vcv.G.obs
    vp <- VP * 0
    vp[not.missing, not.missing] <- qg$vcv.P.obs
    list(mu.obs = mu_full, va.obs = va, vp.obs = vp)
  }) |>
    purrr::list_transpose(simplify = FALSE)
}, mc.cores = chains) |>
  purrr::list_transpose(simplify = FALSE)

save.image('Model_outputs/last_qgparams_stan.rdata')


# Diagnostics ----------------------------------------------------------------

# Stan-native diagnostics (divergences, Rhat, n_eff)
print(
  post,
  pars = c('sire_vcov', 'dam_vcov', 'int_vcov', 'resid_vcov',
           'overall_mean', 'sire_mean', 'dam_mean', 'int_mean'),
  digits = 3
)

check_hmc_diagnostics(post)


# Posterior Predictive Check -------------------------------------------------

ppc <- bind_rows(
  data.frame(metric = rep('Trunk',        nrow(stan.data$length_obs))),
  data.frame(metric = rep('Tail',         nrow(stan.data$length_obs))),
  data.frame(metric = rep('Hatch',        stan.data$n_hatch)),
  data.frame(metric = rep('Settle|Hatch', stan.data$n_settle_hatch))
) |>
  mutate(id = row_number(), .by = metric)

ppc <- cbind(ppc, t(sapply(seq_len(nrow(ppc)), function(i) {
  m  <- ppc$metric[i]
  id <- ppc$id[i]
  x <- switch(
    m,
    Hatch          = list(obs = stan.data$hatch[id],        ppd = p$hatch.ppd[id, ]),
    'Settle|Hatch' = list(obs = stan.data$settle_hatch[id], ppd = p$settle_hatch.ppd[id, ]),
    list(obs = stan.data$length_obs[id, match(m, names.2)], ppd = p$length.ppd[id, m, ])
  )
  c(pct.gte.obs = mean(x$obs >= x$ppd, na.rm = TRUE),
    mean.diff   = mean(x$obs  - x$ppd, na.rm = TRUE))
})))

ppc.smry <- smrzPPC(ppc)


# Save final workspace -------------------------------------------------------

save.image(format(end.time, 'Model_outputs/posterior_stan_%Y%m%d_%H%M.rdata', tz = 'GMT'))


# Trace / density plots (bayesplot or base rstan) ----------------------------

pdf(format(end.time, 'Model_outputs/diagnostics_stan_%Y%m%d_%H%M.pdf', tz = 'GMT'))

stan_trace(post, pars = c('sire_vcov', 'dam_vcov', 'int_vcov', 'resid_vcov'))
stan_dens( post, pars = c('sire_vcov', 'dam_vcov', 'int_vcov', 'resid_vcov'),
           separate_chains = TRUE)

ggplot(ppc) +
  geom_histogram(aes(pct.gte.obs), binwidth = 0.05) +
  geom_vline(aes(xintercept = median.pct),  data = ppc.smry, color = 'red') +
  geom_vline(aes(xintercept = lower.pct),   data = ppc.smry, linetype = 'dashed', color = 'red') +
  geom_vline(aes(xintercept = upper.pct),   data = ppc.smry, linetype = 'dashed', color = 'red') +
  facet_wrap(~ metric, scales = 'free') +
  labs(x = 'Percent of PPD \u2265 Observed', y = 'Count')

ggplot(ppc) +
  geom_histogram(aes(mean.diff), bins = 50) +
  geom_vline(aes(xintercept = median.diff), data = ppc.smry, color = 'red') +
  geom_vline(aes(xintercept = lower.diff),  data = ppc.smry, linetype = 'dashed', color = 'red') +
  geom_vline(aes(xintercept = upper.diff),  data = ppc.smry, linetype = 'dashed', color = 'red') +
  facet_wrap(~ metric, scales = 'free') +
  labs(x = 'Metric Difference (Observed \u2212 PPD)', y = 'Count')

dev.off()


cat(
  'Run start:  ', format(start.time, tz = 'GMT'), '\n',
  'Run end:    ', format(end.time,   tz = 'GMT'), '\n',
  'Run elapsed:', format(difftime(end.time, start.time)), '\n',
  sep = ''
)
