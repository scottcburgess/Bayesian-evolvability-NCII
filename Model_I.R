rm(list = ls())
library(tidyverse)
library(runjags)

# MCMC parameters
chains <- 10
adapt <- 100
burnin <- 100000
total.sample <- 10000 
thin <- 1000

# Load data
df <- readRDS("Data/trunk_tail_data.rds") 

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
    length.range = cbind(round(range(df$trunk)), round(range(df$tail))),
    length1 = cbind(df$trunk, df$tail),
    length2 = cbind(df$trunk, df$tail),
    length3 = cbind(df$trunk, df$tail)
  ),
  model = "model {
    # for each t-trait...
    for(t in 1:2) {  
      # ---- prior for overall mean and variance ----
      mean.overall[t] ~ dunif(length.range[1, t], length.range[2, t])
      var.overall[t] ~ dunif(0, 1e5)
      
      # ---- prior for block and effect means ----
      for(b in 1:n.blocks) {
        overall.block.mean[b, t] ~ dunif(length.range[1, t], length.range[2, t])
        overall.block.var[b, t] ~ dunif(0, 1e5)
        block.mean[b, t] ~ dunif(length.range[1, t], length.range[2, t])
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
  
    # ---- prior for sire covariance ----
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
    
    # ---- prior for residual covariance ----
    resid.corr ~ dunif(-1, 1)
    resid.vcov[1, 2] <- resid.corr * sqrt(resid.vcov[1, 1] * resid.vcov[2, 2])
    resid.vcov[2, 1] <- resid.vcov[1, 2]
    
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
      for(t in 1:2) {
        # likelihood of overall mean for computing evolvability
        length1[l, t] ~ dnorm(mean.overall[t], 1 / var.overall[t])
        length2[l, t] ~ dnorm(overall.block.mean[block[l], t], 1 / overall.block.var[block[l], t])
        
        # expected mean for the l-th larvae and t-th trait
        mu[l, t] <- block.mean[block[l], t] + 
          sire.eff[sire[l], t] + 
          dam.eff[dam[l], t] +
          interaction.eff[interaction[l], t]
      }
      # likelihood of l-th larvae for both traits from multivariate normal
      length3[l, ] ~ dmnorm.vcov(mu[l, ], resid.vcov)
      
      # draw of posterior predictive distribution
      length.ppd[l, 1:2] ~ dmnorm.vcov(mu[l, ], resid.vcov)
    }
  }",
  monitor = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov', 
    'resid.vcov', 'mean.overall', 'overall.block.mean', 'length.ppd'
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
rownames(p$mean.overall) <- 
  dimnames(p$overall.block.mean)[[2]] <- 
  dimnames(p$length.ppd)[[2]] <- c("Trunk", "Tail")
dimnames(p$sire.vcov)[1:2] <- 
  dimnames(p$dam.vcov)[1:2] <- 
  dimnames(p$interaction.vcov)[1:2] <- 
  dimnames(p$resid.vcov)[1:2] <- 
  list(rownames(p$mean.overall), rownames(p$mean.overall))

# Add QG metrics to list
p$VA <- 4 * p$sire.vcov
p$VM <- p$dam.vcov - p$sire.vcov
p$VD <- 4 * p$interaction.vcov
p$VP <- p$VA + p$VM + p$VD + p$resid.vcov
p$H <- p$VA / p$VP
p$E <- rbind(
  trunk = p$VA[1, 1, ] / (p$mean.overall[1, ] ^ 2),
  tail = p$VA[2, 2, ] / (p$mean.overall[2, ] ^ 2)
)

# Calculate average evolvability parameters of the G-matrix
e.params_means <- do.call(
  rbind,
  parallel::mclapply(1:dim(p$VA)[3], function(i) {
    evolvability::evolvabilityMeans(
      G = as.vector(p$VA[, , i]),
      means = p$mean.overall[, i]
    )
  }, mc.cores = 14) 
)

# Calculate posterior distribution of evolvability parameters 
# from a random set of selection gradients  
e.params_BetaMCMC <- evolvability::evolvabilityBetaMCMC(
  G_mcmc = evolvability::meanStdGMCMC(
    t(apply(p$VA, 3, as.vector)),
    t(p$mean.overall)
  ),
  Beta = evolvability::randomBeta(1000, 2),
  post.dist = TRUE
)

# Calculate evolvability parameters 
# along a specific set of selection gradients
B <- matrix(
  c(
    c(0, 1), # strong selection for long tails only, 
    c(-1, -1), # strong selection for short trunks and short tails
    c(1, -1) # strong selection for large trunks and small tails
  ), 
  nrow = 2, 
  ncol = 3
)

e.params_beta <- do.call(
  rbind,
  parallel::mclapply(1:dim(p$VA)[3], function(i) {
    do.call(
      rbind,
      lapply(1:ncol(B), function(j) {
        tmp <- evolvability::evolvabilityBeta(
          G = p$VA[, , i],
          Beta = B[, j],
          means = p$mean.overall[, i]
        )
        data.frame(
          sample = i,           
          Beta_index = j,       
          e = tmp$e,
          r = tmp$r,
          c = tmp$c,
          a = tmp$a,
          i = tmp$i
        )
      })
    )
  }, mc.cores = 14)
)
rownames(e.params_beta) <- NULL


# CODA summary ------------------------------------------------------------

post.smry <- summary(
  post,
  vars = c('deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov', 'resid.vcov') 
) |>  
  as.data.frame() |> 
  rownames_to_column('metric') |>
  select(metric, SSeff:psrf) |> 
  pivot_longer(-metric, names_to = 'diag', values_to = 'values') 

diag.smry <- post.smry |> 
  group_by(diag) |> 
  summarize(
    median = median(values),
    lower = unname(quantile(values, 0.025)),
    upper = unname(quantile(values, 0.975)),
    .groups = 'drop'
  )


# Posterior Predictive Check ----------------------------------------------

length.obs <- cbind(Trunk = df$trunk, Tail = df$tail)

ppc <- expand_grid(
  metric = colnames(length.obs),
  id = 1:nrow(length.obs)
) |> 
  mutate(metric = factor(metric, colnames(length.obs)))

ppc$pct.gte.obs <- sapply(1:nrow(ppc), function(i) {
  obs <- length.obs[ppc$id[i], ppc$metric[i]]
  ppd <- p$length.ppd[ppc$id[i], ppc$metric[i], ]
  mean(obs >= ppd)
})

ppc$mean.diff <- sapply(1:nrow(ppc), function(i) {
  obs <- length.obs[ppc$id[i], ppc$metric[i]]
  ppd <- p$length.ppd[ppc$id[i], ppc$metric[i], ]
  mean(obs - ppd)
})

ppc.smry <- ppc |> 
  group_by(metric) |> 
  summarize(
    median.pct = median(pct.gte.obs),
    lower.pct = unname(quantile(pct.gte.obs, 0.025)),
    upper.pct = unname(quantile(pct.gte.obs, 0.975)),   
    median.diff = median(mean.diff),
    lower.diff = unname(quantile(mean.diff, 0.025)),
    upper.diff = unname(quantile(mean.diff, 0.975)),
    .groups = 'drop'
  )


# Save all objects
save.image(format(end, "Model_outputs/Model_I_posterior_%Y%m%d_%H%M.rdata"))


# Plot posterior distributions
plot(
  post, 
  vars = c('deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov', 'resid.vcov'),
  file = format(end, "Model_outputs/Model_I_plots_%Y%m%d_%H%M.pdf")
)


# Plot diagnostics
pdf(format(end, "Model_outputs/Model_I_diagnostics_%Y%m%d_%H%M.pdf"))

ggplot(post.smry) +
  geom_histogram(aes(values), bins = 20) +
  facet_wrap(~diag, scales = 'free_x')

ggplot(ppc) +
  geom_histogram(aes(pct.gte.obs), binwidth = 0.05) +
  geom_vline(aes(xintercept = median.pct), data = ppc.smry, color = 'red') +
  geom_vline(aes(xintercept = lower.pct), data = ppc.smry, linetype = 'dashed', color = 'red') +
  geom_vline(aes(xintercept = upper.pct), data = ppc.smry, linetype = 'dashed', color = 'red') +
  facet_wrap(~ metric) +
  labs(x = 'Percent of PPD >= Observed', y = 'Count')

ggplot(ppc) +
  geom_histogram(aes(mean.diff), bins = 50) +
  geom_vline(aes(xintercept = median.diff), data = ppc.smry, color = 'red') +
  geom_vline(aes(xintercept = lower.diff), data = ppc.smry, linetype = 'dashed', color = 'red') +
  geom_vline(aes(xintercept = upper.diff), data = ppc.smry, linetype = 'dashed', color = 'red') +
  facet_wrap(~ metric, scales = 'free_x') +
  labs(x = 'Length Difference (Observed - PPD)', y = 'Count')

dev.off()


print(elapsed)