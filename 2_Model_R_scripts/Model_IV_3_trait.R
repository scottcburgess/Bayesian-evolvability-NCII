rm(list = ls())
library(tidyverse)
library(runjags)

# MCMC parameters
chains <- 5 #10
adapt <- 100
burnin <- 1000 #50000
total.sample <- 10000 #50000 
thin <- 10 #100

# Load data
head_tail.df <- readRDS("../1_Data/head_tail_data.rds") 
hatch_settle.df <- readRDS("../1_Data/hatch_settle_data.rds") 

blocks.to.keep <- table(
  block = hatch_settle.df$block,
  metric = hatch_settle.df$metric
) |> 
  as.data.frame() |> 
  filter(Freq > 0) |> 
  group_by(block) |> 
  summarize(n = n(), .groups = "drop") |> 
  filter(n == 2) |> 
  pull(block) |> 
  as.character() |> 
  as.integer()

hatch_settle.df <- hatch_settle.df |> 
  filter(block %in% blocks.to.keep & metric == "settling") |> 
  select(-metric)

interactions <- intersect(head_tail.df$interaction, hatch_settle.df$interaction)

head_tail.df <- filter(head_tail.df, interaction %in% interactions)
hatch_settle.df <- filter(hatch_settle.df, interaction %in% interactions)

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
    n.blocks = length(unique(head_tail.df$block)),
    n.sires = length(unique(head_tail.df$sire)),
    n.dams = length(unique(head_tail.df$dam)),
    n.interactions = length(unique(head_tail.df$interaction)),
    n.larvae = nrow(head_tail.df),
    block = as.numeric(factor(head_tail.df$block)),
    sire = as.numeric(factor(head_tail.df$sire)),
    dam = as.numeric(factor(head_tail.df$dam)),
    interaction = as.numeric(factor(head_tail.df$interaction)),
    block.mean.range = cbind(
      round(range(head_tail.df$head)), 
      round(range(head_tail.df$tail)),
      qlogis(c(0.2, 0.8)) #qlogis(c(min(block.settle$lower), max(block.settle$upper)))
    ),
    trait = cbind(head_tail.df$head, head_tail.df$tail),
    beta.var = 1,
    n.settle = nrow(hatch_settle.df),
    settle.block = as.numeric(factor(hatch_settle.df$block)),
    settle.sire = as.numeric(factor(hatch_settle.df$sire)),
    settle.dam = as.numeric(factor(hatch_settle.df$dam)),
    settle.interaction = as.numeric(factor(hatch_settle.df$interaction)),
    settle = hatch_settle.df$outcome,
    i = c(1, 1, 2),
    j = c(2, 3, 3)
  ),
  model = "model {
    # for each t-trait...
    for(t in 1:3) {  
      # ---- priors for block and effect means ----
      for(b in 1:n.blocks) {
        block.mean[b, t] ~ dunif(block.mean.range[1, t], block.mean.range[2, t])
      }
      additive.mean[t] ~ dnorm(0, 1e-5)
      dam.mean[t] ~ dnorm(0, 1e-5)
      interaction.mean[t] ~ dnorm(0, 1e-5)
      
      # ---- prior for variances ----
      additive.vcov[t, t] ~ dunif(0, 1e3)
      maternal.vcov[t, t] ~ dunif(0, 1e3)
      interaction.vcov[t, t] ~ dunif(0, 1e3)
      resid.vcov[t, t] ~ dunif(0, 1e3)
    }
    
    
    # ---- priors for settling linear model coefficients ----
    intercept ~ dnorm(0, 1 / 100)
    block.beta ~ dnorm(0, 1 / beta.var)
    sire.beta ~ dnorm(0, 1 / beta.var)
    maternal.beta ~ dnorm(0, 1 / beta.var)
    int.beta ~ dnorm(0, 1 / beta.var)
  
    
    # ---- correlation priors for covariances ----
    additive.corr ~ dunif(-1, 1)
    maternal.corr ~ dunif(-1, 1)
    interaction.corr ~ dunif(-1, 1)
    resid.corr ~ dunif(-1, 1)
    
    
    # ---- construct variance/covariance matrices ----
    for(k in 1:3) {
      additive.vcov[i[k], j[k]] <- additive.corr * sqrt(additive.vcov[i[k], i[k]] * additive.vcov[j[k], j[k]])
      additive.vcov[j[k], i[k]] <- additive.vcov[i[k], j[k]]
      maternal.vcov[i[k], j[k]]  <- maternal.corr * sqrt(maternal.vcov[i[k], i[k]] * maternal.vcov[j[k], j[k]])
      maternal.vcov[j[k], i[k]] <- maternal.vcov[i[k], j[k]]
      interaction.vcov[i[k], j[k]] <- interaction.corr * sqrt(interaction.vcov[i[k], i[k]] * interaction.vcov[j[k], j[k]])
      interaction.vcov[j[k], i[k]] <- interaction.vcov[i[k], j[k]]
      resid.vcov[i[k], j[k]] <- resid.corr * sqrt(resid.vcov[i[k], i[k]] * resid.vcov[j[k], j[k]])
      resid.vcov[j[k], i[k]] <- resid.vcov[i[k], j[k]]
    }
    
    # ---- prior for additive sire effect (for each sire) ----
    for(s in 1:n.sires) {
      additive.eff[s, 1:3] ~ dmnorm.vcov(additive.mean, additive.vcov)
    }
    
    # ---- prior for maternal effect (for each dam) ----
    for(d in 1:n.dams) {
      maternal.eff[d, 1:3] ~ dmnorm.vcov(dam.mean - additive.mean, maternal.vcov)
    }
    
    # ---- prior for interaction effect (for each sire x dam interaction) ----
    for(int in 1:n.interactions) {
      interaction.eff[int, 1:3] ~ dmnorm.vcov(interaction.mean, interaction.vcov)
    }
    
    
    # ---- head/tail likelihood ----
    for(l in 1:n.larvae) {
      for(t in 1:2) {
        # expected mean for the l-th larvae and t-th trait
        mu[l, t] <- block.mean[block[l], t] + 
          (2 * additive.eff[sire[l], t]) + 
          maternal.eff[dam[l], t] +
          interaction.eff[interaction[l], t]
      }
      # likelihood of l-th larvae for both traits from multivariate normal
      trait[l, ] ~ dmnorm.vcov(mu[l, ], resid.vcov[1:2, 1:2])
    }
    
    
    # ---- likelihood of settling ----
    for(s in 1:n.settle) {
      logit(pr.settle[s]) <- intercept + 
        (block.beta * block.mean[settle.block[s], 3]) +
        (sire.beta * additive.eff[settle.sire[s], 3]) +
        (maternal.beta * maternal.eff[settle.dam[s], 3]) +
        (int.beta * interaction.eff[settle.interaction[s], 3])
      settle[s] ~ dbern(pr.settle[s])
    }
  }",
  monitor = c(
    "deviance", "additive.vcov", "maternal.vcov", "interaction.vcov", 
    "resid.vcov", 'intercept', 'block.eff', 'additive.eff', 'maternal.eff', 
    'interaction.eff', 'block.beta', 'sire.beta', 'maternal.beta',
    'int.beta', 'pr.settle'
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
  method = "parallel",
  summarise = FALSE
)
end <- Sys.time()
elapsed <- swfscMisc::autoUnits(post$timetaken)


# 
# # Extract posterior to list of arrays - p
# p <- swfscMisc::runjags2list(post)
# rownames(p$mean.overall) <- c("head", "tail")
# dimnames(p$additive.vcov)[1:2] <- 
#   dimnames(p$maternal.vcov)[1:2] <- 
#   dimnames(p$interaction.vcov)[1:2] <- 
#   dimnames(p$resid.vcov)[1:2] <- 
#   list(rownames(p$mean.overall), rownames(p$mean.overall))
# 
# # Add QG metrics to list
# p$VA <- 4 * p$additive.vcov
# p$VM <- p$maternal.vcov
# p$VD <- 4 * p$interaction.vcov
# p$VP <- 2 * p$additive.vcov + p$VM + p$interaction.vcov + p$resid.vcov
# p$H <- p$VA / p$VP
# p$E <- rbind(
#   head = p$VA[1, 1, ] / (p$mean.overall[1, ] ^ 2),
#   tail = p$VA[2, 2, ] / (p$mean.overall[2, ] ^ 2)
# )
# 
# # Compute evolvability
# e.params_means <- do.call(
#   rbind,
#   parallel::mclapply(1:dim(p$VA)[3], function(i) {
#     evolvability::evolvabilityMeans(
#       G = as.vector(p$VA[, , i]),
#       means = p$mean.overall[, i]
#     )
#   }, mc.cores = 14) 
# )
#   
# e.params_BetaMCMC <- evolvability::evolvabilityBetaMCMC(
#   G_mcmc = evolvability::meanStdGMCMC(
#     t(apply(p$VA, 3, as.vector)),
#     t(p$mean.overall)
#   ),
#   Beta = evolvability::randomBeta(1000, 2),
#   post.dist = TRUE
# )
# 
# # Save all objects and plot posterior summaries
# save.image(format(end, "../3_Model_outputs/Model_I_posterior_%Y%m%d_%H%M.rdata"))
# 
# plot(
#   post, 
#   file = format(end, "../3_Model_outputs/Model_I_plots_%Y%m%d_%H%M.pdf")
# )
# 
# print(elapsed)
