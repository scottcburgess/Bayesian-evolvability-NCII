rm(list = ls())
library(tidyverse)
library(runjags)
library(parallel)
source("0_misc_funcs.R")

chains <- 14
adapt <- 100
burnin <- 500000
total.sample <- 1400
thin <- 100

df <- readRDS("hatch_settle_data.rds") 

blocks.to.keep <- table(block = df$block, metric = df$metric) %>% 
  as.data.frame() %>% 
  filter(Freq > 0) %>% 
  group_by(block) %>% 
  summarize(n = n(), .groups = "drop") %>% 
  filter(n == 2) %>% 
  pull(block) %>% 
  as.character() %>% 
  as.integer()

settle.df <- df %>% 
  filter(block %in% blocks.to.keep & metric == "settling") %>% 
  select(-metric)

length.df <- readRDS("head_tail_data.rds") 

interactions <- intersect(settle.df$interaction, length.df$interaction)

settle.df <- filter(settle.df, interaction %in% interactions)
length.df <- filter(length.df, interaction %in% interactions)

int.df <- length.df %>% 
  group_by(interaction) %>% 
  summarize(
    block = unique(block),
    sire = unique(sire),
    dam = unique(dam),
    .groups = "drop"
  ) %>% 
  mutate(
    block.num = as.numeric(factor(block)),
    sire.num = as.numeric(factor(sire)),
    dam.num = as.numeric(factor(dam)),
    interaction.num = as.numeric(factor(interaction))
  ) %>% 
  select(-block, -sire, -dam)

settle.df <- left_join(settle.df, int.df, by = "interaction")
length.df <- left_join(length.df, int.df, by = "interaction")

model.data <- list(
  beta.var = 1,
  max.var = 10, 
  n.blocks = max(int.df$block.num),
  n.sires = max(int.df$sire.num),
  n.dams = max(int.df$dam.num),
  n.interactions = nrow(int.df),
  block = int.df$block.num,
  sire = int.df$sire.num,
  dam = int.df$dam.num,
  n.length = nrow(length.df),
  length.range = cbind(range(length.df$head), range(length.df$tail)),
  length = as.matrix(length.df[, c("head", "tail")]),
  length.interaction = length.df$interaction.num,
  n.settle = nrow(settle.df),
  settle = cbind(settle.df$outcome, settle.df$outcome),
  settle.interaction = settle.df$interaction.num
)

post <- run.jags(
  data = model.data,
  model = "model {
    for(t in 1:2) {
      # ---- variance ----
      additive.var[t] ~ dunif(0, max.var)
      maternal.var[t] ~ dunif(0, max.var)
      block.var[t] ~ dunif(0, max.var)
      
      # ---- effects ----
      for(s in 1:n.sires) {
        additive.sire.eff[s, t] ~ dnorm(0, 1 / additive.var[t])
      }
  
      for(d in 1:n.dams) {
        additive.dam.eff[d, t] ~ dnorm(0, 1 / additive.var[t])
        maternal.eff[d, t] ~ dnorm(0, 1 / maternal.var[t])
      }
    
      for(b in 1:n.blocks) {
        block.eff[b, t] ~ dnorm(0, 1 / block.var[t])
      }

      # ---- logistic linear equation priors ----
      intercept[t] ~ dnorm(0, 1 / 100)
      int.beta[t] ~ dnorm(0, 1 / beta.var)
      sire.beta[t] ~ dnorm(0, 1 / beta.var)
      maternal.beta[t] ~ dnorm(0, 1 / beta.var)
      block.beta[t] ~ dnorm(0, 1 / beta.var)
  
      for(i in 1:n.interactions) {  
        interaction.mean[i, t] ~ dunif(length.range[1, t], length.range[2, t])
        interaction.var[i, t] ~ dunif(0, max.var)
      
        # expected mean length in interaction
        mean.length[i, t] <- interaction.mean[i, t] + 
          additive.sire.eff[sire[i], t] + 
          additive.dam.eff[dam[i], t] +
          maternal.eff[dam[i], t] +
          block.eff[block[i], t]
      
        # linear equation for expected probability of settling
        logit(pr.settle[i, t]) <- intercept[t] + 
          (int.beta[t] * interaction.mean[i, t]) +
          (sire.beta[t] * additive.sire.eff[sire[i], t]) +
          (maternal.beta[t] * maternal.eff[dam[i], t]) +
          (block.beta[t] * block.eff[block[i], t])
      }
      
      # likelihood of head or tail length
      for(l in 1:n.length) {
        length[l, t] ~ dnorm(
          mean.length[length.interaction[l], t], 
          1 / interaction.var[length.interaction[l], t]
        )
      }
    
      # likelihood of settling
      for(l in 1:n.settle) {
        settle[l, t] ~ dbern(pr.settle[settle.interaction[l], t])
      }
    }
  }",
  monitor = c(
    "deviance", "intercept", "pr.settle",
    "int.beta", "sire.beta", "maternal.beta", "block.beta",
    "interaction.mean", "additive.sire.eff", "maternal.eff", "block.eff"
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

p <- swfscMisc::runjags2list(post)
dimnames(p$intercept)[[1]] <- 
  dimnames(p$int.beta)[[1]] <-
  dimnames(p$sire.beta)[[1]] <-
  dimnames(p$maternal.beta)[[1]] <-
  dimnames(p$block.beta)[[1]] <-
  dimnames(p$pr.settle)[[2]] <-
  dimnames(p$interaction.mean)[[2]] <-
  dimnames(p$additive.sire.eff)[[2]] <-
  dimnames(p$maternal.eff)[[2]] <-
  dimnames(p$block.eff)[[2]] <-
  c("head", "tail")

mean.x.vec <- function(mat, n = 100) {
  apply(mat, 2, function(x) seq(min(x), max(x), length.out = n))
}

int.mean.x <- mean.x.vec(p$interaction.mean)
sire.eff.x <- mean.x.vec(p$additive.sire.eff)
maternal.eff.x <- mean.x.vec(p$maternal.eff)
block.eff.x <- mean.x.vec(p$block.eff)

pred.pr.settle <- lapply(colnames(sire.eff.x), function(m) {
  int <- mclapply(int.mean.x[, m], function(x) {
    sapply(1:model.data$n.interactions, function(i) {
      p$intercept[m, ] + 
        (p$int.beta[m, ] * x) +
        (p$sire.beta[m, ] * t(p$additive.sire.eff[model.data$sire[i], m, ])) +
        (p$maternal.beta[m, ] * t(p$maternal.eff[model.data$dam[i], m, ])) +
        (p$block.beta[m, ] * t(p$block.eff[model.data$block[i], m, ]))
    }) %>% 
      as.vector() %>% 
      vecSmry() %>% 
      plogis() %>% 
      c(effect = x)
  }, mc.cores = 14) %>% 
    do.call(rbind, .) %>% 
    as.data.frame() %>% 
    mutate(metric = m, effect.type = "Interaction Mean")
  
  sire <- mclapply(sire.eff.x[, m], function(x) {
    sapply(1:model.data$n.interactions, function(i) {
      p$intercept[m, ] + 
        (p$int.beta[m, ] * t(p$interaction.mean[i, m, ])) +
        (p$sire.beta[m, ] * x) +
        (p$maternal.beta[m, ] * t(p$maternal.eff[model.data$dam[i], m, ])) +
        (p$block.beta[m, ] * t(p$block.eff[model.data$block[i], m, ]))
    }) %>% 
      as.vector() %>% 
      vecSmry() %>% 
      plogis() %>% 
      c(effect = x)
  }, mc.cores = 14) %>% 
    do.call(rbind, .) %>% 
    as.data.frame() %>% 
    mutate(metric = m, effect.type = "Additive Sire")
  
  maternal <- mclapply(maternal.eff.x[, m], function(x) {
    sapply(1:model.data$n.interactions, function(i) {
      p$intercept[m, ] + 
        (p$int.beta[m, ] * t(p$interaction.mean[i, m, ])) +
        (p$sire.beta[m, ] * t(p$additive.sire.eff[model.data$sire[i], m, ])) +
        p$maternal.beta[m, ] * x +
        (p$block.beta[m, ] * t(p$block.eff[model.data$block[i], m, ]))
    }) %>% 
      as.vector() %>% 
      vecSmry() %>% 
      plogis() %>% 
      c(effect = x)
  }, mc.cores = 14) %>% 
    do.call(rbind, .) %>% 
    as.data.frame() %>% 
    mutate(metric = m, effect.type = "Maternal")
  
  block <- mclapply(block.eff.x[, m], function(x) {
    sapply(1:model.data$n.interactions, function(i) {
      p$intercept[m, ] + 
        (p$int.beta[m, ] * t(p$interaction.mean[i, m, ])) +
        (p$sire.beta[m, ] * t(p$additive.sire.eff[model.data$sire[i], m, ])) +
        (p$maternal.beta[m, ] * t(p$maternal.eff[model.data$dam[i], m, ])) +
        (p$block.beta[m, ] * x)
    }) %>% 
      as.vector() %>% 
      vecSmry() %>% 
      plogis() %>% 
      c(effect = x)
  }, mc.cores = 14) %>% 
    do.call(rbind, .) %>% 
    as.data.frame() %>% 
    mutate(metric = m, effect.type = "Block")
  
  bind_rows(int, sire, maternal, block)
}) %>% 
  bind_rows()

save.image(format(end, "Model_IV_posterior_%Y%m%d_%H%M.rdata"))

plot(
  post, 
  vars = c("deviance", "intercept", "int.beta", "sire.beta", "maternal.beta", "block.beta"),
  file = format(end, "Model_IV_plots_%Y%m%d_%H%M.pdf")
)

print(elapsed)
