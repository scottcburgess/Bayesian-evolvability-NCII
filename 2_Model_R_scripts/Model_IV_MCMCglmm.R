# This is a re-imagined Model IV, based on reviewer comments.
# The model is constructed using MCMCglmm.
# The goal is to fit trunk length, tail length, and probability of settlement as a single G-matrix
# Trunk and tail length measured on the same larvae, 
# but settlement was measured on different larvae from the same full-sib family (=interaction),
# so cannot estimate residual covariances among all traits at the individual level
# (i.e, residual covariance between settle and head/tail is not estimable)
# So need to constrain the residual covariance matrix.
# Otherwise, settle can covary with head and tail at the sire,dam, and interaction level


rm(list = ls())
library(MCMCglmm)
library(tidyverse)
library(patchwork)
# library(parallel)
# library(coda)

# Load data
df_head_tail <- readRDS("1_Data/head_tail_data.rds") 
df_hatch_settle <- readRDS("1_Data/hatch_settle_data.rds") 

blocks_to_use <- df_hatch_settle %>% 
  filter(metric == "settling") %>% 
  pull(block) %>%
  unique()

df <- df_hatch_settle %>% 
  filter(metric == "settling") %>% 
  bind_rows(df_head_tail) %>% 
  select(animal, block, interaction, sire, dam, outcome, head, tail) %>% 
  rename(settle = outcome) %>% 
  filter(block %in% blocks_to_use)

df$block <- factor(df$block)
df$sire <- factor(df$sire) # each sire has a unique id
df$dam <- factor(df$dam) # each dam has a unique id
df$interaction <- factor(df$interaction)


# Set a weak prior
# Parameter expanded prior to help with mixing for small variances components.
# Proper prior since nu = number of traits
# In R, fix = 3 means makes the residual for the 3rd trait (settle)
# to 1 on the latent scale (probit) and forces its residual covariances with the other traits to 0

prior <- list(
  G = list(
    G1 = list(V = diag(3), nu = 3, alpha.mu = rep(0, 3), alpha.V = diag(3) * 1e06),  # sire
    G2 = list(V = diag(3), nu = 3, alpha.mu = rep(0, 3), alpha.V = diag(3) * 1e06),  # dam
    G3 = list(V = diag(3), nu = 3, alpha.mu = rep(0, 3), alpha.V = diag(3) * 1e06)   # interaction
  ),
  R = list(V = diag(3), nu = 2, fix = 3)  # residuals — no parameter expansion used here
)

nitt = 500000
burnin = 5000
thin = 50
(nitt - burnin) / thin # aiming for about 10,000

set.seed(999)
start.time <- Sys.time()
model <- MCMCglmm(
  fixed = cbind(head, tail, settle) ~ trait - 1 + trait:block,
  random = ~ us(trait):sire + 
    us(trait):dam + 
    us(trait):interaction,
  rcov = ~ us(trait):units,
  family = c("gaussian", "gaussian", "categorical"),
  data = df,
  prior = prior,
  nitt = nitt, burnin = burnin, thin = thin
)
end.time <- Sys.time()
end.time - start.time # ~ 1.5hrs

# saveRDS(model, "3_Model_outputs/Model_IV_MCMCglmm.rdata")

# Check diagnostics
effectiveSize(model$VCV) # Goal: >1000 good, >10000 preferred
autocorr.diag(model$VCV, lags = 1) # <0.1 for the first Lag is reasonable
heidel.diag(model$VCV)

# Check density and trace visually
# Extract G-structure columns
G_cols <- grep("sire|dam|interaction", colnames(model$VCV), value = TRUE)
G_df <- model$VCV[, G_cols]

# Convert to long format for plotting
G_long <- as.data.frame(G_df) %>%
  mutate(iter = 1:nrow(.)) %>%
  pivot_longer(-iter, names_to = "parameter", values_to = "value")

# Create PDF
pdf("3_Model_outputs/Model_IV_MCMCglmm_density_trace.pdf", width = 7, height = 8)

for (param in unique(G_long$parameter)) {
  
  df_param <- G_long %>% filter(parameter == param)
  
  # Density plot
  mode_est <- modeest::mlv(df_param$value, method = "venter")
  median_est <- median(df_param$value)
  p_density <- ggplot(df_param, aes(x = value)) +
    geom_density(fill = "darkblue", alpha = 0.5) +
    geom_vline(xintercept = mode_est, col = "darkblue", linetype = 'dashed') +
    geom_vline(xintercept = median_est, col = "grey", linetype = 'dotted') +
    labs(title = paste0(param, " posterior density"),
         x = "Value", y = "Density") +
    theme_minimal()
  
  # Trace plot
  p_trace <- ggplot(df_param, aes(x = iter, y = value)) +
    geom_line(color = "black") +
    labs(title = paste0(param, " trace"),
         x = "Iteration", y = "Value") +
    theme_minimal()
  
  # Combine vertically
  print(p_density / p_trace)
}

dev.off()
