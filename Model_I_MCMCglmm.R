# This is Model I converted to MCMCglmm, for comparison
rm(list = ls())
library(MCMCglmm)

# Load data
df <- readRDS("Data/trunk_tail_data.rds") 

df$block <- factor(df$block)
df$sire <- factor(df$sire) # each sire has a unique id
df$dam <- factor(df$dam) # each dam has a unique id
df$interaction <- factor(df$interaction)

# Set a weak prior
# Parameter expanded prior to help with mixing for small variances components.
# Proper prior since nu = number of traits
prior <- list(
  G = list(
    G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1e06),  # sire
    G2 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1e06),  # dam
    G3 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1e06)   # interaction
  ),
  R = list(V = diag(2), nu = 2)  # residuals — no parameter expansion used here
)

nitt = 500000
burnin = 5000
thin = 50
(nitt - burnin) / thin # aiming for about 10,000

set.seed(999)
start.time <- Sys.time()
model <- MCMCglmm(
  fixed = cbind(trunk, tail) ~ trait - 1 + trait:block,
  random = ~ us(trait):sire + 
            us(trait):dam + 
            us(trait):interaction,
  rcov = ~ us(trait):units,
  family = c("gaussian", "gaussian"),
  data = df,
  prior = prior,
  nitt = nitt, burnin = burnin, thin = thin
)
end.time <- Sys.time()
end.time - start.time # ~ 7 mins

# Check diagnostics
effectiveSize(model$VCV) # >1000 good, >10000 preferred
autocorr.diag(model$VCV, lags = 1) # <0.1 for the first Lag is reasonable
heidel.diag(model$VCV)

# load these after running MCMCglmm, because they seem to conflict
library(tidyverse)
library(patchwork)

# Check density and trace visually
# Extract G-structure columns
V_cols <- grep("sire|dam|interaction|units", colnames(model$VCV), value = TRUE)
V_df <- model$VCV[, V_cols]

# Convert to long format for plotting
V_long <- as.data.frame(V_df) %>% 
  mutate(iter = 1:nrow(.)) %>%
  pivot_longer(-iter, names_to = "parameter", values_to = "value")

# Create PDF
pdf("Model_outputs/Model_I_MCMCglmm_density_trace.pdf", width = 7, height = 8)

for (param in unique(V_long$parameter)) {
  
  df_param <- V_long %>% filter(parameter == param)
  
  # Density plot
  mode_est <- modeest::mlv(df_param$value, method = "venter")
  median_est <- median(df_param$value)
  p_density <- ggplot(df_param, aes(x = value)) +
    geom_histogram(bins = 200, fill = "darkblue", alpha = 0.5) +
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


# Function to extract posterior mode and HDI for VA, VM, VD
extract_matrix <- function(vcvs, labels, scale = 1, subtract = NULL) {
  get_component <- function(label) {
    x <- model$VCV[, label]
    if (!is.null(subtract)) {
      x <- x - model$VCV[, gsub(names(subtract), subtract, label)]
    }
    x <- x * scale
  }
  mat <- sapply(labels, get_component)
  mat <- as.data.frame(mat)
  return(mat)
}

summarize_matrix <- function(mat){
  foo <- apply(mat, 2, function(x) {
    mode = round(modeest::mlv(x, method = "venter"),2)
    median = round(median(x),2)
    HDI = round(HDInterval::hdi(x),2)
    c(mode = mode, median = median, lwd = HDI[1], upr = HDI[2])})
  return(foo)}
  
# Labels
labels <- c(
  Trunk = "traittrunk:traittrunk.%s",
  Tail = "traittail:traittail.%s",
  Covariance = "traittrunk:traittail.%s"
)

# Format labels
format_labels <- function(source) {
  setNames(lapply(labels, sprintf, source), names(labels))
}

# Matrices
VAmat <- extract_matrix(model$VCV, format_labels("sire"), scale = 4)
VMmat <- extract_matrix(model$VCV, format_labels("dam"), subtract = c("dam" = "sire"))
VDmat <- extract_matrix(model$VCV, format_labels("interaction"), scale = 4)
VRmat <- extract_matrix(model$VCV, format_labels("units"))
VPmat <- VAmat + VMmat + VDmat + VRmat
 
# Reported in Table S1
summarize_matrix(VAmat)
summarize_matrix(VMmat)
summarize_matrix(VDmat)
summarize_matrix(VRmat)
summarize_matrix(VPmat)
