# This is Model I converted to MCMCglmm, for comparison

rm(list = ls())
library(MCMCglmm)
library(parallel)
library(coda)

# Load data
df <- readRDS("1_Data/head_tail_data.rds") 

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
  fixed = cbind(head, tail) ~ trait - 1 + trait:block,
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
end.time - start.time # ~ 5 mins

# Check diagnostics
effectiveSize(model$VCV) # >1000 good, >10000 preferred
autocorr.diag(model$VCV, lags = 1) # <0.1 for the first Lag is reasonable
heidel.diag(model$VCV)

# Check chains
plot(model$VCV[,"traithead:traithead.sire"])
plot(model$VCV[,"traittail:traittail.sire"])
plot(model$VCV[,"traithead:traittail.sire"])
plot(model$VCV[,"traithead:traithead.dam"])
plot(model$VCV[,"traittail:traittail.dam"])
plot(model$VCV[,"traithead:traittail.dam"])
plot(model$VCV[,"traithead:traithead.interaction"])
plot(model$VCV[,"traittail:traittail.interaction"])
plot(model$VCV[,"traithead:traittail.interaction"])

# Function to extract posterior mode and HDI for VA, VM, VD
extract_matrix <- function(vcvs, labels, scale = 1, subtract = NULL) {
  get_component <- function(label) {
    x <- model$VCV[, label]
    if (!is.null(subtract)) {
      x <- x - model$VCV[, gsub(names(subtract), subtract, label)]
    }
    x <- x * scale
    c(
      # posterior.mode(x),
      modeest::mlv(x, method = "venter"),
      HDInterval::hdi(x),
      median(x)
    )
  }
  
  mat <- sapply(labels, get_component)
  mat <- as.data.frame(mat)
  rownames(mat) <- c("mode", "lwr", "upr", "median")
  return(mat)
}

# Labels
labels <- c(
  Head = "traithead:traithead.%s",
  Tail = "traittail:traittail.%s",
  Covariance = "traithead:traittail.%s"
)

# Format labels
format_labels <- function(source) {
  setNames(lapply(labels, sprintf, source), names(labels))
}

# Matrices
VAmat <- round(extract_matrix(model$VCV, format_labels("sire"), scale = 4),2)
VMmat <- round(extract_matrix(model$VCV, format_labels("dam"), subtract = c("dam" = "sire")),2)
VDmat <- round(extract_matrix(model$VCV, format_labels("interaction"), scale = 4),2)

VAmat; VMmat; VDmat


# Get VP
VCV <- model$VCV

# Extract variance components by source
vcv_head <- rowSums(VCV[, grep("traithead:traithead", colnames(VCV))])
vcv_tail <- rowSums(VCV[, grep("traittail:traittail", colnames(VCV))])
vcv_covariance <- rowSums(VCV[, grep("traithead:traittail", colnames(VCV))])

#VP trunk
round(modeest::mlv(vcv_head, method = "venter"),2)
round(median(vcv_head), 2)
round(HDInterval::hdi(vcv_head),2)

#VP tail
round(modeest::mlv(vcv_tail, method = "venter"),2)
round(median(vcv_tail),2)
round(HDInterval::hdi(vcv_tail),2)

#VP head
round(modeest::mlv(vcv_covariance, method = "venter"),2)
round(median(vcv_covariance),2)
round(HDInterval::hdi(vcv_covariance),2)
