rm(list = ls())
library(tidyverse)
library(ellipse)
source("0_misc_funcs.R")

load("Model_outputs/Model_I_posterior_20250930_0019.rdata")

# Mean standardize, calculate mode, and get the ellipse
get_ellipse <- function(mat){
  foo <- matrix(c(
    vecSmry(mat[,1])['mode'],
    vecSmry(mat[,2])['mode'],
    vecSmry(mat[,2])['mode'],
    vecSmry(mat[,4])['mode']
  ), nrow = 2, byrow = TRUE)
  as.data.frame(ellipse(foo, level = 0.95))  
}

ell_VA <- get_ellipse(evolvability::meanStdGMCMC(
  t(apply(p$VA, 3, as.vector)),
  t(p$mean.overall)))
ell_VM <- get_ellipse(evolvability::meanStdGMCMC(
  t(apply(p$VM, 3, as.vector)),
  t(p$mean.overall)))
ell_VD <- get_ellipse(evolvability::meanStdGMCMC(
  t(apply(p$VD, 3, as.vector)),
  t(p$mean.overall)))
ell_VR <- get_ellipse(evolvability::meanStdGMCMC(
  t(apply(p$resid.vcov, 3, as.vector)),
  t(p$mean.overall)))
ell_VP <- get_ellipse(evolvability::meanStdGMCMC(
  t(apply(p$VP, 3, as.vector)),
  t(p$mean.overall)))

# Combine for plotting
ellipses <- bind_rows(
  ell_VD %>% mutate(type = "VD"),
  ell_VM %>% mutate(type = "VM"),
  ell_VA %>% mutate(type = "VA"),
  ell_VR %>% mutate(type = "VR"),
  ell_VP %>% mutate(type = "VP")
)

# Set preferred order
ellipses$type <- factor(ellipses$type, 
                        levels = c("VA", 
                                   "VM", 
                                   "VD", 
                                   "VR",
                                   "VP"))

# Custom colors 
Matrix <- c(
  VD = "#ffd166",
  VM = "#118ab2",
  VA = "#ef476f",
  VR = "#06d6a0",
  VP = "grey30"
)

# Custom labels
legend_labels <- c(
  VA = expression(paste("Additive genetic (G)", ~ Pr(G[12] > 0), " = 59%")),
  VM = expression(paste("Maternal effect (M)", ~ Pr(M[12] > 0), " = 95%")),
  VD = expression(paste("Dominance (D)", ~ Pr(D[12] > 0), " = 83%")),
  VR = expression(paste("Residual (R)", ~ Pr(R[12] > 0), " = 100%")),
  VP = expression(paste("Phenotypic (P)", ~ Pr(P[12] > 0), " = 100%"))
)


# Plot
p1 <- ggplot(ellipses, aes(x = x, y = y, color = type)) +
  geom_hline(yintercept = 0, size = 1, color = "gray") +
  geom_vline(xintercept = 0, size = 1, color = "gray") +
  geom_path(size = 1.3, alpha = 0.7) +
  scale_color_manual(values = Matrix, labels = legend_labels) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(x = "Trunk length", y = "Tail length", color = "")
p1

ggsave("Figures and Tables/Visualize_matrices.pdf", 
       plot = p1, 
       height = 4, 
       width = 8)
dev.off()
