rm(list=ls())
library('tidyverse')
library('ggridges')
source('0_misc_funcs.R')

options(scipen = 999)


# load and prepare ----
load("Model_outputs/Model_IV_posterior_20251001_0640.rdata") 


## G matrix  
G_matrix <- p$VA[1:2,1:2,]

# Trunk selection differential  
cov_trunk_p <- vcv.obs$vcv.G.obs[1,3,]
mean_p <- vcv.obs$mean.obs[, "Settling"]
cov_trunk_p_relative <- cov_trunk_p / mean_p # selection differential

# Tail selection differential
cov_tail_p <- vcv.obs$vcv.G.obs[2,3,]
mean_p <- vcv.obs$mean.obs[, "Settling"]
cov_tail_p_relative <- cov_tail_p / mean_p # selection differential

## Selection gradients
sg <- rbind(cov_trunk_p_relative, cov_tail_p_relative)  # 2 x n_iter

n_iter <- dim(G_matrix[,,])[3]

beta <- matrix(NA, nrow = 2, ncol = n_iter,
               dimnames = list(c("Trunk", "Tail"), NULL))

for (i in seq_len(n_iter)) {
  G <- G_matrix[, , i]
  S_vec <- sg[, i]
  beta[, i] <- solve(G) %*% S_vec
}

# The % change in trunk and tail length per generation (assuming fitness is settlement)   
# Extract means for convenience
mean_trunk <- p$mean.overall[1,]
mean_tail  <- p$mean.overall[2,]

# Extract variances and covariance
VA11 <- p$VA[1,1,]  # Var(Trunk)
VA22 <- p$VA[2,2,]  # Var(Tail)
VA12 <- p$VA[1,2,]  # Cov(Trunk, Tail)

# Compute trunk and tail evolvability
E_matrix <- rbind(
  Trunk = (VA11 / mean_trunk^2) * beta[1, ] + (VA12 / (mean_trunk * mean_tail)) * beta[2, ],
  Tail  = (VA22 / mean_tail^2)  * beta[2, ] + (VA12 / (mean_trunk * mean_tail)) * beta[1, ]
)


# Make Figure ----

## Plotting parameters ----
brksA<- seq(-0.003,0.01,0.001)
brksB<- seq(-0.003,0.01,0.0005)
lmtsA <- c(-0.002,0.003)
lmtsB <- c(-0.001,0.0035)
x_text_size <- 5
y_text_size <- 7
axis_label_size <- 7
title_label_size <- 7
point_size <- 2
segment_size <- 0.5
alp <- 0.4


## Panel A ----
d <- E_matrix['Trunk', ] * 100

df <- data.frame(
  sample = names(d),       
  value = as.numeric(d))

summaries <- as.data.frame(t(vecSmry(df$value)))
summaries$y <- 0

panelA <- ggplot(df, 
                 aes(x = value)) +
  geom_density(alpha = alp,
               fill = 'lightgrey',
               linewidth = 0.1) +
  labs(x = "% change per generation",
       y = "Probability density",
       title = "a) Trunk length") +
  theme_ridges() +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size, angle = 45),
        axis.text.y = element_text(size = y_text_size),
        axis.title = element_text(size = axis_label_size),
        plot.title = element_text(size=title_label_size, face = "plain", hjust = 0)) +
  geom_point(data = summaries, 
             aes(x = mode,
                 y = y),
             color = "grey30",
             size = point_size) +
  geom_segment(data = summaries, 
               aes(x = lower.hdi, 
                   xend = upper.hdi,
                   y = y,
                   yend = y),
               color = "grey30",
               linetype = "solid", 
               size = segment_size) +
  scale_x_continuous(breaks = brksA, limits = lmtsA) 



## Panel B ----
d <- E_matrix['Tail', ] * 100

df <- data.frame(
  sample = names(d),       
  value = as.numeric(d))

summaries <- as.data.frame(t(vecSmry(df$value)))
summaries$y <- 0

panelB <- ggplot(df, 
                 aes(x = value)) +
  geom_density(alpha = alp,
               fill = 'lightgrey',
               linewidth = 0.1) +
  labs(x = "% change per generation",
       y = "Probability density",
       title = "b) Tail length") +
  theme_ridges() +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size, angle = 45),
        axis.text.y = element_text(size = y_text_size),
        axis.title = element_text(size = axis_label_size),
        plot.title = element_text(size=title_label_size, face = "plain", hjust = 0)) +
  geom_point(data = summaries, 
             aes(x = mode,
                 y = y),
             color = "grey30",
             size = point_size) +
  geom_segment(data = summaries, 
               aes(x = lower.hdi, 
                   xend = upper.hdi,
                   y = y,
                   yend = y),
               color = "grey30",
               linetype = "solid", 
               size = segment_size) +
  scale_x_continuous(breaks = brksB, limits = lmtsB) 



fig6 <- gridExtra::grid.arrange(panelA,
                                panelB,
                                nrow = 1,
                                ncol = 2)
ggsave("Figures and Tables/Figure 6.pdf", 
       plot = fig6, 
       height = 2, 
       width = 5)
dev.off()