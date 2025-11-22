rm(list=ls())
library('tidyverse')
library('ggridges')
source('0_misc_funcs.R')

# options(scipen = 999)


# load and prepare ----
load("Model_outputs/Model_IV_posterior_20251001_0640.rdata") 


# Trunk selection differential  
cov_trunk_p <- vcv.obs$vcv.G.obs[1,3,]
mean_p <- vcv.obs$mean.obs[, "Settling"]
cov_trunk_p_relative <- cov_trunk_p / mean_p # selection differential

# Tail selection differential
cov_tail_p <- vcv.obs$vcv.G.obs[2,3,]
mean_p <- vcv.obs$mean.obs[, "Settling"]
cov_tail_p_relative <- cov_tail_p / mean_p # selection differential

# Extract means
mean_z <- rbind(p$mean.overall[1,], p$mean.overall[2,])

# Compute trunk and tail evolvability
sg <- rbind(cov_trunk_p_relative, cov_tail_p_relative)  # 2 x n_iter

n_iter <- dim(sg)[2]

E_matrix <- matrix(NA, nrow = 2, ncol = n_iter,
                   dimnames = list(c("Trunk", "Tail"), NULL))

for (i in seq_len(n_iter)) {
  S_vec <- sg[, i]
  mean_vec <- mean_z[, i]
  E_matrix[, i] <- S_vec / mean_vec
}


# Make Figure ----

## Plotting parameters ----
brksA<- seq(-0.5,1,0.1)
brksB<- seq(-0.5,1,0.1)
lmtsA <- c(-0.5,1)
lmtsB <- c(-0.5,1)
x_text_size <- 5
y_text_size <- 7
axis_label_size <- 7
title_label_size <- 7
point_size <- 2
segment_size <- 0.5
alp <- 0.4


## Panel A ----
d <- E_matrix['Trunk', , drop = F] * 100

df <- data.frame(
  sample = rownames(d),       
  value = as.numeric(d))

summaries <- as.data.frame(t(vecSmry(df$value)))
summaries$y <- 0

panelA <- ggplot(df, 
                 aes(x = value)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
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
d <- E_matrix['Tail', , drop = F] * 100

df <- data.frame(
  sample = rownames(d),       
  value = as.numeric(d))

summaries <- as.data.frame(t(vecSmry(df$value)))
summaries$y <- 0

panelB <- ggplot(df, 
                 aes(x = value)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
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