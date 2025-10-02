rm(list=ls())
library('tidyverse')
library('ggridges')
source('0_misc_funcs.R')

options(scipen = 999)


# trunk Tail: load and prepare ----
Data1 <- new.env()
load("Model_outputs/Model_I_posterior_20250930_0019.rdata", envir = Data1) 

# Prepare data 
# Get the posterior samples by averaging across all of the random beta's
eB_posterior_trunk_tail <- apply(Data1$e.params_BetaMCMC$post.dist$eB, 1, mean)
rB_posterior_trunk_tail <- apply(Data1$e.params_BetaMCMC$post.dist$rB, 1, mean)
cB_posterior_trunk_tail <- apply(Data1$e.params_BetaMCMC$post.dist$cB, 1, mean)
# Check
# e.params_BetaMCMC$summary # median (called e_mean) should be the same as
# vecSmry(eB_posterior_trunk_tail) # median here (but we're using the mode)

E_trunk_tail <- 
  data.frame(name = c(rep("eB",length(eB_posterior_trunk_tail)),
                      rep("rB",length(rB_posterior_trunk_tail)),
                      rep("cB",length(cB_posterior_trunk_tail))),
             value = c(eB_posterior_trunk_tail,
                       rB_posterior_trunk_tail,
                       cB_posterior_trunk_tail))

E_trunk_tail <- E_trunk_tail |> 
  mutate(name = factor(name, levels = c("eB","rB","cB")))


# Hatching Settling : load and prepare ----
Data2 <- new.env()
load("Model_outputs/Model_III_posterior_20250930_1554.rdata", envir = Data2) 

## Prepare data 
eB_posterior_hatch_settle <- apply(Data2$e.params_BetaMCMC$post.dist$eB, 1, mean)
rB_posterior_hatch_settle <- apply(Data2$e.params_BetaMCMC$post.dist$rB, 1, mean)
cB_posterior_hatch_settle <- apply(Data2$e.params_BetaMCMC$post.dist$cB, 1, mean)
# Check
# e.params_BetaMCMC$summary # median (called e_mean) should be the same as
# vecSmry(eB_posterior_hatch_settle) # median here (but we're using the mode)

E_hatch_settle <- data.frame(name = c(rep("eB",length(eB_posterior_hatch_settle)),
                                      rep("rB",length(rB_posterior_hatch_settle)),
                                      rep("cB",length(cB_posterior_hatch_settle))),
                             value = c(eB_posterior_hatch_settle,
                                       rB_posterior_hatch_settle,
                                       cB_posterior_hatch_settle))
    
E_hatch_settle <- E_hatch_settle |> 
  mutate(name = factor(name, levels = c("eB","rB","cB")))





# Make Figure ----
vc_color <- data.frame(color = c("#0077b6",
                                 "#00b4d8",
                                 "#73e8ff"))


## Plotting parameters ----
brksA <- seq(0,0.002,0.0001)
brksB <- seq(0,0.1,0.01)
lmtsA <- c(0,0.0008)
lmtsB <- c(0,0.06)
x_text_size <- 5
y_text_size <- 7
axis_label_size <- 9
title_label_size <- 9
point_size <- 2
segment_size <- 0.5
alp <- 0.4
bw1 <- 0.000015
bw2 <- 0.001
sc <- 0.9

## Panel A ----
summaries <- E_trunk_tail |>
  group_by(name) |>
  summarise(summary_values = list(vecSmry(value)), .groups = 'drop') |>
  unnest_wider(summary_values, names_repair = "unique")


panelA <- ggplot(E_trunk_tail, 
                 aes(x = value, y = name, fill = name)) +
  geom_density_ridges(scale = sc,
                      alpha = alp,
                      bandwidth = bw1,
                      color = 'lightgrey',
                      linewidth = 0.1) +
  theme_ridges() + 
  labs(x = "Evolvability",
       y = "Metric",
       title = "a) Trunk-Tail length") +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size, angle = 45),
        axis.text.y = element_text(size = y_text_size),
        axis.title = element_text(size = axis_label_size),
        plot.title = element_text(size=title_label_size, face = "plain")) +
  geom_point(data = summaries, 
             aes(x = mode, 
                 y = as.numeric(name),
                 color = name), 
             size = point_size) +
  geom_segment(data = summaries, 
               aes(x = lower.hdi, 
                   xend = upper.hdi, 
                   y = as.numeric(name), 
                   yend = as.numeric(name),
                   color = name), 
               linetype = "solid", 
               size = segment_size) +
  scale_x_continuous(breaks = brksA,
                     limits = lmtsA) +
  scale_fill_manual(values = c("eB" = vc_color[1,], 
                               "rB" = vc_color[2,],
                               "cB" = vc_color[3,])) +
  scale_y_discrete(labels = c("eB" = expression(paste("e(",beta,")")), 
                              "rB" = expression(paste("r(",beta,")")),
                              "cB" = expression(paste("c(",beta,")")))) + 
  scale_color_manual(values = c("eB" = vc_color[1,], 
                                "rB" = vc_color[2,],
                                "cB" = vc_color[3,]))

## Panel B ----
summaries <- E_hatch_settle |>
  group_by(name) |>
  summarise(summary_values = list(vecSmry(value)), .groups = 'drop') |>
  unnest_wider(summary_values, names_repair = "unique")


panelB <- ggplot(E_hatch_settle, 
                 aes(x = value, y = name, fill = name)) +
  geom_density_ridges(scale = sc,
                      alpha = alp,
                      bandwidth = bw2,
                      color = 'lightgrey',
                      linewidth = 0.1) +
  theme_ridges() + 
  labs(x = "Evolvability",
       y = "Metric",
       title = "b) Hatch-Settle") +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size, angle = 45),
        axis.text.y = element_text(size = y_text_size),
        axis.title = element_text(size = axis_label_size),
        plot.title = element_text(size=title_label_size, face = "plain")) +
  geom_point(data = summaries, 
             aes(x = mode, 
                 y = as.numeric(name),
                 color = name), 
             size = point_size) +
  geom_segment(data = summaries, 
               aes(x = lower.hdi, 
                   xend = upper.hdi, 
                   y = as.numeric(name), 
                   yend = as.numeric(name),
                   color = name), 
               linetype = "solid", 
               size = segment_size) +
  scale_x_continuous(breaks = brksB,
                     limits = lmtsB) +
  scale_fill_manual(values = c("eB" = vc_color[1,], 
                               "rB" = vc_color[2,],
                               "cB" = vc_color[3,])) +
  scale_y_discrete(labels = c("eB" = expression(paste("e(",beta,")")), 
                              "rB" = expression(paste("r(",beta,")")),
                              "cB" = expression(paste("c(",beta,")")))) + 
  scale_color_manual(values = c("eB" = vc_color[1,], 
                                "rB" = vc_color[2,],
                                "cB" = vc_color[3,]))



# Save plot ----
fig4 <- gridExtra::grid.arrange(panelA,
                                panelB,
                        nrow = 1,
                        ncol = 2)
ggsave("Figures and Tables/Figure 4.pdf", 
       plot = fig4, 
       height = 2, 
       width = 5)
dev.off()



# Summaries ----
tmp <- E_hatch_settle |> filter(name=="eB") |> select("value")
round(vecSmry(tmp$value),4) # trunk VA (proportion)

tmp <- E_trunk_tail |> filter(name=="eB") |> select("value")
round(vecSmry(tmp$value),6) # Tail VA (proportion)
