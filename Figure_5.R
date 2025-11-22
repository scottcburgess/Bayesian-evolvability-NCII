rm(list=ls())
library('tidyverse')
library('ggridges')
source('0_misc_funcs.R')

options(scipen = 999)


# trunk Tail: load and prepare ----
load("Model_outputs/Model_I_posterior_20250930_0019.rdata") 

# Make Figure ----
vc_color <- data.frame(color = c("#0077b6",
                                 "#00b4d8",
                                 "#73e8ff"))


## Plotting parameters ----
brksA <- seq(0,0.01,0.0005)
brksB <- seq(0,0.01,0.0005)
brksC <- seq(0,0.01,0.0005)
lmtsA <- c(0,0.0025)
lmtsB <- c(0,0.0025)
lmtsC <- c(0,0.0025)
x_text_size <- 5
y_text_size <- 7
axis_label_size <- 7
title_label_size <- 6
point_size <- 2
segment_size <- 0.5
alp <- 0.4
bw1 <- 0.00005
bw2 <- 0.0004
bw3 <- 0.0004
sc <- 0.9

## Panel A ----
d <- e.params_beta |>
  filter(Beta_index == "1") |>
  select("e", "r", "c") |>
  pivot_longer(
    cols = c(e, r, c),
    names_to = "name",
    values_to = "value") |>
  mutate(name = factor(name, levels = c("e","r","c")))

summaries <- d |>
  group_by(name) |>
  summarise(summary_values = list(vecSmry(value)), .groups = 'drop') |>
  unnest_wider(summary_values, names_repair = "unique")


panelA <- ggplot(d, 
                 aes(x = value, y = name, fill = name)) +
  geom_density_ridges(scale = sc,
                      alpha = alp,
                      bandwidth = bw1,
                      color = 'lightgrey',
                      linewidth = 0.1) +
  theme_ridges() + 
  labs(x = "Evolvability (%)",
       y = "Metric",
       title = "a) Selection for longer tails") +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size, angle = 45),
        axis.text.y = element_text(size = y_text_size),
        axis.title = element_text(size = axis_label_size),
        plot.title.position = "plot",
        plot.title = element_text(size=title_label_size, face = "plain", hjust = 0)) +
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
                     limits = lmtsA,
                     labels = ~.x*100) +
  scale_fill_manual(values = c("e" = vc_color[1,], 
                               "r" = vc_color[2,],
                               "c" = vc_color[3,])) +
  scale_y_discrete(labels = c("e" = expression(paste("e(",beta,")")), 
                              "r" = expression(paste("r(",beta,")")),
                              "c" = expression(paste("c(",beta,")")))) + 
  scale_color_manual(values = c("e" = vc_color[1,], 
                                "r" = vc_color[2,],
                                "c" = vc_color[3,]))

## Panel B ----
d <- e.params_beta |>
  filter(Beta_index == "2") |>
  select("e", "r", "c") |>
  pivot_longer(
    cols = c(e, r, c),
    names_to = "name",
    values_to = "value") |>
  mutate(name = factor(name, levels = c("e","r","c")))

summaries <- d |>
  group_by(name) |>
  summarise(summary_values = list(vecSmry(value)), .groups = 'drop') |>
  unnest_wider(summary_values, names_repair = "unique")


panelB <- ggplot(d, 
                 aes(x = value, y = name, fill = name)) +
  geom_density_ridges(scale = sc,
                      alpha = alp,
                      bandwidth = bw2,
                      color = 'lightgrey',
                      linewidth = 0.1) +
  theme_ridges() + 
  labs(x = "Evolvability (%)",
       y = "Metric",
       title = "b) Selection for short trunks, short tails") +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size, angle = 45),
        axis.text.y = element_text(size = y_text_size),
        axis.title = element_text(size = axis_label_size),
        plot.title.position = "plot",
        plot.title = element_text(size=title_label_size, face = "plain", hjust = 0)) +
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
                     limits = lmtsB,
                     labels = ~.x*100) +
  scale_fill_manual(values = c("e" = vc_color[1,], 
                               "r" = vc_color[2,],
                               "c" = vc_color[3,])) +
  scale_y_discrete(labels = c("e" = expression(paste("e(",beta,")")), 
                              "r" = expression(paste("r(",beta,")")),
                              "c" = expression(paste("c(",beta,")")))) + 
  scale_color_manual(values = c("e" = vc_color[1,], 
                                "r" = vc_color[2,],
                                "c" = vc_color[3,]))

## Panel C ----
d <- e.params_beta |>
  filter(Beta_index == "3") |>
  select("e", "r", "c") |>
  pivot_longer(
    cols = c(e, r, c),
    names_to = "name",
    values_to = "value") |>
  mutate(name = factor(name, levels = c("e","r","c")))

summaries <- d |>
  group_by(name) |>
  summarise(summary_values = list(vecSmry(value)), .groups = 'drop') |>
  unnest_wider(summary_values, names_repair = "unique")


panelC <- ggplot(d, 
                 aes(x = value, y = name, fill = name)) +
  geom_density_ridges(scale = sc,
                      alpha = alp,
                      bandwidth = bw3,
                      color = 'lightgrey',
                      linewidth = 0.1) +
  theme_ridges() + 
  labs(x = "Evolvability (%)",
       y = "Metric",
       title = "c) Selection for long trunks, short tails") +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size, angle = 45),
        axis.text.y = element_text(size = y_text_size),
        axis.title = element_text(size = axis_label_size),
        plot.title.position = "plot",
        plot.title = element_text(size=title_label_size, face = "plain", hjust = 0)) +
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
  scale_x_continuous(breaks = brksC,
                     limits = lmtsC,
                     labels = ~.x*100) +
  scale_fill_manual(values = c("e" = vc_color[1,], 
                               "r" = vc_color[2,],
                               "c" = vc_color[3,])) +
  scale_y_discrete(labels = c("e" = expression(paste("e(",beta,")")), 
                              "r" = expression(paste("r(",beta,")")),
                              "c" = expression(paste("c(",beta,")")))) + 
  scale_color_manual(values = c("e" = vc_color[1,], 
                                "r" = vc_color[2,],
                                "c" = vc_color[3,]))


fig5 <- gridExtra::grid.arrange(panelA,
                                panelB,
                                panelC,
                                nrow = 1,
                                ncol = 3)
ggsave("Figures and Tables/Figure 5.pdf", 
       plot = fig5, 
       height = 2, 
       width = 6)
dev.off()
