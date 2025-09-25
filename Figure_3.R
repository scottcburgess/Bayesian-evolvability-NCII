rm(list=ls())
library('tidyverse')
library(ggridges)

source('0_misc_funcs.R')


# Trunk Tail: load and prepare ----
load("Model_outputs/Model_I_posterior_20250923_2343.rdata") 

# quick check 
# vecSmry(p$H[1,1,]) # is ~similar to
# vecSmry(p$VA[1,1,] / p$VP[1,1,])
# vecSmry(vcv.obs$VA$vcv.G.obs[1,1,] / vcv.obs$VA$vcv.P.obs[1,1,])

# Prepare data 
tmp <- data.frame(
  TrunkVA = p$VA[1,1,] / p$VP[1,1,],
  TailVA = p$VA[2,2,] / p$VP[2,2,],
  TrunkVM = p$VM[1,1,] / p$VP[1,1,],
  TailVM = p$VM[2,2,] / p$VP[2,2,],
  TrunkVD = p$VD[1,1,] / p$VP[1,1,],
  TailVD = p$VD[2,2,] / p$VP[2,2,],
  # Residual = VP - (VA + VM + VD)
  TrunkVR = (p$VP[1,1,] - (p$VA[1,1,] + p$VM[1,1,] + p$VD[1,1,])) / p$VP[1,1,],
  TailVR = (p$VP[2,2,] - (p$VA[2,2,] + p$VM[2,2,] + p$VD[2,2,])) / p$VP[2,2,]
  )

d_Trunktail <- tmp |> 
  pivot_longer(cols = everything(), cols_vary = 'slowest') |> 
  as.data.frame()




# Trunk Tail ratio: load and prepare ----
load("Model_outputs/Model_II_posterior_20250921_0942.rdata") 

## Prepare data
tmp <- data.frame(
  ratioVA = var.obs$VA$var.a.obs / var.obs$VA$var.obs,
  ratioVM = var.obs$VM$var.a.obs / var.obs$VM$var.obs,
  ratioVD = var.obs$VD$var.a.obs / var.obs$VD$var.obs)

d_Trunktail_ratio <- tmp |> 
  pivot_longer(cols = everything(), cols_vary = 'slowest') |> 
  as.data.frame()




# Hatching Settling : load and prepare ----
load("Model_outputs/Model_III_posterior_20250923_0204.rdata") 

# quick check for hatching
# vecSmry(p$H[1,]) # is ~similar to
# vecSmry(vcv.obs$VA$vcv.G.obs[1,1,] / vcv.obs$VA$vcv.P.obs[1,1,])

## Prepare data 
tmp <- data.frame(
  HatchVA = vcv.obs$VA$vcv.G.obs[1,1,] / vcv.obs$VA$vcv.P.obs[1,1,],
  SettleVA = vcv.obs$VA$vcv.G.obs[2,2,] / vcv.obs$VA$vcv.P.obs[2,2,],
  HatchVM = vcv.obs$VM$vcv.G.obs[1,1,] / vcv.obs$VM$vcv.P.obs[1,1,],
  SettleVM = vcv.obs$VM$vcv.G.obs[2,2,] / vcv.obs$VM$vcv.P.obs[2,2,],
  HatchVD = vcv.obs$VD$vcv.G.obs[1,1,] / vcv.obs$VD$vcv.P.obs[1,1,],
  SettleVD = vcv.obs$VD$vcv.G.obs[2,2,] / vcv.obs$VD$vcv.P.obs[2,2,]
  )


d_hatchsettle <- tmp |> 
  pivot_longer(cols = everything(), cols_vary = 'slowest') |> 
  as.data.frame()





# Make Figure ----
vc_color <- data.frame(color = c("grey",
                                 "#2a9d8f",
                                 "#E9C46A",
                                 "#E76F51"))

## Plotting parameters ----
brks <- seq(0,1,0.1)
lmts <- c(0,0.8)
x_text_size <- 6
y_text_size <- 8
axis_label_size <- 10
title_label_size <- 10
point_size <- 2
segment_size <- 0.5
alp <- 0.4
bw <- 0.005
sc <- 0.9

## Panel A ----
names_to_use <- c("TrunkVA", "TrunkVM", "TrunkVD", "TrunkVR")

d <- d_Trunktail |>
  filter(name %in% names_to_use) |> 
  mutate(name = factor(name, levels = rev(names_to_use)))

summaries <- d |>
  group_by(name) |>
  summarise(summary_values = list(vecSmry(value)), .groups = 'drop') |>
  unnest_wider(summary_values, names_repair = "unique")


panelA <- ggplot(d, 
                 aes(x = value, y = name, fill = name)) +
  geom_density_ridges(scale = sc,
                      alpha = alp,
                      bandwidth = bw,
                      color = 'lightgrey',
                      linewidth = 0.1) +
  theme_ridges() + 
  labs(x = "Proportion of phenotypic variance",
       y = "Variance component",
       title = "a) Trunk length") +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size),
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
  scale_x_continuous(breaks = brks,
                     limits = lmts) +
  scale_fill_manual(values = c("TrunkVA" = vc_color[4,], 
                               "TrunkVM" = vc_color[3,],
                               "TrunkVD" = vc_color[2,],
                               "TrunkVR" = vc_color[1,])) +
  scale_y_discrete(labels = c("TrunkVA" = expression(V[A]), 
                              "TrunkVM" = expression(V[M]),
                              "TrunkVD" = expression(V[D]),
                              "TrunkVR" = expression(V[R]))) + 
  scale_color_manual(values = c("TrunkVA" = vc_color[4,], 
                                "TrunkVM" = vc_color[3,],
                                "TrunkVD" = vc_color[2,],
                                "TrunkVR" = vc_color[1,]))

## Panel B ----
names_to_use <- c("TailVA", "TailVM", "TailVD", "TailVR")

d <- d_Trunktail |>
  filter(name %in% names_to_use) |> 
  mutate(name = factor(name, levels = rev(names_to_use)))

summaries <- d |>
  group_by(name) |>
  summarise(summary_values = list(vecSmry(value)), .groups = 'drop') |>
  unnest_wider(summary_values, names_repair = "unique")


panelB <- ggplot(d, 
                 aes(x = value, y = name, fill = name)) +
  geom_density_ridges(scale = sc,
                      alpha = alp,
                      bandwidth = bw,
                      color = 'lightgrey',
                      linewidth = 0.1) +
  # stat_density_ridges(quantile_lines = TRUE, 
  #                     alpha = alp,
  #                     scale = sc,
  #                     quantiles = c(0.05, 0.5, 0.95)) +
  theme_ridges() + 
  labs(x = "Proportion of phenotypic variance",
       y = "Variance component",
       title = "b) Tail length") +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size),
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
  scale_x_continuous(breaks = brks,
                     limits = lmts) +
scale_fill_manual(values = c("TailVA" = vc_color[4,], 
                               "TailVM" = vc_color[3,],
                               "TailVD" = vc_color[2,],
                               "TailVR" = vc_color[1,])) +
  scale_y_discrete(labels = c("TailVA" = expression(V[A]), 
                              "TailVM" = expression(V[M]),
                              "TailVD" = expression(V[D]),
                              "TailVR" = expression(V[R]))) + 
  scale_color_manual(values = c("TailVA" = vc_color[4,], 
                                "TailVM" = vc_color[3,],
                                "TailVD" = vc_color[2,],
                                "TailVR" = vc_color[1,]))

## Panel C ----
names_to_use <- c("ratioVA", "ratioVM", "ratioVD")

d <- d_Trunktail_ratio |>
  filter(name %in% names_to_use) |> 
  mutate(name = factor(name, levels = rev(names_to_use)))

summaries <- d |>
  group_by(name) |>
  summarise(summary_values = list(vecSmry(value)), .groups = 'drop') |>
  unnest_wider(summary_values, names_repair = "unique")


panelC <- ggplot(d, 
                 aes(x = value, y = name, fill = name)) +
  geom_density_ridges(scale = sc,
                      alpha = alp,
                      bandwidth = bw,
                      color = 'lightgrey',
                      linewidth = 0.1) +
  theme_ridges() + 
  labs(x = "Proportion of phenotypic variance",
       y = "Variance component",
       title = "c) Trunk:Tail ratio") +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size),
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
  scale_x_continuous(breaks = brks,
                     limits = lmts) +
  scale_fill_manual(values = c("ratioVA" = vc_color[4,], 
                               "ratioVM" = vc_color[3,],
                               "ratioVD" = vc_color[2,])) +
  scale_y_discrete(labels = c("ratioVA" = expression(V[A]), 
                              "ratioVM" = expression(V[M]),
                              "ratioVD" = expression(V[D]))) + 
  scale_color_manual(values = c("ratioVA" = vc_color[4,], 
                                "ratioVM" = vc_color[3,],
                                "ratioVD" = vc_color[2,]))


## Panel D ----
names_to_use <- c("HatchVA", "HatchVM", "HatchVD")

d <- d_hatchsettle |>
  filter(name %in% names_to_use) |> 
  mutate(name = factor(name, levels = rev(names_to_use)))

summaries <- d |>
  group_by(name) |>
  summarise(summary_values = list(vecSmry(value)), .groups = 'drop') |>
  unnest_wider(summary_values, names_repair = "unique")


panelD <- ggplot(d, 
                 aes(x = value, y = name, fill = name)) +
  geom_density_ridges(scale = sc,
                      alpha = alp,
                      bandwidth = bw,
                      color = 'lightgrey',
                      linewidth = 0.1) +
  theme_ridges() + 
  labs(x = "Proportion of phenotypic variance",
       y = "Variance component",
       title = "d) Hatching probability") +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size),
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
  scale_x_continuous(breaks = brks,
                     limits = lmts) +
  scale_fill_manual(values = c("HatchVA" = vc_color[4,], 
                               "HatchVM" = vc_color[3,],
                               "HatchVD" = vc_color[2,])) +
  scale_y_discrete(labels = c("HatchVA" = expression(V[A]), 
                              "HatchVM" = expression(V[M]),
                              "HatchVD" = expression(V[D]))) + 
  scale_color_manual(values = c("HatchVA" = vc_color[4,], 
                                "HatchVM" = vc_color[3,],
                                "HatchVD" = vc_color[2,]))
  
                                
## Panel E ----
names_to_use <- c("SettleVA", "SettleVM", "SettleVD")

d <- d_hatchsettle |>
  filter(name %in% names_to_use) |> 
  mutate(name = factor(name, levels = rev(names_to_use)))

summaries <- d |>
  group_by(name) |>
  summarise(summary_values = list(vecSmry(value)), .groups = 'drop') |>
  unnest_wider(summary_values, names_repair = "unique")


panelE <- ggplot(d, 
                 aes(x = value, y = name, fill = name)) +
  geom_density_ridges(scale = sc,
                      alpha = alp,
                      bandwidth = bw,
                      color = 'lightgrey',
                      linewidth = 0.1) +
  theme_ridges() + 
  labs(x = "Proportion of phenotypic variance",
       y = "Variance component",
       title = "e) Settlement probability") +
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(size = x_text_size),
        axis.text.y = element_text(size = y_text_size),
        axis.title = element_text(size = axis_label_size),
        plot.title = element_text(size=title_label_size,face = "plain")) +
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
  scale_x_continuous(breaks = brks,
                     limits = lmts) +
  scale_fill_manual(values = c("SettleVA" = vc_color[4,], 
                               "SettleVM" = vc_color[3,],
                               "SettleVD" = vc_color[2,])) +
  scale_y_discrete(labels = c("SettleVA" = expression(V[A]), 
                              "SettleVM" = expression(V[M]),
                              "SettleVD" = expression(V[D]))) + 
  scale_color_manual(values = c("SettleVA" = vc_color[4,], 
                                "SettleVM" = vc_color[3,],
                                "SettleVD" = vc_color[2,]))



# Save plot ----
fig3 <- gridExtra::grid.arrange(panelA,
                                panelB,
                                panelC,
                                panelD,
                                panelE,
                        nrow = 2,
                        ncol = 3)
ggsave("Figures and Tables/Figure 3.pdf", 
       plot = fig3, 
       height = 4, 
       width = 8)
dev.off()



# Summaries ----
tmp <- d_Trunktail |> filter(name=="TrunkVA") |> select("value")
round(vecSmry(tmp$value),4) # Trunk VA (proportion)
tmp <- d_Trunktail |> filter(name=="TrunkVM") |> select("value")
round(vecSmry(tmp$value),4) # Trunk VM (proportion)
tmp <- d_Trunktail |> filter(name=="TrunkVD") |> select("value")
round(vecSmry(tmp$value),4) # Trunk VM (proportion)



tmp <- d_Trunktail |> filter(name=="TailVA") |> select("value")
round(vecSmry(tmp$value),4) # Trunk VA (proportion)
tmp <- d_Trunktail |> filter(name=="TailVM") |> select("value")
round(vecSmry(tmp$value),4) # Trunk VM (proportion)
tmp <- d_Trunktail |> filter(name=="TailVD") |> select("value")
round(vecSmry(tmp$value),4) # Trunk VM (proportion)


tmp <- d_hatchsettle |> filter(name=="HatchVA") |> select("value")
round(vecSmry(tmp$value),4) # Hatch VA (proportion)
tmp <- d_hatchsettle |> filter(name=="HatchVM") |> select("value")
round(vecSmry(tmp$value),4) # Hatch VM (proportion)
tmp <- d_hatchsettle |> filter(name=="HatchVD") |> select("value")
round(vecSmry(tmp$value),4) # Hatch VD (proportion)

tmp <- d_hatchsettle |> filter(name=="SettleVA") |> select("value")
round(vecSmry(tmp$value),4) # Settle VA (proportion)
tmp <- d_hatchsettle |> filter(name=="SettleVM") |> select("value")
round(vecSmry(tmp$value),4) # Settle VM (proportion)
tmp <- d_hatchsettle |> filter(name=="SettleVD") |> select("value")
round(vecSmry(tmp$value),4) # Settle VD (proportion)


