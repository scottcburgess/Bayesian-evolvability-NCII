rm(list=ls())
library('tidyverse')


# Load data ----
df_head_tail <- readRDS("1_Data/head_tail_data.rds") 
df_hatch_settle <- readRDS("1_Data/hatch_settle_data.rds") 

# Prepare data ----
# Head - Tail
# Remove the fixed effects of block, 
# then calculate residuals, for plotting
df_head_tail <- df_head_tail |> 
  group_by(block) |> 
  reframe(sire = sire,
          head.resid = head - mean(head, na.rm=T),
          tail.resid = tail - mean(tail, na.rm=T))

# Calculate the sire averages
sire_means_head_tail <- df_head_tail |> 
  group_by(sire) |> 
  summarize(mean.head.resid = mean(head.resid),
            mean.tail.resid = mean(tail.resid))

# Hatch - Settle
# Calculate the family averages
family_means_hatch_settle <- df_hatch_settle |> 
  group_by(block, interaction, metric) |> 
  summarize(mean = mean(outcome)) |> 
  pivot_wider(names_from = metric,
              values_from = mean)

# Remove block effects 
family_means_hatch_settle <- family_means_hatch_settle |> 
  group_by(block) |> 
  reframe(hatch.resid = hatching - mean(hatching, na.rm=T),
          settle.resid = settling - mean(settling, na.rm=T))

# Make plot ----
panelA <- ggplot() +
  geom_point(data = df_head_tail,
             aes(x = head.resid,
                 y = tail.resid),
             alpha = 0.1) +
  # geom_point(data = sire_means_head_tail,
  #            aes(x = mean.head.resid,
  #                y = mean.tail.resid),
  #            color="blue",
  #            alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "Trunk length\n(residual)",
       y = "Tail length\n(residual)",
       title = "a)") +
  theme_bw()


panelB <- ggplot() +
  geom_point(data = family_means_hatch_settle,
             aes(x = hatch.resid,
                 y = settle.resid),
             alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "Hatching proportion\n(residual)",
       y = "Settling proportion\n(residual)",
       title = "b)") +
  theme_bw()


# Save plot ----
fig2 <- gridExtra::grid.arrange(panelA,
                                panelB,
                                nrow = 1,
                                ncol = 2)
ggsave("5_Figure_outputs/Figure 2.pdf", 
       plot = fig2, 
       height = 2.5, 
       width = 5)
dev.off()

