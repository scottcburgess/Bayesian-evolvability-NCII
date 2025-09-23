rm(list = ls())
library(tidyverse)
library(runjags)


# Model III ---------------------------------------------------------------

load('Model_outputs/Model_III_posterior_20250923_0204.rdata')

smry.iii <- summary(
  post,
  vars = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov',
    'pr.overall', 'pr.block'
  ) 
) 
view(smry.iii)

plot(
  post, 
  vars = c(
    'sire.vcov', 'dam.vcov', 'interaction.vcov'
  ) ,
  plot.type = 'trace',
  separate.chains = TRUE
)

df |> 
  group_by(sire, metric) |> 
  summarize(pr = mean(outcome), .groups = 'drop') |> 
  ggplot() +
  geom_histogram(aes(pr)) +
  facet_wrap(~ metric, ncol = 1)




# Model IV ----------------------------------------------------------------

load('Model_out

smry.iv <- summary(
  post,
  vars = c(
    'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov'
  ) 
) 
view(smry.iv)
