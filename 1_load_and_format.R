rm(list = ls())
library(tidyverse)

df <- read.csv("Moccidentalis All Traits.csv") %>% 
  mutate(
    interaction = paste0("B", block, ".S", sire, "xD", dam),
    trait = ifelse(trait == "h.length", "head", trait),
    trait = ifelse(trait == "n.length", "tail", trait)
  ) %>% 
  pivot_wider(names_from = "trait", values_from = "value") %>% 
  # only keep blocks with more than 1 sire and dam
  group_by(block) %>% 
  mutate(to.keep = n_distinct(sire) > 1 & n_distinct(dam) > 1) %>%
  ungroup() %>% 
  filter(to.keep) %>% 
  select(-to.keep)  

df %>%
  group_by(block) %>% 
  summarize(
    num.larvae = n_distinct(animal),
    num.sires = n_distinct(sire),
    num.dams = n_distinct(dam),
    num.lengths = sum(!is.na(head) & !is.na(tail)),
    num.hatching = sum(!is.na(hatching)),
    num.settling = sum(!is.na(settling)),
    .groups = "drop"
  )

df %>% 
  select(animal, block, interaction, sire, dam, head, tail) %>% 
  filter(!(is.na(head) & is.na(tail))) %>% 
  saveRDS("head_tail_data.rds")

df %>% 
  select(animal, block, interaction, sire, dam, hatching, settling) %>% 
  filter(!(is.na(hatching) & is.na(settling))) %>% 
  pivot_longer(c("hatching", "settling"), names_to = "metric", values_to = "outcome") %>% 
  filter(!is.na(outcome)) %>% 
  saveRDS("hatch_settle_data.rds")