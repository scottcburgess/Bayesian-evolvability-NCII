rm(list=ls())
library('dplyr')

# Table 1 ----
# Load data
head_tail_data <- readRDS("1_Data/head_tail_data.rds") 
hatch_settle_data <- readRDS("1_Data/hatch_settle_data.rds") 

# Summarize trunk and tail measurements
# No. blocks, No. sires, No. dams, No. full sib families
head_tail_data %>% 
  summarise(n.blocks=n_distinct(block),
            n.sire=n_distinct(sire),
            n.dam=n_distinct(dam),
            n.fam=n_distinct(interaction))
# Mean (±sd) no. of measurements per full sib family
head_tail_data %>% 
  group_by(interaction) %>% 
  summarize(n.trunk=length(head),
            n.tail=length(tail),
            mean.trunk=mean(n.trunk)) %>% 
  summarize(mean.trunk=mean(n.trunk),
            stdev.trunk=sd(n.trunk),
            total.trunk=sum(n.trunk))

head_tail_data %>% 
  group_by(block) %>%
  summarise(
    n.sire=n_distinct(sire),
    n.dam=n_distinct(dam),
    n.fam=n_distinct(interaction))

# Summarize hatching and settlement measurements
# No. blocks, No. sires, No. dams, No. full sib families
hatch_settle_data %>% 
  group_by(metric) %>% 
  summarise(n.blocks=n_distinct(block),
            n.sire=n_distinct(sire),
            n.dam=n_distinct(dam),
            n.fam=n_distinct(interaction))
# Mean (±sd) no. of measurements per full sib family
hatch_settle_data %>% 
  group_by(interaction, metric) %>% 
  count() %>% 
  group_by(metric) %>% 
  summarize(mean=mean(n),
            stdev=sd(n),
            total.n=sum(n))


hatch_settle_data %>% 
  group_by(block) %>%
  summarise(
    n.sire=n_distinct(sire),
    n.dam=n_distinct(dam),
    n.fam=n_distinct(interaction))

