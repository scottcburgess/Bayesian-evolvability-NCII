library('dplyr')
source('0_misc_funcs.R')

### Table 1 ###

# Load data
trunk_tail <- readRDS("head_tail_data.rds") 
hatch_settle <- readRDS("hatch_settle_data.rds") 

# Summarize trunk and tail measurements
# No. blocks, No. sires, No. dams, No. full sib families
trunk_tail %>% 
  summarise(n.blocks=n_distinct(block),
            n.sire=n_distinct(sire),
            n.dam=n_distinct(dam),
            n.fam=n_distinct(interaction))
# Mean (±sd) no. of measurements per full sib family
trunk_tail %>% 
  group_by(interaction) %>% 
  summarize(n.trunk=length(head),
            n.tail=length(tail),
            mean.trunk=mean(n.trunk)) %>% 
  summarize(mean.trunk=mean(n.trunk),
            stdev.trunk=sd(n.trunk),
            total.trunk=sum(n.trunk))

# Summarize hatching and settlement measurements
# No. blocks, No. sires, No. dams, No. full sib families
hatch_settle %>% 
  group_by(metric) %>% 
  summarise(n.blocks=n_distinct(block),
                  n.sire=n_distinct(sire),
                  n.dam=n_distinct(dam),
                  n.fam=n_distinct(interaction))
# Mean (±sd) no. of measurements per full sib family
hatch_settle %>% 
  group_by(interaction, metric) %>% 
  count() %>% 
  group_by(metric) %>% 
  summarize(mean=mean(n),
            stdev=sd(n),
            total.n=sum(n))
###############



### Table 2 ###
options(digits=3,scipen=999)
load("Model_I_posterior_20230116_1944.rdata") 
Head.H <- vecSmry(p$H[1,1,])
Tail.H <- vecSmry(p$H[2,2,])
Head.E <- vecSmry(p$E[1,])
Tail.E <- vecSmry(p$E[2,])

load("Model_II_posterior_20230117_0603.rdata") 
Head.Tail.H <- vecSmry(p$H)
Head.Tail.E <- vecSmry(p$E)

load("Model_III_posterior_20230119_1441.rdata") 
Hatch.H <- vecSmry(p$H[1,])
Settle.H <- vecSmry(p$H[2,])
Hatch.E <- vecSmry(p$E[1,])
Settle.E <- vecSmry(p$E[2,])

cbind.data.frame(Head.H,
                 Tail.H,
                 Head.Tail.H,
                 Hatch.H,
                 Settle.H)

cbind.data.frame(Head.E,
                 Tail.E,
                 Head.Tail.E,
                 Hatch.E,
                 Settle.E)

###############