rm(list=ls())
library('dplyr')
source('0_misc_funcs.R')

setwd('/Users/scottburgess/Google Drive/Data/GitHub/Bayesian-evolvability-NCII/')

# Model I Results ----
load("Model_I_posterior_20230116_1944.rdata") 

# Phenotypic correlation between trunk and tail length
round(vecSmry(p$VP[1,2,] / (sqrt(p$VP[1,1,]) * sqrt(p$VP[2,2,]))),2)

# Mean trunk length
round(vecSmry(p$mean.overall[1,]),2)
# Variance in trunk length
round(vecSmry(p$VP[1,1,]),2)

# Mean tail length
round(vecSmry(p$mean.overall[2,]),2)
# Variance in tail length
round(vecSmry(p$VP[2,2,]),2)


vecSmry(p$VM[1,1,]) / vecSmry(p$VA[1,1,])
vecSmry(p$VM[2,2,]) / vecSmry(p$VA[2,2,])


# Model II Results ----
load("Model_II_posterior_20230117_0903.rdata") 

# Mean overall
round(vecSmry(p$mean.overall),3)
exp(-0.791) * 100

round(vecSmry(vcv.obs$VA$var.a.obs / vcv.obs$VA$var.obs),4)
round(vecSmry(vcv.obs$VM$var.a.obs / vcv.obs$VM$var.obs),4)
round(vecSmry(vcv.obs$VD$var.a.obs / vcv.obs$VD$var.obs),4)



# Model III Results ----
df_hatch_settle <- readRDS("hatch_settle_data.rds") 

# load("Model_III_posterior_20230119_1741.rdata") 
# 
# vecSmry(p$pr.overall[1,])

# Mean settlement probability per family
# Calculate the family averages
family_means_hatch_settle <- df_hatch_settle %>% 
  group_by(block, interaction, metric) %>% 
  summarize(mean = mean(outcome)) %>% 
  pivot_wider(names_from = metric,
              values_from = mean)

# Phenotypic correlation between hatching and settling among families
with(family_means_hatch_settle, cor.test(hatching, settling))

# Average hatching per familty
mean(family_means_hatch_settle$hatching,na.rm=T)
quantile(family_means_hatch_settle$hatching,na.rm=T, c(0.025,0.975))

# Average settlement per familty
mean(family_means_hatch_settle$settling,na.rm=T)
quantile(family_means_hatch_settle$settling,na.rm=T, c(0.025,0.975))




# Evolvability ----
load("Model_I_posterior_20230116_1944.rdata") 
# Prepare data 
# Get the posterior samples by averaging across all of the random beta's
eB_posterior_head_tail <- apply(e.params_BetaMCMC$post.dist$eB, 1, mean)
rB_posterior_head_tail <- apply(e.params_BetaMCMC$post.dist$rB, 1, mean)
cB_posterior_head_tail <- apply(e.params_BetaMCMC$post.dist$cB, 1, mean)
round(vecSmry(eB_posterior_head_tail),5)
round(vecSmry(rB_posterior_head_tail),5)
round(vecSmry(cB_posterior_head_tail),5)

load("Model_III_posterior_20230119_1741.rdata") 
## Prepare data 
eB_posterior_hatch_settle <- apply(e.params_BetaMCMC$post.dist$eB, 1, mean)
rB_posterior_hatch_settle <- apply(e.params_BetaMCMC$post.dist$rB, 1, mean)
cB_posterior_hatch_settle <- apply(e.params_BetaMCMC$post.dist$cB, 1, mean)
round(vecSmry(eB_posterior_hatch_settle),5)
round(vecSmry(rB_posterior_hatch_settle),5)
round(vecSmry(cB_posterior_hatch_settle),5)



# Model IV ----



# Table 1 ----
# Load data
load("Model_I_posterior_20230116_1944.rdata") 
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





# Table 2 ----
# Head Tail
load("Model_I_posterior_20230116_1944.rdata") 
# VA
round(vecSmry(p$VA[1,1,]),2) # Head
round(vecSmry(p$VA[2,2,]),2) # Tail
round(vecSmry(p$VA[1,2,]),2) # Covariance
# VM
round(vecSmry(p$VM[1,1,]),2) # Head
round(vecSmry(p$VM[2,2,]),2) # Tail
round(vecSmry(p$VM[1,2,]),2) # Covariance
# VD
round(vecSmry(p$VD[1,1,]),2) # Head
round(vecSmry(p$VD[2,2,]),2) # Tail
round(vecSmry(p$VD[1,2,]),2) # Covariance
# VP
round(vecSmry(p$VP[1,1,]),2) # Head
round(vecSmry(p$VP[2,2,]),2) # Tail
round(vecSmry(p$VP[1,2,]),2) # Covariance

# Head Tail ratio
load("Model_II_posterior_20230117_0903.rdata") 
round(vecSmry(p$VA),5) 
round(vecSmry(p$VM),5) 
round(vecSmry(p$VD),5) 
round(vecSmry(p$VP),5) 

# Hatch - Settle
load("Model_III_posterior_20230119_1741.rdata") 
round(vecSmry(vcv.obs$VA$vcv.G.obs[1,1,]),4) # Hatch
round(vecSmry(vcv.obs$VA$vcv.G.obs[2,2,]),4) # Settle
round(vecSmry(vcv.obs$VA$vcv.G.obs[1,2,]),4) # Covariance

round(vecSmry(vcv.obs$VM$vcv.G.obs[1,1,]),4) # Hatch
round(vecSmry(vcv.obs$VM$vcv.G.obs[2,2,]),4) # Settle
round(vecSmry(vcv.obs$VM$vcv.G.obs[1,2,]),4) # Covariance

round(vecSmry(vcv.obs$VD$vcv.G.obs[1,1,]),4) # Hatch
round(vecSmry(vcv.obs$VD$vcv.G.obs[2,2,]),4) # Settle
round(vecSmry(vcv.obs$VD$vcv.G.obs[1,2,]),4) # Covariance

round(vecSmry(vcv.obs$VP$vcv.G.obs[1,1,]),4) # Hatch
round(vecSmry(vcv.obs$VP$vcv.G.obs[2,2,]),4) # Settle
round(vecSmry(vcv.obs$VP$vcv.G.obs[1,2,]),4) # Covariance





### Table 3 ###
options(digits=3,scipen=999)
load("Model_I_posterior_20230116_1944.rdata") 
Head.H <- vecSmry(p$H[1,1,])
Tail.H <- vecSmry(p$H[2,2,])
Head.E <- vecSmry(p$E[1,])
Tail.E <- vecSmry(p$E[2,])
# Other summaries for the results section
vecSmry(p$mean.overall[1,]) # Mean trunk length
# mean(df$head)
vecSmry(p$VP[1,1,]) # Trunk P variance

vecSmry(p$mean.overall[2,]) # Mean tail length
# mean(df$tail)
vecSmry(p$VP[2,2,]) # Tail P variance

vecSmry(p$VP[1,2,]) # Phenotypic covariance

vecSmry(p$VA[1,1,] / p$VP[1,1,]) * 100 # VA / VP
vecSmry(p$VM[1,1,] / p$VP[1,1,]) * 100 # VM / VP
vecSmry(p$VD[1,1,] / p$VP[1,1,]) * 100 # VD / VP

vecSmry(p$VA[2,2,] / p$VP[2,2,]) * 100 # VA / VP
vecSmry(p$VM[2,2,] / p$VP[2,2,]) * 100 # VM / VP
vecSmry(p$VD[2,2,] / p$VP[2,2,]) * 100 # VD / VP


load("Model_II_posterior_20230117_0903.rdata") 
Head.Tail.H <- vecSmry(p$H)
Head.Tail.E <- vecSmry(p$E)

# Other summaries for the results section
vecSmry(vcv.obs$VA$mean.obs) # Mean ratio
# mean(df$head/df$tail)
vecSmry(vcv.obs$VA$var.obs) # ratio P variance

vecSmry(vcv.obs$VA$var.a.obs / vcv.obs$VA$var.obs) * 100 # VA / VP
vecSmry(vcv.obs$VM$var.a.obs / vcv.obs$VM$var.obs) * 100 # VM / VP
vecSmry(vcv.obs$VD$var.a.obs / vcv.obs$VD$var.obs) * 100 # VD / VP


load("Model_III_posterior_20230119_1741.rdata") 
Hatch.H <- vecSmry(p$H[1,])
Settle.H <- vecSmry(p$H[2,])
Hatch.E <- vecSmry(p$E[1,])
Settle.E <- vecSmry(p$E[2,])

# Other summaries for the results section
# df %>% group_by(metric) %>% summarize(mean=mean(outcome))

vecSmry(vcv.obs$VA$mean.obs[,'hatching']) # Mean hatching probability
vecSmry(vcv.obs$VA$vcv.P.obs[1,1,]) # Hatching P variance

vecSmry(vcv.obs$VA$mean.obs[,'settling']) # Mean hatching probability
vecSmry(vcv.obs$VA$vcv.P.obs[2,2,]) # Settling P variance

vecSmry(vcv.obs$VA$vcv.P.obs[1,2,]) # Phenotypic covariance

vecSmry(vcv.obs$VA$vcv.G.obs[1,1,] / vcv.obs$VA$vcv.P.obs[1,1,]) * 100 # VA / VP
vecSmry(vcv.obs$VM$vcv.G.obs[1,1,] / vcv.obs$VM$vcv.P.obs[1,1,]) * 100 # VM / VP
vecSmry(vcv.obs$VD$vcv.G.obs[1,1,] / vcv.obs$VD$vcv.P.obs[1,1,]) * 100 # VD / VP

vecSmry(vcv.obs$VA$vcv.G.obs[2,2,] / vcv.obs$VA$vcv.P.obs[2,2,]) * 100 # VA / VP
vecSmry(vcv.obs$VM$vcv.G.obs[2,2,] / vcv.obs$VM$vcv.P.obs[2,2,]) * 100 # VM / VP
vecSmry(vcv.obs$VD$vcv.G.obs[2,2,] / vcv.obs$VD$vcv.P.obs[2,2,]) * 100 # VD / VP


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



# Load data
load('Model_IV_posterior_20230117_0038.rdata')

# Full sib families
vecSmry(p$int.beta[1,]) # head
vecSmry(p$int.beta[2,]) # tail

# Maternal effect
vecSmry(p$maternal.beta[1,]) # head 
vecSmry(p$maternal.beta[2,]) # tail

# Additive Sire 
vecSmry(p$sire.beta[1,]) # head 
vecSmry(p$sire.beta[2,]) # tail
