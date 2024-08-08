rm(list=ls())
library('dplyr')
source('0_misc_funcs.R')

# Load data
load("Model_II_posterior_20230117_0903.rdata") 

######## Summary stats of the posteriors for variance components ##############
summary_dat <- as.data.frame(rbind(vecSmry(p$additive.var),
                                   vecSmry(p$maternal.var),
                                   vecSmry(p$interaction.var),
                                   vecSmry(p$resid.var),
                                   vecSmry(var.obs$VA$var.a.obs),
                                   vecSmry(var.obs$VM$var.a.obs),
                                   vecSmry(var.obs$VD$var.a.obs),
                                   vecSmry(p$H),
                                   vecSmry(p$E)),
                             row.names=c('G',
                                         'M',
                                         'D',
                                         'res',
                                         'VA',
                                         'VM',
                                         'VD',
                                         'H',
                                         'E'))
round(summary_dat,7)
############################################


################# Make Figure ################
vc.color <- data.frame(color = c("#E76F51",
                                 "#f0a794",
                                 "#E9C46A",
                                 "#f3dead",
                                 "#2a9d8f",
                                 "#82ded3"))

# Use line 41 or 42 just for viewing, but comment out when saving final pdf.
# When saving pdf for manuscript, uncomment lines 44 and 86
# windows(width=3,height=6) # use on PC
quartz(width=3,height=6) # use on Mac

# pdf('Figure Model II.pdf', width = 3, height = 6)

par(mfrow = c(3,1), mar=c(3,0,2,0),oma=c(1,3,2,1))

var_lims <- c(0,0.0007)
bar.width <- 0.000003

# a)
plot.posterior(x=var.obs$VA$var.a.obs,
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[1:2,])
mtext("a)",side=3,adj=0)
mtext("Trunk tail ratio ",side=3,adj=0.5,line=1)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,1,by=0.00005),line=0)
mtext("Additive genetic variance, VA",side=1,line=2.5,cex=0.8)
# b)
plot.posterior(x=var.obs$VM$var.a.obs,
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[3:4,])
mtext("b)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,1,by=0.00005))
mtext("Maternal effect variance, VM",side=1,line=2.5,cex=0.8)
# c)
plot.posterior(x=var.obs$VD$var.a.obs,
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[5:6,])
mtext("c)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,1,by=0.00005),line=0)
mtext("Dominance variance, VD",side=1,line=2.5,cex=0.8)
############################

# dev.off()