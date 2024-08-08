rm(list=ls())
library('dplyr')
source('0_misc_funcs.R')

# Load data
load("Model_III_posterior_20230119_1741.rdata") 

######## Summary stats of the posteriors for variance components ##############
summary_dat <- as.data.frame(rbind(vecSmry(p$additive.vcov[1,1,]),
                                   vecSmry(p$additive.vcov[2,2,]),
                                   vecSmry(p$additive.vcov[1,2,]),
                                   vecSmry(p$maternal.vcov[1,1,]),
                                   vecSmry(p$maternal.vcov[2,2,]),
                                   vecSmry(p$maternal.vcov[1,2,]),
                                   vecSmry(p$interaction.vcov[1,1,]),
                                   vecSmry(p$interaction.vcov[2,2,]),
                                   vecSmry(p$interaction.vcov[1,2,]),
                                   vecSmry(vcv.obs$VA$vcv.G.obs[1,1,]),
                                   vecSmry(vcv.obs$VA$vcv.G.obs[2,2,]),
                                   vecSmry(vcv.obs$VA$vcv.G.obs[1,2,]),
                                   vecSmry(vcv.obs$VM$vcv.G.obs[1,1,]),
                                   vecSmry(vcv.obs$VM$vcv.G.obs[2,2,]),
                                   vecSmry(vcv.obs$VM$vcv.G.obs[1,2,]),
                                   vecSmry(vcv.obs$VD$vcv.G.obs[1,1,]),
                                   vecSmry(vcv.obs$VD$vcv.G.obs[2,2,]),
                                   vecSmry(vcv.obs$VD$vcv.G.obs[1,2,]),
                                   vecSmry(p$H[1,]),
                                   vecSmry(p$H[2,]),
                                   vecSmry(p$E[1,]),
                                   vecSmry(p$E[2,])),
                             row.names=c('Hatch.G','Settle.G','Hatch.Settle.cG',
                                         'Hatch.M','Settle.M','Hatch.Settle.cM',
                                         'Hatch.D','Settle.D','Hatch.Settle.cD',
                                         'Hatch.VA','Settle.VA','Hatch.Settle.cVA',
                                         'Hatch.VM','Settle.VM','Hatch.Settle.cVM',
                                         'Hatch.VD','Settle.VD','Hatch.Settle.cVD',
                                         'Hatch.H','Settle.H',
                                         'Hatch.E','Settle.E'))
round(summary_dat,6)
############################################


################# Make Figure ################
vc.color <- data.frame(color = c("#E76F51",
                                 "#f0a794",
                                 "#E9C46A",
                                 "#f3dead",
                                 "#2a9d8f",
                                 "#82ded3"))

# Use line 118 or 119 just for viewing, but comment out when saving final pdf.
# When saving pdf for manuscript, uncomment lines 121 and 232
# windows(width=6,height=6) # use on PC
quartz(width=6,height=6) # use on Mac

# pdf('Figure Model III.pdf', width = 6, height = 6)

par(mfrow = c(3,3), mar=c(4,0,2,1),oma=c(1,3,2,1))

var_lims <- c(0,0.16)
covar_lims <- c(-0.05,0.05)
bar.width = 0.0001

# a)
plot.posterior(x=vcv.obs$VA$vcv.G.obs[1,1,],
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[1:2,])
mtext("a)",side=3,adj=0)
mtext("Hatching",side=3,adj=0.5,line=1)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,4,by=0.01))
mtext("Additive genetic variance, VA",side=1,line=2.5,cex=0.8)
# b) 
plot.posterior(x=vcv.obs$VA$vcv.G.obs[2,2,],
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[1:2,])
mtext("b)",side=3,adj=0)
mtext("Settlement",side=3,adj=0.5,line=1)
axis(side=1,at=seq(0,4,by=0.01))
mtext("Additive genetic variance, VA",side=1,line=2.5,cex=0.8)
# c) 
plot.posterior(x=vcv.obs$VA$vcv.G.obs[1,2,],
               bar.width=bar.width,
               xlims=covar_lims,
               cols=vc.color[1:2,])
mtext("c) ",side=3,adj=0)
mtext("Hatching-Settlement",side=3,adj=0.5,line=1)
axis(side=1,at=seq(-30,30,by=0.01))
mtext("Additive genetic covariance",side=1,line=2.5,cex=0.8)

# d)
plot.posterior(x=vcv.obs$VM$vcv.G.obs[1,1,],
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[3:4,])
mtext("d)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,4,by=0.01))
mtext("Maternal effect variance, VM",side=1,line=2.5,cex=0.8)
# e) 
plot.posterior(x=vcv.obs$VM$vcv.G.obs[2,2,],
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[3:4,])
mtext("e)",side=3,adj=0)
axis(side=1,at=seq(0,4,by=0.01))
mtext("Maternal effect variance, VM",side=1,line=2.5,cex=0.8)
# f) 
plot.posterior(x=vcv.obs$VM$vcv.G.obs[1,2,],
               bar.width=bar.width,
               xlims=covar_lims,
               cols=vc.color[3:4,])
mtext("f)",side=3,adj=0)
axis(side=1,at=seq(-30,30,by=0.01))
mtext("Maternal effect covariance",side=1,line=2.5,cex=0.8)

# g)
plot.posterior(x=vcv.obs$VD$vcv.G.obs[1,1,],
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[5:6,])
mtext("g)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,4,by=0.01))
mtext("Dominance variance, VD",side=1,line=2.5,cex=0.8)
# h) 
plot.posterior(x=vcv.obs$VD$vcv.G.obs[2,2,],
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[5:6,])
mtext("h)",side=3,adj=0)
axis(side=1,at=seq(0,4,by=0.01))
mtext("Dominance variance, VD",side=1,line=2.5,cex=0.8)
# i) 
plot.posterior(x=vcv.obs$VD$vcv.G.obs[1,2,],
               bar.width=bar.width,
               xlims=covar_lims,
               cols=vc.color[5:6,])
mtext("i)",side=3,adj=0)
axis(side=1,at=seq(-30,30,by=0.01))
mtext("Dominance covariance",side=1,line=2.5,cex=0.8)
############################

# dev.off()