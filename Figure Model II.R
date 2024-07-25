library('dplyr')
source('0_misc_funcs.R')

# Load data
load("Model_II_posterior_20230117_0903.rdata") 

######## Summary stats of the posteriors for variance components ##############
summary_dat <- as.data.frame(rbind(vecSmry(p$additive.var),
                                   vecSmry(p$maternal.var),
                                   vecSmry(p$interaction.var),
                                   vecSmry(p$resid.var),
                                   vecSmry(qgparams.post$var.a.obs),
                                   vecSmry(p$H),
                                   vecSmry(p$E)),
                             row.names=c('G',
                                         'M',
                                         'D',
                                         'res',
                                         'VA',
                                         'H',
                                         'E'))
summary_dat <- round(summary_dat,7)
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
# windows(width=3,height=6) # use on PC
quartz(width=3,height=6) # use on Mac

# pdf('Figure Model II.pdf', width = 3, height = 6)

par(mfrow = c(3,1), mar=c(3,0,2,0),oma=c(1,3,2,1))

var_lims <- c(0,0.003)
bar.width <- 0.00001

# a)
plot.posterior(x=p$additive,
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[1:2,],
               trait.plot="G",
               trait.summary="G")
mtext("a)",side=3,adj=0)
mtext("Trunk tail ratio ",side=3,adj=0.5,line=1)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,1,by=0.0005),line=0)
mtext("Additive sire variance, G",side=1,line=2.5,cex=0.8)
# b)
plot.posterior(x=p$maternal.var,
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[3:4,],
               trait.plot="M",
               trait.summary="M")
mtext("b)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,1,by=0.0005))
mtext("Maternal effect variance, M",side=1,line=2.5,cex=0.8)
# c)
plot.posterior(x=p$interaction.var,
               bar.width=bar.width,
               xlims=var_lims,
               cols=vc.color[5:6,],
               trait.plot="D",trait.summary="D")
mtext("c)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,1,by=0.0005),line=0)
mtext("Interaction variance, D",side=1,line=2.5,cex=0.8)
############################

# dev.off()