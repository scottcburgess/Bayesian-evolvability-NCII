rm(list=ls())
library('dplyr')
source('0_misc_funcs.R')

# Load data
load("Model_I_posterior_20230116_1944.rdata") 
pE_head_tail <- p$E
e.params_BetaMCMC_head_tail <- e.params_BetaMCMC

load("Model_III_posterior_20230119_1741.rdata") 
pE_hatch_settle <- p$E
e.params_BetaMCMC_hatch_settle <- e.params_BetaMCMC

# Get the posterior samples by averaging across all of the random beta's
eB_posterior_head_tail <- apply(e.params_BetaMCMC_head_tail$post.dist$eB, 1, mean)
# Check
# e.params_BetaMCMC_head_tail$summary # median (called e_mean) should be the same as
# vecSmry(eB_posterior_head_tail) # median here (but we're using the mode)
# vecSmry(pE_head_tail[1,]) # Should be roughly similar to the sinlge trait estimate, since low covariance
rB_posterior_head_tail <- apply(e.params_BetaMCMC_head_tail$post.dist$rB, 1, mean)
cB_posterior_head_tail <- apply(e.params_BetaMCMC_head_tail$post.dist$cB, 1, mean)
eB_posterior_hatch_settle <- apply(e.params_BetaMCMC_hatch_settle$post.dist$eB, 1, mean)
# Check
# e.params_BetaMCMC_hatch_settle$summary # median (called e_mean) should be the same as
# vecSmry(eB_posterior_hatch_settle) # median here (but we're using the mode)
# vecSmry(pE_hatch_settle[1,]) # Should be roughly similar to the sinlge trait estimate, since low covariance
rB_posterior_hatch_settle <- apply(e.params_BetaMCMC_hatch_settle$post.dist$rB, 1, mean)
cB_posterior_hatch_settle <- apply(e.params_BetaMCMC_hatch_settle$post.dist$cB, 1, mean)

vecSmry(eB_posterior_hatch_settle / eB_posterior_head_tail)
vecSmry(rB_posterior_hatch_settle / rB_posterior_head_tail)
vecSmry(cB_posterior_hatch_settle / cB_posterior_head_tail)

################# Make Figure ################
vc.color <- data.frame(color = c("#0077b6",
                                 "#9dddff",
                                 "#00b4d8",
                                 "#73e8ff",
                                 "#90e0ef",
                                 "#d3f3f9"))

# Use line 118 or 119 just for viewing, but comment out when saving final pdf.
# When saving pdf for manuscript, uncomment lines 121 and 232
# windows(width=6,height=4) # use on PC
quartz(width=6,height=4) # use on Mac

# pdf('Figure Model Evolvability.pdf', width = 6, height = 4)

par(mfrow = c(2,3), mar=c(5,0,2,1),oma=c(1,3,2,1))


# Trunk - tail
e_lims <- c(0,0.002)
bar.width <- 0.00001
# a)
plot.posterior(x=eB_posterior_head_tail,
               bar.width=bar.width,
               xlims=e_lims,
               cols=vc.color[1:2,])
mtext("a)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,1,by=0.0001))
mtext(expression(paste("Evolvability, e",beta)),side=1,line=2.5,cex=0.8)
# b) 
plot.posterior(x=rB_posterior_head_tail,
               bar.width=bar.width,
               xlims=e_lims,
               cols=vc.color[3:4,])
mtext("b)",side=3,adj=0)
mtext("Trunk-Tail",side=3,adj=0.5,line=1.5)
axis(side=1,at=seq(0,1,by=0.0001))
mtext(expression(paste("Respondability, r",beta)),side=1,line=2.5,cex=0.8)
# c) 
plot.posterior(x=cB_posterior_head_tail,
               bar.width=bar.width,
               xlims=e_lims,
               cols=vc.color[5:6,])
mtext("c)",side=3,adj=0)
axis(side=1,at=seq(0,1,by=0.0001))
mtext(expression(paste("Conditional evolvability, r",beta)),side=1,line=2.5,cex=0.8)

# Hatch - settle
e_lims <- c(0,0.1)
bar.width <- 0.0005
# d)
plot.posterior(x=eB_posterior_hatch_settle,
               bar.width=bar.width,
               xlims=e_lims,
               cols=vc.color[1:2,])
mtext("d)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,2,by=0.01))
mtext(expression(paste("Evolvability, e",beta)),side=1,line=2.5,cex=0.8)
# e) 
plot.posterior(x=rB_posterior_hatch_settle,
               bar.width=bar.width,
               xlims=e_lims,
               cols=vc.color[3:4,])
mtext("e)",side=3,adj=0)
mtext("Hatch-Settle",side=3,adj=0.5,line=1.5)
axis(side=1,at=seq(0,2,by=0.01))
mtext(expression(paste("Respondability, r",beta)),side=1,line=2.5,cex=0.8)
# f) 
plot.posterior(x=cB_posterior_hatch_settle,
               bar.width=bar.width,
               xlims=e_lims,
               cols=vc.color[5:6,])
mtext("f)",side=3,adj=0)
axis(side=1,at=seq(0,2,by=0.01))
mtext(expression(paste("Conditional evolvability, r",beta)),side=1,line=2.5,cex=0.8)

# dev.off()
