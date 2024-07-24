library('dplyr')
source('0_misc_funcs.R')

# Load data
load("hatch_settle_posterior_20230119_1741.rdata") 

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
                                   vecSmry(qgparams.post$hatching$var.a.obs),
                                   vecSmry(qgparams.post$settling$var.a.obs),
                                   vecSmry(qgparams.post$hatching$h2.obs),
                                   vecSmry(qgparams.post$settling$h2.obs),
                                   vecSmry(qgparams.post$hatching$E),
                                   vecSmry(qgparams.post$settling$E)),
                             row.names=c('Hatch.G','Settle.G','Hatch.Settle.cG',
                                         'Hatch.M','Settle.M','Hatch.Settle.cM',
                                         'Hatch.D','Settle.D','Hatch.Settle.cD',
                                         'Hatch.VA','Settle.VA',
                                         'Hatch.H','Settle.H',
                                         'Hatch.E','Settle.E'))
summary_dat <- round(summary_dat,4)
############################################




####### Function to make plot of posterior
plot.posterior <- function(x,bar.width,xlims,cols,trait.plot,trait.summary){
  # x <- p$additive.vcov[1,1,] # data to plot
  # bar.width = 0.01 # width of the bars in the histogram
  # xlims <- c(0,1.5)
  # cols = "dodgerblue" # color of the bars to plot
  # trait.plot = "Hatch.G" # What trait to plot from summary_dat
  # trait.summary = "Hatch.G" # What trait to summarize in legend
  blims <- max(abs(x))+10
  brks <- seq(-blims,blims,by=bar.width)
  foo <- hist(x, breaks=brks,plot=F)
  d <- data.frame(x=foo$mids, y=foo$density)
  d <- d %>% filter(y>0)
  ylims <- c(0,max(d$y))
  
  plot(xlims,ylims,type="n",bty="l",ylab="",xlab="",axes=F)

  interval <- as.numeric(summary_dat[trait.plot,c("lower.hdi","upper.hdi")])
  
  with(d,segments(x,
                  rep(0,length(x)),
                  x,
                  y,
                  lend=2,
                  col=ifelse(d$x > interval[1] & d$x < interval[2],
                             cols[1],cols[2])))
  
  abline(v=0)
# Add summaries to plot
  median.text <- eval(paste("median = ",summary_dat[trait.summary,"median"]))  
  mode.text <- eval(paste("mode = ",summary_dat[trait.summary,"mode"]))  
  int <- summary_dat[trait.summary, c("lower.hdi","upper.hdi")]
  int.text <- eval(paste("95% hdi = ",int[1], " - ", int[2],sep=""))  
  
  # title <- ifelse(grepl("G",trait.summary),"G",
  #                 ifelse(grepl("M",trait.summary),"M","D"))
  

  legend(xlims[2]*0.8,ylims[2],
         legend=rbind(median.text,mode.text,int.text),
         bty="n",
         adj=1,
         # title=title,
         # title.adj=0
         )
}
###################################################





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

hatch_lims <- c(0,3)
settle_lims <- c(0,3)
covar_lims <- c(-0.9,0.9)
bar.width = 0.01

# a)
plot.posterior(x=p$additive.vcov[1,1,],
               bar.width=bar.width,
               xlims=hatch_lims,
               cols=vc.color[1:2,],
               trait.plot="Hatch.G",
               trait.summary="Hatch.G")
mtext("a)",side=3,adj=0)
mtext("Hatching",side=3,adj=0.5,line=1)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,20,by=0.5))
mtext("Additive sire variance, G",side=1,line=2.5,cex=0.8)
# b) 
plot.posterior(x=p$additive.vcov[2,2,],
               bar.width=bar.width,
               xlims=settle_lims,
               cols=vc.color[1:2,],
               trait.plot="Settle.G",
               trait.summary="Settle.G")
mtext("b)",side=3,adj=0)
mtext("Settlment",side=3,adj=0.5,line=1)
axis(side=1,at=seq(0,20,by=0.5))
mtext("Additive sire variance, G",side=1,line=2.5,cex=0.8)
# c) 
plot.posterior(x=p$additive.vcov[1,2,],
               bar.width=bar.width,
               xlims=covar_lims,
               cols=vc.color[1:2,],
               trait.plot="Hatch.Settle.cG",
               trait.summary="Hatch.Settle.cG")
mtext("c) ",side=3,adj=0)
mtext("Hatching-Settlement",side=3,adj=0.5,line=1)
axis(side=1,at=seq(-30,30,by=0.2))
mtext("Additive sire covariance",side=1,line=2.5,cex=0.8)

# d)
plot.posterior(x=p$maternal.vcov[1,1,],
               bar.width=bar.width,
               xlims=hatch_lims,
               cols=vc.color[3:4,],
               trait.plot="Hatch.M",
               trait.summary="Hatch.M")
mtext("d)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,20,by=0.5))
mtext("Maternal effect variance, M",side=1,line=2.5,cex=0.8)
# e) 
plot.posterior(x=p$maternal.vcov[2,2,],
               bar.width=bar.width,
               xlims=settle_lims,
               cols=vc.color[3:4,],
               trait.plot="Settle.M",
               trait.summary="Settle.M")
mtext("e)",side=3,adj=0)
axis(side=1,at=seq(0,20,by=0.5))
mtext("Maternal effect variance, M",side=1,line=2.5,cex=0.8)
# f) 
plot.posterior(x=p$maternal.vcov[1,2,],
               bar.width=bar.width,
               xlims=covar_lims,
               cols=vc.color[3:4,],
               trait.plot="Hatch.Settle.cM",
               trait.summary="Hatch.Settle.cM")
mtext("f)",side=3,adj=0)
axis(side=1,at=seq(-30,30,by=0.2))
mtext("Maternal effect covariance",side=1,line=2.5,cex=0.8)

# g)
plot.posterior(x=p$interaction.vcov[1,1,],
               bar.width=bar.width,
               xlims=hatch_lims,
               cols=vc.color[5:6,],
               trait.plot="Hatch.D",
               trait.summary="Hatch.D")
mtext("g)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,20,by=0.5))
mtext("Interaction variance, D",side=1,line=2.5,cex=0.8)
# h) 
plot.posterior(x=p$interaction.vcov[2,2,],
               bar.width=bar.width,
               xlims=settle_lims,
               cols=vc.color[5:6,],
               trait.plot="Settle.D",
               trait.summary="Settle.D")
mtext("h)",side=3,adj=0)
axis(side=1,at=seq(0,20,by=0.5))
mtext("Interaction variance, D",side=1,line=2.5,cex=0.8)
# i) 
plot.posterior(x=p$interaction.vcov[1,2,],
               bar.width=bar.width,
               xlims=covar_lims,
               cols=vc.color[5:6,],
               trait.plot="Hatch.Settle.cD",
               trait.summary="Hatch.Settle.cD")
mtext("i)",side=3,adj=0)
axis(side=1,at=seq(-30,30,by=0.2))
mtext("Interaction covariance",side=1,line=2.5,cex=0.8)
############################

# dev.off()