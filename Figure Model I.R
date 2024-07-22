library('dplyr')

# Load data
load("head_tail_posterior_20230116_1944.rdata") 

# Function to summarize posteriors
vecSmry <- function(x) {      
  setNames(
    c(mean(x),median(x), modeest::mlv(x, method = "venter"), HDInterval::hdi(x)),
    c("mean", "median", "mode", "lower.hdi", "upper.hdi")
  )[c("lower.hdi", "mean", "median", "mode", "upper.hdi")]
}

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
                                   vecSmry(p$resid.vcov[1,1,]),
                                   vecSmry(p$resid.vcov[2,2,]),
                                   vecSmry(p$resid.vcov[1,2,]),
                                   vecSmry(p$VA[1,1,]),
                                   vecSmry(p$VA[2,2,]),
                                   vecSmry(p$VA[1,2,]),
                                   vecSmry(p$VM[1,1,]),
                                   vecSmry(p$VM[2,2,]),
                                   vecSmry(p$VM[1,2,]),
                                   vecSmry(p$VD[1,1,]),
                                   vecSmry(p$VD[2,2,]),
                                   vecSmry(p$VD[1,2,]),
                                   vecSmry(p$VP[1,1,]),
                                   vecSmry(p$VP[2,2,]),
                                   vecSmry(p$VP[1,2,]),
                                   vecSmry(p$H[1,1,]),
                                   vecSmry(p$H[2,2,]),
                                   vecSmry(p$E[1,]),
                                   vecSmry(p$E[2,])),
                             row.names=c('Head.G','Tail.G','Head.Tail.cG',
                                         'Head.M','Tail.M','Head.Tail.cM',
                                         'Head.D','Tail.D','Head.Tail.cD',
                                         'Head.res','Tail.res','Head.Tail.cres',
                                         'Head.VA','Tail.VA','Head.Tail.cVA',
                                         'Head.VM','Tail.VM','Head.Tail.cVM',
                                         'Head.VD','Tail.VD','Head.Tail.cVD',
                                         'Head.VP','Tail.VP','Head.Tail.cVP',
                                         'Head.H','Tail.H',
                                         'Head.E','Tail.E'))
summary_dat <- round(summary_dat,2)
############################################




####### Function to make plot of posterior
plot.posterior <- function(x,bar.width,xlims,cols,trait.plot,trait.summary){
  # For testing
  # x <- p$additive.vcov[1,1,] # data to plot
  # bar.width = 0.5 # width of the bars in the histogram
  # xlims <- c(0,120)
  # cols = "dodgerblue" # color of the bars to plot
  # trait.plot = "Head.G" # What trait to plot from summary_dat (e.g., G or VA)
  # trait.summary = "Head.G" # What trait to summarize in legend (e.g., G or VA)
  blims <- max(abs(x))*2
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
  median.text <- eval(paste("median = ",summary_dat[trait.summary, "median"]))
  mode.text <- eval(paste("mode = ",summary_dat[trait.summary,"mode"]))  
  int <- summary_dat[trait.summary,c("lower.hdi","upper.hdi")]
  int.text <- eval(paste("95% hdi = ",int[1], " - ", int[2],sep=""))  
  
  # title <- ifelse(grepl("VA",trait.summary),"VA",
  #                 ifelse(grepl("VM",trait.summary),"VM","VD"))
    
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

# pdf('Figure Model I.pdf', width = 6, height = 6)

par(mfrow = c(3,3), mar=c(4,0,2,1),oma=c(1,3,2,1))

head_lims <- c(0,100)
tail_lims <- c(0,600)
covar_lims <- c(-60,120)
bar.width <- 0.5

# a)
plot.posterior(x=p$additive.vcov[1,1,],
               bar.width=bar.width,
               xlims=head_lims,
               cols=vc.color[1:2,],
               trait.plot="Head.G",
               trait.summary="Head.G")
mtext("a)",side=3,adj=0)
mtext("Trunk length",side=3,adj=0.5,line=1)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,800,by=10))
mtext("Additive sire variance, G",side=1,line=2,cex=0.8)

# b) 
plot.posterior(x=p$additive.vcov[2,2,],
               bar.width=bar.width,
               xlims=tail_lims,
               cols=vc.color[1:2,],
               trait.plot="Tail.G",
               trait.summary="Tail.G")
mtext("b)",side=3,adj=0)
mtext("Tail length",side=3,adj=0.5,line=1)
axis(side=1,at=seq(0,800,by=50))
mtext("Additive sire variance, G",side=1,line=2,cex=0.8)
# c) 
plot.posterior(x=p$additive.vcov[1,2,],
               bar.width=bar.width,
               xlims=covar_lims,
               cols=vc.color[1:2,],
               trait.plot="Head.Tail.cG",
               trait.summary="Head.Tail.cG")
mtext("c) ",side=3,adj=0)
mtext("Trunk-Tail",side=3,adj=0.5,line=1)
axis(side=1,at=seq(-300,300,by=20))
mtext("Additive sire covariance",side=1,line=2,cex=0.8)

# d)
plot.posterior(x=p$maternal.vcov[1,1,],
               bar.width=bar.width,
               xlims=head_lims,
               cols=vc.color[3:4,],
               trait.plot="Head.M",
               trait.summary="Head.M")
mtext("d)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,200,by=10))
mtext("Maternal effect variance, M",side=1,line=2,cex=0.8)
# e) 
plot.posterior(x=p$maternal.vcov[2,2,],
               bar.width=bar.width,
               xlims=tail_lims,
               cols=vc.color[3:4,],
               trait.plot="Tail.M",
               trait.summary="Tail.M")
mtext("e)",side=3,adj=0)
axis(side=1,at=seq(0,800,by=50))
mtext("Maternal effect variance, M",side=1,line=2,cex=0.8)
# f) 
plot.posterior(x=p$maternal.vcov[1,2,],
               bar.width=bar.width,
               xlims=covar_lims,
               cols=vc.color[3:4,],
               trait.plot="Head.Tail.cM",
               trait.summary="Head.Tail.cM")
mtext("f)",side=3,adj=0)
axis(side=1,at=seq(-300,300,by=20))
mtext("Maternal effect covariance",side=1,line=2,cex=0.8)

# g)
plot.posterior(x=p$interaction.vcov[1,1,],
               bar.width=bar.width,
               xlims=head_lims,
               cols=vc.color[5:6,],
               trait.plot="Head.D",
               trait.summary="Head.D")
mtext("g)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,800,by=10))
mtext("Interaction variance, D",side=1,line=2,cex=0.8)
# h) 
plot.posterior(x=p$interaction.vcov[2,2,],
               bar.width=bar.width,
               xlims=tail_lims,
               cols=vc.color[5:6,],
               trait.plot="Tail.D",
               trait.summary="Tail.D")
mtext("h)",side=3,adj=0)
axis(side=1,at=seq(0,800,by=50))
mtext("Interaction variance, D",side=1,line=2,cex=0.8)
# i) 
plot.posterior(x=p$interaction.vcov[1,2,],
               bar.width=bar.width,
               xlims=covar_lims,
               cols=vc.color[5:6,],
               trait.plot="Head.Tail.cD",
               trait.summary="Head.Tail.cD")
mtext("i)",side=3,adj=0)
axis(side=1,at=seq(-300,300,by=20))
mtext("Interaction covariance",side=1,line=2,cex=0.8)
############################

# dev.off()