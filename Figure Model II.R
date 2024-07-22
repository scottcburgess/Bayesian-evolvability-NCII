library('dplyr')

# Load data
load("ratio_posterior_20230117_0903.rdata") 

# Function to summarize posteriors
vecSmry <- function(x) {      
  setNames(
    c(mean(x),median(x), modeest::mlv(x, method = "venter"), HDInterval::hdi(x)),
    c("mean", "median", "mode", "lower.hdi", "upper.hdi")
  )[c("lower.hdi", "mean", "median", "mode", "upper.hdi")]
}

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


####### Function to make plot of posterior
plot.posterior <- function(x,bar.width,xlims,cols,trait.plot,trait.summary){
  # x <- p$additive # data to plot
  # x <- qgparams.post$var.a.obs # data to plot
  # bar.width = 0.00001 # width of the bars in the histogram
  # cols = "dodgerblue" # color of the bars to plot
  # trait.plot = "VA" # What trait to plot from summary_dat
  # trait.summary = "VA" # What trait to summarize in the legend
  # xlims = c(0,0.003)
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
  
  # title <- ifelse(grepl("A",trait.summary),"VA",
  #                 ifelse(grepl("D",trait),"VD","M"))

  title <- trait.summary
  
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