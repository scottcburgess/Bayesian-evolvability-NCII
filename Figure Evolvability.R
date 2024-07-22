library('dplyr')

# Load data
load("head_tail_posterior_20230116_1944.rdata") 
pE_head_tail <- p$E
e.params_BetaMCMC_head_tail <- e.params_BetaMCMC

load("hatch_settle_posterior_20230119_1741.rdata") 
pE_hatch_settle <- p$E
e.params_BetaMCMC_hatch_settle <- e.params_BetaMCMC

# Get the posterior samples by averaging across all of the random beta's
eB_posterior_head_tail <- apply(e.params_BetaMCMC_head_tail$post.dist$eB, 1, mean)
rB_posterior_head_tail <- apply(e.params_BetaMCMC_head_tail$post.dist$rB, 1, mean)
cB_posterior_head_tail <- apply(e.params_BetaMCMC_head_tail$post.dist$cB, 1, mean)
eB_posterior_hatch_settle <- apply(e.params_BetaMCMC_hatch_settle$post.dist$eB, 1, mean)
rB_posterior_hatch_settle <- apply(e.params_BetaMCMC_hatch_settle$post.dist$rB, 1, mean)
cB_posterior_hatch_settle <- apply(e.params_BetaMCMC_hatch_settle$post.dist$cB, 1, mean)

# Function to summarize vectors
vecSmry <- function(x) {      
  setNames(
    c(median(x), modeest::mlv(x, method = "venter"), HDInterval::hdi(x)),
    c("median", "mode", "lower.hdi", "upper.hdi")
  )[c("lower.hdi", "median", "mode", "upper.hdi")]
}

####### Function to make plot of posterior
plot.posterior <- function(x,bar.width,xlims,cols,trait.plot){
  # For testing
  # x <- eB_posterior_head_tail # data to plot
  # bar.width = 0.0001 # width of the bars in the histogram
  # xlims <- c(0,0.001)
  # cols = "dodgerblue" # color of the bars to plot
  # trait.plot = "e_mean" # What trait to plot from summary_dat
  summary.matrix = vecSmry(x)
  blims <- max(abs(x))*2
  brks <- seq(-blims,blims,by=bar.width)
  foo <- hist(x, breaks=brks,plot=F)
  d <- data.frame(x=foo$mids, y=foo$density)
  d <- d %>% filter(y>0)
  ylims <- c(0,max(d$y))

  plot(xlims,ylims,type="n",bty="l",ylab="",xlab="",axes=F)

  interval <- as.numeric(round(summary.matrix[c("lower.hdi","upper.hdi")],6))
  
  with(d,segments(x,
                  rep(0,length(x)),
                  x,
                  y,
                  lend=2,
                  col=ifelse(d$x > interval[1] & d$x < interval[2],
                             cols[1],cols[2])))
  
  abline(v=0)
# Add summaries to plot
  median.text <- eval(paste("median = ",round(summary.matrix["median"],6)))
  mode.text <- eval(paste("mode = ",round(summary.matrix["mode"],6)))  
  int.text <- eval(paste("95% hdi = ",interval[1], " - ", interval[2],sep=""))  
  
  # title <- ifelse(grepl("e_",trait.plot),expression(paste("e",beta)),
  #                 ifelse(grepl("r_",trait.plot),expression(paste("r",beta)),
  #                        expression(paste("c",beta))))
  

  legend(xlims[2]*0.9,ylims[2],
         legend=rbind(median.text,mode.text,int.text),
         bty="n",
         adj=1,
         # title=title,
         # title.adj=0,
         xpd=T)
}
###################################################





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
bar.width <- 0.000001
# a)
plot.posterior(x=eB_posterior_head_tail,
               bar.width=bar.width,
               xlims=e_lims,
               cols=vc.color[1:2,],
               trait.plot = "e_mean")
mtext("a)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,1,by=0.0001))
mtext(expression(paste("Evolvability, e",beta)),side=1,line=2.5,cex=0.8)
# b) 
plot.posterior(x=rB_posterior_head_tail,
               bar.width=bar.width,
               xlims=e_lims,cols=vc.color[3:4,],
               trait.plot = "r_mean")
mtext("b)",side=3,adj=0)
mtext("Trunk-Tail",side=3,adj=0.5,line=1.5)
axis(side=1,at=seq(0,1,by=0.0001))
mtext(expression(paste("Respondability, r",beta)),side=1,line=2.5,cex=0.8)
# c) 
plot.posterior(x=cB_posterior_head_tail,
               bar.width=bar.width,
               xlims=e_lims,
               cols=vc.color[5:6,],
               trait.plot = "c_mean")
mtext("c)",side=3,adj=0)
axis(side=1,at=seq(0,1,by=0.0001))
mtext(expression(paste("Conditional evolvability, r",beta)),side=1,line=2.5,cex=0.8)

# Hatch - settle
e_lims <- c(0,2)
bar.width <- 0.005
# d)
plot.posterior(x=eB_posterior_hatch_settle,
               bar.width=bar.width,
               xlims=e_lims,
               cols=vc.color[1:2,],
               trait.plot = "e_mean")
mtext("d)",side=3,adj=0)
mtext("Probability density",side=2,line=1,cex=0.8)
axis(side=1,at=seq(0,2,by=0.1))
mtext(expression(paste("Evolvability, e",beta)),side=1,line=2.5,cex=0.8)
# e) 
plot.posterior(x=rB_posterior_hatch_settle,
               bar.width=bar.width,
               xlims=e_lims,
               cols=vc.color[3:4,],
               trait.plot = "r_mean")
mtext("e)",side=3,adj=0)
mtext("Hatch-Settle",side=3,adj=0.5,line=1.5)
axis(side=1,at=seq(0,2,by=0.1))
mtext(expression(paste("Respondability, r",beta)),side=1,line=2.5,cex=0.8)
# f) 
plot.posterior(x=cB_posterior_hatch_settle,
               bar.width=bar.width,
               xlims=e_lims,
               cols=vc.color[5:6,],
               trait.plot = "c_mean")
mtext("f)",side=3,adj=0)
axis(side=1,at=seq(0,2,by=0.1))
mtext(expression(paste("Conditional evolvability, r",beta)),side=1,line=2.5,cex=0.8)

# dev.off()
