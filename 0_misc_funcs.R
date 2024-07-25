vecSmry <- function(x) {      
  setNames(
    c(median(x), modeest::mlv(x, method = "venter"), HDInterval::hdi(x)),
    c("median", "mode", "lower.hdi", "upper.hdi")
  )[c("lower.hdi", "median", "mode", "upper.hdi")]
}

comp.smry <- function(x, gt = 0) {
  library(tidyverse)
  df <- if(is.null(dim(x))) {
    as.data.frame(rbind(vecSmry(x)))
  } else if(length(dim(x)) == 2) {
    as.data.frame(t(apply(x, 1, vecSmry)))
  } else {
    metrics <- dimnames(x)[[1]]
    t(sapply(c(lapply(metrics, rep, times = 2), list(metrics)), function(i) {
      xi <- x[i[1], i[2], ]
      result <- vecSmry(xi)
      setNames(
        c(result, mean(xi > gt)), 
        c(names(result), paste0("pct.gt.", gt))
      )
    })) %>% 
      as.data.frame() %>% 
      mutate(measure = c(metrics, "cov")) %>% 
      column_to_rownames("measure")
  }
  pander::pander(df, split.tables = Inf, keep.line.breaks = TRUE)
}


plot.metric <- function(x) {
  library(tidyverse)
  if(is.null(dim(x))) {
    x %>% 
      enframe() %>%
      ggplot(aes(value)) +
      geom_histogram(bins = 100) +
      labs(x = "Value", y = "Count")
  } else {
    metrics <- dimnames(x)[[1]]
    sapply(c(lapply(metrics, rep, times = 2), list(metrics)), function(i) {
      x[i[1], i[2], ]
    }) %>% 
      as.data.frame() %>% 
      setNames(c(metrics, "cov")) %>% 
      pivot_longer(everything()) %>% 
      mutate(name = factor(name, levels = c(metrics, "cov"))) %>% 
      ggplot(aes(value)) +
      geom_histogram(bins = 100) +
      labs(x = "Value", y = "Count") +
      facet_wrap(~ name, scales = "free")
  }
}


# Function to make plot of posterior
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