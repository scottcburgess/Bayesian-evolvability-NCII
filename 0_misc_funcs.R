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
