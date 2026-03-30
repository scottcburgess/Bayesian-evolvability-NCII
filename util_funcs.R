addQGmetrics <- function(p) {
  p$VA <- 4 * p$sire.vcov
  p$VM <- p$dam.vcov - p$sire.vcov
  p$VD <- 4 * p$interaction.vcov
  p$VP <- p$VA + p$VM + p$VD
  p$H <- p$VA / p$VP
  if(!is.null(p$resid.vcov)) {
    p$VP <- p$VP + p$resid.vcov
    p$VR <- p$resid.vcov
  }
  p
}


# Use QGglmm to extract full variance/covariance matrix on observed scale
convertVCVscale.II <- function(metric, p) {
  vcv <- parallel::mclapply(1:dim(p[[metric]])[3], function(i) {
    QGglmm::QGmvparams(
      vcv.G = p[[metric]][, , i],
      vcv.P = p$VP[, , i],
      predict = qlogis(p$pr.block[, , i]),
      models = c('binom1.logit', 'binom1.logit'),
      verbose = FALSE
    )
  }, mc.cores = 10) |> 
    purrr::list_transpose()
  
  sapply(vcv, function(x) {
    if(is.null(dim(x[[1]]))) {
      x <- do.call(rbind, x)
      colnames(x) <- dimnames(p[[metric]])[[1]]
      x
    } else {
      do.call(
        abind::abind, 
        c(x, list(along = 3, new.names = dimnames(p[[metric]])))
      )
    }
  })
}

# Use QGglmm to extract full variance/covariance matrix on observed scale
convertVCVscale.III <- function(p) {
  # Run QGmvparams across iterations
  vcv <- parallel::mclapply(
    X = seq_len(dim(p$VA)[3]),
    FUN = function(i) {
      QGglmm::QGmvparams(
        vcv.G   = p$VA[, , i],
        vcv.P   = p$VP[, , i],
        predict = p$overall.block.mean[, , i],
        models  = c("Gaussian", "Gaussian", "binom1.logit"),
        verbose = FALSE
      )
    },
    mc.cores = 10
  ) |> purrr::list_transpose()
  
  # Collapse results into matrices or arrays
  lapply(vcv, function(x) {
    if(is.null(dim(x[[1]]))) {
      out <- do.call(rbind, x)
      dimnames(out) <- list(
        iter  = dimnames(p$VA)[[3]],
        trait = dimnames(p$VA)[[1]]
      )
    } else {
      out <- abind::abind(x, along = 3)
      dimnames(out)[[3]] <- dimnames(p$VA)[[3]]
    }
    out
  })
}

smrzPost <- function(post, v) {
  post.smry <- summary(
    post,
    vars = c(
      'deviance', 'sire.vcov', 'dam.vcov', 'interaction.vcov',
      'pr.overall', 'pr.block'
    ) 
  ) |>  
    as.data.frame() |> 
    rownames_to_column('metric') |>
    select(metric, SSeff:psrf) |> 
    pivot_longer(-metric, names_to = 'diag', values_to = 'values') 
  
  diag.smry <- post.smry |> 
    group_by(diag) |> 
    summarize(
      median = median(values),
      lower = unname(quantile(values, 0.025)),
      upper = unname(quantile(values, 0.975)),
      .groups = 'drop'
    )
  
  list(post = post.smry, diag = diag.smry)
}

smrzPPC <- function(ppc) {
  ppc |> 
    group_by(metric) |> 
    summarize(
      median.pct = median(pct.gte.obs),
      lower.pct = unname(quantile(pct.gte.obs, 0.025)),
      upper.pct = unname(quantile(pct.gte.obs, 0.975)),   
      median.diff = median(mean.diff),
      lower.diff = unname(quantile(mean.diff, 0.025)),
      upper.diff = unname(quantile(mean.diff, 0.975)),
      .groups = 'drop'
    )
}