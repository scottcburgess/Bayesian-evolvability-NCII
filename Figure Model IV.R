library(tidyverse)
library(gridExtra)

# Load data
load('predict_settling_posterior_20230117_0038.rdata')

pr.breaks <- seq(0, 1, length.out = 100)

# Use line 10 or 11 just for viewing, but comment out when saving final pdf.
# windows(width=4,height=4) # if on PC
quartz(width=6,height=6) # if on Mac

a.b <- sapply(dimnames(p$interaction.mean)[[2]], function(m) {
  int.breaks <- seq(
    floor(min(p$interaction.mean[, m, ])), 
    ceiling(max(p$interaction.mean[, m, ])),
    length.out = length(pr.breaks)
  )
  
  pr.settle <- parallel::mclapply(int.breaks, function(x) {
    sapply(1:model.data$n.interactions, function(i) {
      p$intercept[m, ] + 
        (p$int.beta[m, ] * x) +
        (p$sire.beta[m, ] * t(p$additive.sire.eff[model.data$sire[i], m, ])) +
        (p$maternal.beta[m, ] * t(p$maternal.eff[model.data$dam[i], m, ])) +
        (p$block.beta[m, ] * t(p$block.eff[model.data$block[i], m, ]))
    }) %>%
      as.vector() %>%
      plogis() %>%
      cut(pr.breaks, include.lowest = TRUE) %>%
      table()
  }, mc.cores = 10) %>% 
    do.call(cbind, .) %>% 
    as.data.frame() %>% 
    remove_rownames() %>% 
    setNames(int.breaks) %>% 
    mutate(pr.settle = apply(cbind(pr.breaks[-length(pr.breaks)], pr.breaks[-1]), 1, mean)) %>% 
    pivot_longer(-pr.settle, names_to = 'interaction.mean', values_to = 'freq') %>% 
    mutate(
      interaction.mean = as.numeric(interaction.mean),
      interaction.mean.lik = dnorm(interaction.mean, mean(p$interaction.mean[, m, ]), sd(p$interaction.mean[, m, ])),
      wt = freq * interaction.mean.lik
    )
  
  ggplot(pr.settle, aes(interaction.mean, pr.settle)) +
    geom_tile(aes(fill = wt)) +
    geom_hline(yintercept = 0.5, color = 'white', alpha = 0.6, linetype = 'dashed', linewidth = 1) +
    scale_fill_viridis_c(option = 'inferno') +
    annotate(
      'label', 
      x = -Inf, 
      y = -Inf, 
      label = paste(
        paste0('median = ', round(median(p$int.beta[m, ]), 3), '\n'),
        paste0('mode = ', round(modeest::venter(p$int.beta[m, ]), 3), '\n'),
        paste0('95% = ', paste(round(HDInterval::hdi(p$int.beta[m, ]), 3), collapse = ' - ')),
        sep = '',
        collapse = ''
      ),
      color="white",
      fill=NA,
      label.size=NA,
      size=3,
      hjust = 0,
      vjust = 0
    ) +
    labs(
      x = paste(ifelse(m == 'head', 'Head', 'Tail'), 'length (mean per full-sib family)'), 
      y = 'Probability of settling',
      title = ifelse(m == 'head', 'a)', 'b)') 
    ) +
    coord_cartesian(xlim = range(int.breaks), ylim = c(0, 1), expand = FALSE) +
    theme_minimal() +
    theme(
      legend.position = 'none',
      panel.grid = element_blank()
    )
}, simplify = FALSE)



effect.colors <- c('Additive Sire' = '#0081a7', 'Maternal' = '#f07167')
c.d <- sapply(dimnames(p$interaction.mean)[[2]], function(m) {
  smry <- lapply(p[c('sire.beta', 'maternal.beta')], function(x) {
    paste(
      paste0('median = ', round(median(x[m, ]), 3), '\n'),
      paste0('mode = ', round(modeest::venter(x[m, ]), 3), '\n'),
      paste0('95% = ', paste(round(HDInterval::hdi(x[m, ]), 3), collapse = ' - ')),
      sep = '',
      collapse = ''
    )
  })
  
  pred.pr.settle |> 
    filter(metric == m & effect.type %in% c('Additive Sire', 'Maternal')) |>  
    ggplot(aes(x = effect)) +
    geom_hline(yintercept = mean(model.data$settle[, 1]), linetype = 'dashed') +
    # geom_ribbon(aes(ymin = lower.hdi, ymax = upper.hdi, fill = effect.type), alpha = 0.5) +
    geom_line(aes(y = median, color = effect.type), linewidth = 2) +  
    scale_color_manual(values = effect.colors) +
    # scale_fill_manual(values = effect.colors) +
    annotate(
      'label', 
      x = -Inf, 
      y = -Inf, 
      label = paste('Additive Sire\n', smry[[1]], sep = '', collapse = ''),
      fill=NA,
      label.size=NA,
      size=3,
      hjust = 0,
      vjust = 0,
      color = effect.colors['Additive Sire']
    ) +
    annotate(
      'label', 
      x = Inf, 
      y = -Inf, 
      label = paste('Maternal\n', smry[[2]], sep = '', collapse = ''),
      fill=NA,
      label.size=NA,
      size=3,
      hjust = 1,
      vjust = 0,
      color = effect.colors['Maternal']
    ) +
    coord_cartesian(xlim = c(-15, 15), ylim = c(0, 1), expand = FALSE) +
    labs(
      x = paste0('Marginal effect (', m, ' length)'), 
      y = 'Probability of settling',
      title = ifelse(m == 'head', 'c)', 'd)')
    ) +
    theme(
      panel.grid = element_blank(),
      legend.position = 'none'
    )
}, simplify = FALSE)


# Uncomment lines 126 and 128 to send figure to a PDF (best for copy/paste inclusion in manuscript)
# pdf('Figure Model IV.pdf', width = 10, height = 10)
do.call(grid.arrange, c(a.b, c.d, ncol = 2))
# dev.off()
