rm(list=ls())
library('tidyverse')
library('gridExtra')

# Load data
load('3_Model_outputs/Model_IV_posterior_20230117_0038.rdata')

pr.breaks <- seq(0, 1, length.out = 100)

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
    }) |>
      as.vector() |>
      plogis() |>
      cut(pr.breaks, include.lowest = TRUE) |>
      table()
  }, mc.cores = 10) |> 
    do.call(cbind, .) |> 
    as.data.frame() |> 
    remove_rownames() |> 
    setNames(int.breaks) |> 
    mutate(pr.settle = apply(cbind(pr.breaks[-length(pr.breaks)], pr.breaks[-1]), 1, mean)) |> 
    pivot_longer(-pr.settle, names_to = 'interaction.mean', values_to = 'freq') |> 
    mutate(
      interaction.mean = as.numeric(interaction.mean),
      interaction.mean.lik = dnorm(interaction.mean, mean(p$interaction.mean[, m, ]), sd(p$interaction.mean[, m, ])),
      wt = freq * interaction.mean.lik
    )
  
  ggplot(pr.settle, aes(interaction.mean, pr.settle)) +
    geom_tile(aes(fill = wt)) +
    geom_hline(yintercept = 0.5, color = 'white', alpha = 0.6, linetype = 'dashed', linewidth = 0.7) +
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
      size=2,
      hjust = 0,
      vjust = 0
    ) +
    labs(
      x = paste(ifelse(m == 'head', 'Trunk', 'Tail'), 'length\n(mean per full-sib family)'), 
      y = 'Probability of settling',
      title = ifelse(m == 'head', 'a)', 'b)') 
    ) +
    coord_cartesian(xlim = range(int.breaks), ylim = c(0, 1), expand = FALSE) +
    theme_minimal() +
    theme(
      legend.position = 'none',
      panel.grid = element_blank(),
      # axis.title.x = element_text(hjust = 0.5),
      # axis.title.y = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 10, face = "plain")
    )
}, simplify = FALSE)


pdf('5_Figure_outputs/Figure 5.pdf', width = 5, height = 2.5)
do.call(grid.arrange, c(a.b, ncol = 2, nrow = 1))
dev.off()
