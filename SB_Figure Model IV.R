library('dplyr')
library('HDInterval') # Used for hdi()
library('modeest') # Used for mlv()
library('gridExtra')

# Load data
load("predict_settling_posterior_20230117_0038.rdata")


# Function to summarize posteriors
vecSmry <- function(x) {      
  setNames(
    c(mean(x),median(x), modeest::mlv(x, method = "venter"), HDInterval::hdi(x)),
    c("mean", "median", "mode", "lower.hdi", "upper.hdi")
  )[c("lower.hdi", "mean", "median", "mode", "upper.hdi")]
}

# Summarize the posteriors for the slope (betas) 
# to add to plots
summary_beta <- as.data.frame(rbind(vecSmry(p$int.beta[1,]), # head
                                    vecSmry(p$int.beta[2,]), # tail
                                    vecSmry(p$sire.beta[1,]), # head sire
                                    vecSmry(p$sire.beta[2,]), # tail sire
                                    vecSmry(p$maternal.beta[1,]), # head maternal
                                    vecSmry(p$maternal.beta[2,])), # tail maternal
                              row.names=c("head","tail",
                                          "sire.beta.head","sire.beta.tail",
                                          "maternal.beta.head","maternal.beta.tail")) 


###### Plot posterior of head or tail length vs probability of settling
# To-do:
# 1) Add legend.text to the bottom left of each plot
# 2) Capitilize "head" and "tail" label on the x-axis

n.pts <- 1000 # use 10 for a quick look, use 1000 for the final version
pr.breaks <- seq(0, 1, length.out = n.pts)

for(m in dimnames(p$interaction.mean)[[2]]) {
  int.breaks <- seq(
    floor(min(p$interaction.mean[, m, ])), 
    ceiling(max(p$interaction.mean[, m, ])),
    length.out = n.pts
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
    mutate(pr.settle = apply(cbind(pr.breaks[-n.pts], pr.breaks[-1]), 1, mean)) %>% 
    pivot_longer(-pr.settle, names_to = "interaction.mean", values_to = "freq") %>% 
    mutate(
      interaction.mean = as.numeric(interaction.mean),
      interaction.mean.lik = dnorm(interaction.mean, mean(p$interaction.mean[, m, ]), sd(p$interaction.mean[, m, ])),
      wt = freq * interaction.mean.lik
    )
  
  # Add summaries to plot
  median.text <- eval(paste("median = ",summary_beta[m, which(names(summary_beta)=="median")]))  
  mode.text <- eval(paste("mode = ",summary_beta[m, which(names(summary_beta)=="mode")]))  
  int <- summary_beta[m, which(names(summary_beta) %in% c("lower.hdi","upper.hdi"))]
  int.text <- eval(paste("95% hdi = ",int[1], " - ", int[2],sep=""))  
  
  xpos <- 0.1
  ypos <- min(int.breaks)+20
  legend.text <- rbind(median.text,mode.text,int.text)
  ##
  
  gg <- ggplot(pr.settle, aes(interaction.mean, pr.settle)) +
    geom_tile(aes(fill = wt)) +
    geom_hline(yintercept = 0.5, color = "white", alpha = 0.6, linetype = "dashed", linewidth = 1) +
    scale_fill_viridis_c(option = "inferno") +
    labs(x = eval(paste(m, "length (mean per full-sib family)")), 
         y = "Probability of settling",
         title=ifelse(m=="head","a) Head length","b) Tail length")) +
    coord_cartesian(xlim = range(int.breaks), ylim = c(0, 1), expand = FALSE) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    )

    if(m=="head") gg_head <- gg
  if(m=="tail") gg_tail <- gg
}


## Plot Marginal Effects
## To do:
# 1) add the lower.hdi and upper.hdi interval as a shaded area.
# 2) only need one legend for c) and d), and place it inside the plot area
# 3) Add the summary text for the betas

gg_marginal_head <- pred.pr.settle %>% 
  filter(metric=="head",effect.type %in% c("Additive Sire","Maternal")) %>% 
  ggplot(aes(x = effect)) + xlim(-10,10) + ylim(0.4,0.8) +
  geom_hline(yintercept = mean(model.data$settle[, 1]), linetype = "dashed") +
  geom_line(aes(y = median, color = effect.type)) +
  labs(x = "Marginal effect", y = "Probability of settling",
       title="c) Head length") +
  theme(panel.grid = element_blank(),
        legend.position = "right")

gg_marginal_tail <- pred.pr.settle %>% 
  filter(metric=="tail",effect.type %in% c("Additive Sire","Maternal")) %>% 
  ggplot(aes(x = effect)) + xlim(-10,10) + ylim(0.4,0.8) +
  geom_hline(yintercept = mean(model.data$settle[, 1]), linetype = "dashed") +
  geom_line(aes(y = median, color = effect.type)) +
  labs(x = "Marginal effect", y = "Probability of settling",
       title="d) Tail length") +
  theme(panel.grid = element_blank(),
        legend.position = "right")

# Make the actual plot
# windows(width=4,height=4) # PC
quartz(width=6,height=6) # Mac
par(oma=c(2,2,2,2))
grid.arrange(gg_head, 
             gg_tail, 
             gg_marginal_head,
             gg_marginal_tail,
             ncol=2)
