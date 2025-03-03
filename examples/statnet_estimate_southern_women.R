#!/usr/bin/env Rscript
##
## File:    statnet_estimate_southern_women.R
## Author:  Alex Stivala
## Created: December 2024
##
## Estimate models for the Davis "Southern Women" two-mode network
## using statnet and the new b1np4c/b2np4c statistics in
## ergm.terms.contrib.
##
## Usage: Rscript statnet_estimate_southern_women.R
##
## Output files in cwd (WARNING overwrites if present):
##
##   southern_women_models.tex
##   southern_women_model1_mcmcdiag.eps
##   southern_women_mddel2_mcmcdiag.eps
##   southern_women_mddel3_mcmcdiag.eps
##   southern_women_mddel4_mcmcdiag.eps
##   southern_women_model1_gof.eps
##   southern_women_mddel2_gof.eps
##   southern_women_mddel3_gof.eps
##   southern_women_mddel4_gof.eps
##   southern_women_model1_cycledist.eps
##   southern_women_mddel2_cycledist.eps
##   southern_women_mddel3_cycledist.eps
##   southern_women_mddel4_cycledist.eps
##
##
## Runtime approx. 2 minutes (8 cores on Lenovo Legion PC, MAX_CYCLELEN = 8)
##         approx. 4 minutes (8 cores on Lenovo Legion PC, MAX_CYCLELEN = 10)
##         cancelled after approx. 45 minutes (8 cores on Lenovo Legion PC,
##         MAX_CYCLELEN = 12), still on model 1 (note models 2 to 4 only
##         took 6 minutes)
## (most time is spent in plot_cycledist() which is the only
## part using parallel mclapply(); reducing MAX_CYCLELEN makes it
## much faster). Set environment variable MC_CORES to choose
## number of cores used.
##

cat("Start:", date(), '\n')

library(statnet)
library(latentnet)
library(ergm.terms.contrib)
library(texreg)
library(parallel)
library(ggplot2)
library(scales)

data(davis) # Southern Women data from latentnet package

setEPS()  # postscript() will use EPS settings



##############################################################################
###
### Functions
###
##############################################################################

obscolour <- 'red' # colour to plot observed graph points/lines
## simulated graph statistics will be boxplot on same plot in default colour

## Using theme_classic() to get no grey background and no gridlines
## as required by some journals e.g. J. Complex Networks
ptheme <- theme_classic() +  theme(legend.position = 'none')

# http://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
orig_scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}
my_scientific_10 <- function(x) {
# also remove + and leading 0 in exponennt
  parse( text=gsub("e", " %*% 10^", gsub("e[+]0", "e", scientific_format()(x))) )
   
}

##
## Plot boxplot of cycle length distributions from networks simulatd
## from model with values from observed network net_obs on smae plot.
## Adapted from EstimNetDirected/scripts/simFitPlots.R
##
plot_cycledist <- function(net_obs, model, epsfilename) {
  
  MAX_CYCLELEN <- 10
  num_sim <- 100

  ## for library(parallel) using simplify2array(mclapply(...)) instead of
  ## sapply(...) and mclapply(...) instead of lapply(...)
  cat('mc.cores =', getOption("mc.cores"), '\n')

  
  ##
  ## Generate simulated graphs from model
  ##
  cat('simulating networks from model...')
  print(system.time(sim_networks <- simulate(model, nsim = num_sim)))
  
  ##
  ## Cycle length distribution (up to MAX_CYCLLEN only, as then gets 
  ## extremely slow and huge numbers of cycles).
  ##
  cat('cycle length distribution, MAX_CYCLELEN = ', MAX_CYCLELEN, '\n')
  if (is.bipartite(net_obs)) {
    cyclelens <- seq(4, MAX_CYCLELEN, 2) #only even length cycles in bipartite
  } else {
    cyclelens <- seq(3, MAX_CYCLELEN)
  }
  cat('computing cycle length distribution in observed graph...')
  print(system.time(obs_cycledist <- simplify2array(mclapply(cyclelens,
                                                             function(x) summary(net_obs ~ cycle(x))))))
  cat('obs_cycledist = ', obs_cycledist, '\n')
  cat('computing cycle length distribution in simulated graphs...')
  print(system.time(sim_cycledist <- mclapply(sim_networks,
                                              function(g) summary(g ~ cycle(cyclelens)))))
  ##print(sim_cycledist)
  cyclelen_df <- data.frame(sim = rep(1:num_sim, each = length(cyclelens)),
                            cyclelen = rep(cyclelens, num_sim),
                            count = NA)
  for (i in 1:num_sim) {
    cyclelen_df[which(cyclelen_df[,"sim"] == i), "count"] <- sim_cycledist[[i]]
  }
  obs_cyclelen_df <- data.frame(cyclelen = cyclelens,
                                count = obs_cycledist)
  print(cyclelen_df)#XXX
  cyclelen_df$cyclelen <- factor(cyclelen_df$cyclelen)
  obs_cyclelen_df$cyclelen <- factor(obs_cyclelen_df$cyclelen)
  p <- ggplot(cyclelen_df, aes(x = cyclelen, y = count)) + geom_boxplot()
  p <- p + geom_line(data = obs_cyclelen_df, aes(x = cyclelen, y = count,
                                                 colour = obscolour, group = 1))
  p <- p + ptheme + xlab('cycle length') + ylab('count')
  p <- p + guides(x = guide_axis(check.overlap = TRUE))
  p <- p + scale_y_log10() + ylab("count (log scale)")
  ## write to separate file (add points for separate plot only,
  p <- p + geom_point(data = obs_cyclelen_df, aes(x = cyclelen, y = count,
                                                    colour = obscolour, group = 1))
  p <- p + scale_y_log10(labels = my_scientific_10)
  ## increase font size to make it readable when included as subfigure
  ## and reduced to smaller panels in LaTeX
  p <- p + theme_classic(base_size = 32) + theme(legend.position = 'none')
  cycledist_outfilename <- epsfilename
  cat("writing cycle length distribution plot to EPS file ", cycledist_outfilename, "\n")
  postscript(cycledist_outfilename, horizontal = FALSE, onefile = FALSE,
             paper = "special", width = 9, height = 6)
  print(p)
  dev.off()  
}



##############################################################################
###
### Main
###
##############################################################################

print(R.version)
print(Sys.info())
print(sessionInfo())

print(davis)


### Model similar to Wang et al. (2009) Model (8.3)
system.time(davis_model1 <- ergm(davis ~ edges + b1star(2) + b2star(2),  control = control.ergm(main.method = "Stochastic-Approximation")))
print(summary(davis_model1))
postscript('southern_women_model1_mcmcdiag.eps')
mcmc.diagnostics(davis_model1, vars.per.page=6)
dev.off()

system.time( davis_model1_gof <- gof(davis_model1, GOF = ~ b1degree + b2degree+dspartners + distance + model) )
print(davis_model1_gof)
postscript('southern_women_model1_gof.eps')
par(mfrow=c(3, 2))
plot(davis_model1_gof)
dev.off()

plot_cycledist(davis, davis_model1, 'southern_women_model1_cycledist.eps')



## Converged model in only a few seconds, OK GoF:
system.time(davis_model2 <- ergm(davis ~ edges + gwb1degree(1, TRUE) + gwb2degree(1, TRUE),  control = control.ergm(main.method = "Stochastic-Approximation")))
print(summary(davis_model2))
postscript('southern_women_model2_mcmcdiag.eps')
mcmc.diagnostics(davis_model2, vars.per.page=6)
dev.off()

system.time( davis_model2_gof <- gof(davis_model2, GOF = ~ b1degree + b2degree+dspartners + distance + model) )
print(davis_model2_gof)
postscript('southern_women_model2_gof.eps')
par(mfrow=c(3, 2))
plot(davis_model2_gof)
dev.off()

plot_cycledist(davis, davis_model2, 'southern_women_model2_cycledist.eps')


#### Model similar to Wang et al. (2009) Model (8.6)
## Does not converge (Unconstrained MCMC sampling did not mix at all. Optimization cannot continue):
##system.time( davis_model3 <- ergm(davis ~ edges + gwb1degree(1, TRUE) + gwb2degree(1, TRUE) + gwb1dsp(log(2), TRUE) + gwb2dsp(log(2), TRUE)) )
##system.time( davis_model3 <- ergm(davis ~ edges + gwb1degree(1, TRUE) + gwb2degree(1, TRUE) + gwb1dsp(0, TRUE) + gwb2dsp(0, TRUE)) )
##system.time( davis_model3 <- ergm(davis ~ edges + gwb1degree(1, TRUE) + gwb2degree(1, TRUE) + gwb1dsp(1, TRUE) + gwb2dsp(1, TRUE)) )
#### Note only converges with stochastic approximation:
system.time( davis_model3 <- ergm(davis ~ edges + gwb1degree(1, TRUE) + gwb2degree(1, TRUE) + gwb1dsp(0.5, TRUE) + gwb2dsp(0.5, TRUE), control = control.ergm(main.method = "Stochastic-Approximation")) )
print(summary(davis_model3))
postscript('southern_women_model3_mcmcdiag.eps')
mcmc.diagnostics(davis_model3, vars.per.page=6)
dev.off()

system.time( davis_model3_gof <- gof(davis_model3, GOF = ~ b1degree + b2degree+dspartners + distance + model) )
print(davis_model3_gof)
postscript('southern_women_model3_gof.eps')
par(mfrow=c(3, 2))
plot(davis_model3_gof)
dev.off()

plot_cycledist(davis, davis_model3, 'southern_women_model3_cycledist.eps')


## Converges (only a minute or so). OK GoF:
system.time( davis_model4 <- ergm(davis ~ edges + gwb1degree(1, TRUE) + gwb2degree(1, TRUE) + b2np4c(1/5), control = control.ergm(main.method = "Stochastic-Approximation")) )
print(summary(davis_model4))
postscript('southern_women_model4_mcmcdiag.eps')
mcmc.diagnostics(davis_model4, vars.per.page=6)
dev.off()

system.time( davis_model4_gof <- gof(davis_model4, GOF = ~ b1degree + b2degree+dspartners + distance + model) )
print(davis_model4_gof)
postscript('southern_women_model4_gof.eps')
par(mfrow=c(3, 2))
plot(davis_model4_gof)
dev.off()

plot_cycledist(davis, davis_model4, 'southern_women_model4_cycledist.eps')


## Output models to LaTeX table
texreg(list(davis_model1, davis_model2, davis_model3, davis_model4), file = 'southern_women_models.tex', table=FALSE)


cat("End:", date())
