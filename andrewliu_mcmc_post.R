source("andrewliu_mcmc.R")

# 1.4: MCMC Visualization
library(MASS)
outs <- list.dirs(path = ".")
outs = outs[substr(outs, 3, 6) == 'out-']
for(out in outs) {
  load(file.path(out, "output.rda"))
  pdf(paste(out, "pdf", sep="."))
  f <- kde2d(results[,1], results[,2], n = 100, h = 2*c(width.SJ(results[,1], method="dpi"), width.SJ(results[,2])), lims=c(50,400,range(results[,2])))
  image(f)
  contour(f, add=TRUE, lwd=1)
  dev.off()
}

# 1.5: Normalizing Constants
impala_Z = integrate(function(x) exp(margN(x, Impala)), 0, Inf)[[1]] # 2.34e-76
waterbucks_Z = sum(exp(margN(0:100000, Waterbucks))) # 1.100903e-223

# 1.6: N > 100
# Analytic
impala_100 = exp(log(integrate(function(x) exp(margN(x, Impala)), 100, Inf)[[1]]) - log(impala_Z)) # 0.359171
waterbucks_100 = exp(log(sum(exp(margN(100:100000, Waterbucks)))) - log(waterbucks_Z)) # 0.9607001

# MCMC
prob_100 = c()
means = c()
stds = c()
for(out in outs) {
  load(file.path(out, "output.rda"))
  prob_100 = c(prob_100, sum(results[,1] > 100) / nrow(results))
  means = c(means, mean(results[,1]))
  stds = c(stds, sd(results[,1]))
}
# Impala has values < 50; Waterbucks has >= 50.
summary(prob_100[prob_100 < 0.5])
summary(prob_100[prob_100 > 0.5])