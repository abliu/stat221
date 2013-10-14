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