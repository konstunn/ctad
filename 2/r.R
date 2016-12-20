#!/usr/bin/Rscript

n = 10^3

x1 = rnorm(n)

library(fitdistrplus)

fitdist(x1, "mle", distr="norm")

fitdist(x1, "mge", distr="norm", gof="KS") 

fitdist(x1, "mge", distr="norm", gof="CvM")

fitdist(x1, "mge", distr="norm", gof="AD")


