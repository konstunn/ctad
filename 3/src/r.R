#!/usr/bin/Rscript --vanilla

N <- 16600

# sample size
n <- 100

# Generate sample
# TODO: customize according to H
# H0: sample is from Cauchy distribution
# H1: sample is from Normal distribution
X <- rcauchy(n*N) # here we set H0 true
X <- split(X , ceiling(seq_along(X)/n))

# number of intervals for grouping (number of bins)
k <- 5

# optimal probabilities vector elements step for asymptotically optimal
# grouping of sample from Cauchy distribution (we assume H0 is true)
dP <- 1/k
P <- seq(from=0, to=1, by=dP)

# TODO: for testing composite hypothesis first estimate distribution parameters 
# from x, assuming x is from Cauchy distribution (we assume H0 is true)
location=0 # Here we test
scale=1    # simple hypothesis

# calculate groups (bins) break points (i.e. quantiles)
x_i <- qcauchy(P, location=location, scale=scale)

# expected probabilities (frequencies)
expected <- rep(dP, k)

# perform sample grouping (binning)
factors <- lapply(X, cut, x_i, include.lowest=TRUE)
groups <- Map(split, X, factors)

# calculate observed probabilities (frequencies)
nestlen <- function(x) lapply(x, length)
observed <- lapply(groups, nestlen)
observed <- matrix(unlist(observed), ncol=k, byrow=TRUE)

# calculate chisq test statistic values
machisq.test <- function(o, p) chisq.test(o, p=p)$statistic
chisq <- apply(observed, 1, machisq.test, p=expected)

df <- data.frame(chisq)
library(ggplot2)
p <- ggplot() + geom_density(aes(chisq), data=df)
print(p)

