
# sample size
n <- 100

# Generate sample
# TODO: customize according to H
# H0: sample is from Cauchy distribution
# H1: sample is from Normal distribution
x <- rcauchy(n) # here we set H0 true

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
factors <- cut(x, x_i, include.lowest=TRUE)
groups <- split(x, factors)

# calculate observed probabilities (frequencies)
observed <- lapply(groups, FUN=length)
observed <- as.numeric(observed)

print(chisq.test(observed, p=expected))

# calculate chi square test statistic by hand
expected <- expected * n
print(sum((observed-expected)^2 / expected))
