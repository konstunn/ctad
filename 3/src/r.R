
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
opt_P <- 1/k

# calculate intervals break points
x_i <- quantile(x, seq(from=0, to=1, by=opt_P))

# TODO: if hypothesis should be composite, then estimate these parameters first
p <- pcauchy(x_i, location=0, scale=1)

# expected probabilities (frequencies)
expected <- zoo::rollapply(p, width=2, function(x) x[2] - x[1])

# perform sample grouping (binning)
factors <- cut(x, x_i, include.lowest=TRUE)
groups <- split(x, factors)

# calculate observed probabilities (frequencies)
counts <- lapply(groups, FUN=length)
observed <- lapply(counts, function(x) x / n)
observed <- as.numeric(observed)

# perform goodness-of-fit test, calculate chi square test statistic
test <- chisq.test(observed, p=expected, rescale.p=TRUE)

print(test$statistic)
