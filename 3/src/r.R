
# sample size
n <- 100

# generate sample from Cauchy distribution
x <- rcauchy(n)

# number of intervals for grouping
k <- 5

# optimal probabilities vector step for AOG of sample from Cauchy distribution
opt_P <- 1/k

# calculate intervals break points
x_i <- quantile(x, seq(from=0, to=1, by=opt_P))

# TODO: if hypothesis should be composite, then estimate these parameters first
p <- pcauchy(x_i, location=0, scale=1)

# calculate expected probabilities
expected <- zoo::rollapply(p, width=2, function(x) x[2] - x[1])

# group sample
f <- cut(x, x_i, include.lowest=TRUE)
groups <- split(x, f)

# calculate observed probabilities
counts <- lapply(groups, FUN=length)
observed <- lapply(counts, function(x) x / n)
observed <- as.numeric(observed)

# perform test
test <- chisq.test(observed, p=expected, rescale.p=TRUE)

print(test$statistic)
