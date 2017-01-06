
n <- 100

x <- rcauchy(n)

# number of bins
k <- 5

x_i <- quantile(x, seq(from=0, to=1, by=1/k))

# if hypothesis is composite, then
# estimate these parameters first
p <- pcauchy(x_i, location=0, scale=1)

expected <- zoo::rollapply(p, width=2, function(x) x[2] - x[1])

# group sample
f <- cut(x, x_i, include.lowest=TRUE)
groups <- split(x, f)
counts <- lapply(groups, FUN=length)
observed <- lapply(counts, function(x) x / n)
observed <- as.numeric(observed)

test <- chisq.test(observed, p=expected, rescale.p=TRUE)

print(test$statistic)
