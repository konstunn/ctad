#!/usr/bin/Rscript --vanilla

# Kolmogorov test statistic
# TODO: calc by ks.test
# distributed as kolmim::pkolm
Sk <- function(X, P) {
	X <- sort(X)
	n <- length(X)
	i <- seq(1, n)
	Dplus <- max(i/n - P(X, theta))
	Dminus <- max(P(X, theta) - (i-1)/n)
	Dn <- max(Dplus, Dminus)
	Sk <- (6*n*Dn + 1)/(6 * sqrt(n))
	return(list(statistic=Sk))
}

# Smirnov test statistic
# distributed as chi-squared, df = 2
smirnov.test <- function(X, P, theta) {
	X <- sort(X)
	n <- length(X)
	i <- seq(1, n)
	Dplus <- max(i/n - P(X, theta))
	Sm <- (6*n*Dplus+1)^2/9*n
	return(list(statistic=Sm))
}

# cvm.cdf
# cvm.test

# ADGofTest::ad.test
# distributed as ad.cdf

statmod.gof2cauchy <- function(n=100, N=16600, Htype, trueH, testS)
{
	sink(stderr())
	message("")
	message(deparse(sys.calls()[[sys.nframe()]]))
	true_H <- paste('r', trueH, sep='')
	message("Random number samples generating...")
	message("Generating big sample...")
	X <- do.call(true_H, args=list(n=n*N))
	message("Splitting big sample into list of samples...")
	X <- split(X , ceiling(seq_along(X)/n))

	library(fitdistrplus)
	library(plyr)

	test <- tolower(testS)
	test <- paste(test, '.test', sep='')

	if (Htype == 'simple') {
		message("Simple hypothesis.")
		estimates <- c(rep(c(location=0,scale=1), N))
		estimates <- split(estimates, ceiling(seq_along(estimates)/2))
	} else if (Htype == 'complex') {
		message("Complex hypothesis.")
		message("Estimating parameters...")
		lfitdist <- function(x) {
			S <- function(theta, x) {
				rez <- do.call(test, args=list(x, 'pcauchy', location=theta[1],
										scale=theta[2]))$statistic
				return(rez)
			}
			rez <- optim(c(0,1), S, x=x)$par
			return(rez)
		}
		estimates <- llply(X, lfitdist, .progress='text')
	} else stop("Invalid Htype value: must be 'simple' or 'complex'")

	# calculate test statistic values
	message("Calculating test statistic array...")
	ltest <- function(x, t) do.call(test, args=list(x, 'cauchy', location=t[1],
													scale=t[2]))$statistic
	# TODO: cast X to matrix
	X <- matrix(X, nrow=N, byrow=TRUE)
	estimates <- matrix(estimates, nrow=N, byrow=TRUE)
	chisq <- maply(X, ltest, t=estimates, .progress='text')
	sink()
	return(data.frame(x=chisq, k=rep(k,N), n=rep(n,N), N=rep(N,N),
					  trueH=rep(trueH,N), Htype=rep(Htype,N)))
}
