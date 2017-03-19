power <- function(dfin, alpha)
{
	# H0
	dfnull <- subset(dfin, H=='norm')

	# H1
	dfalt <- subset(dfin, H=='pgnorm')

	# empirical quantiles
	q <- quantile(dfnull$x, probs=1-alpha)

	F1 <- ecdf(dfalt$x)
	# F1 is function
	beta <- F1(q)

	df <- data.frame(alpha=alpha, power=1-beta)
	return(df)
}
