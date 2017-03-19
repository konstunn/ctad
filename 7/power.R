power <- function(dfin, alpha, h, side='right')
{
	dfnull <- subset(dfin, h0 == h[1])
	dfalt <- subset(dfin, h0 == h[2])

	F1 <- ecdf(dfalt$x)

	if (side == 'left') {
		Sa <- quantile(dfnull$x, probs=alpha)
		beta <- 1 - F1(Sa)
	} else if (side == 'right') {
		Sa <- quantile(dfnull$x, probs=1-alpha)
		beta <- F1(Sa)
	} else if (side == 'both') {
		Sl <- quantile(dfnull$x, probs=alpha/2)
		Sr <- quantile(dfnull$x, probs=1-alpha/2)
		beta <- F1(Sr) - F1(Sl)
	}

	df <- data.frame(alpha=alpha, power=1-beta)
	return(df)
}
