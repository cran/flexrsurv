flexrsurv.optim.control <- function (maxit = 5000, trace=100, REPORT=1, ...) 
{
	rval <- list(maxit = maxit, trace = trace, REPORT = REPORT)
	rval <- c(rval, list(...))
	if (!is.null(rval$fnscale)) 
		message("optim() control parameter fnscale must not be modified")
	rval$fnscale <- -1
	if (is.null(rval$reltol)) 
		rval$reltol <- .Machine$double.eps^(1/1.3)
	rval
}
