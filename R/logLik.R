
## log-likelihood for flexrsurv objects
logLik.flexrsurv <- function(object, ...)
{
	if(!missing(...)) warning("extra arguments discarded")
	val <- object$loglik
	## Note: weights are not implemented
	attr(val, "nobs") <- object$nobs
	attr(val, "df") <-  length(coef(object))
	class(val) <- "logLik"
	val
}

nobs.flexrsurv <- function(object, ...) object$nobs 


vcov.flexrsurv <- function (object, ...) {
	vmat <- object$var
	vname <- names(object$coefficients)
	dimnames(vmat) <- list(vname, vname)
	vmat
}
