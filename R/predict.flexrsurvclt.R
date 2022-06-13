## S3 method for class 'flexrsurv'
predict.flexrsurvclt <- function(object, newdata = NULL,
		type = c("lp", "link", "risk", "hazard", "hazardrate", "rate", "loghazard", "log",
				"lograte", "cumulative.rate", "cumulative.hazard", "cumulative", "cum", "survival", "surv", "netsurv", "clt", "correction"),
		se.fit = FALSE,
		na.action = na.pass, ...){
	
	type <- match.arg(type)
	
	if (type == "clt" | type == "correction" ){
		stop ('type = "clt" or "correction" not yet implemented')
	} else {
		# call to predict.flexrsurv
		# nb parameters for the correction model:
		
		
		ndfexcess <- object$ndf$ndf.excess
		object$coefficients <- coef(object)[1:ndfexcess]
		object$var <- object$var[1:ndfexcess,1:ndfexcess]
		
		return(NextMethod("predict", object, type=type, se.fit=se.fit, na.action=na.action, ...) )
	}
	
}


