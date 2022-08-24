predictCLT <- function(...) UseMethod("predictCLT")

predictCLT.default <- function(knots, degree, newdata, newcoef , ...){

	Spline <- R2bBSplineBasis(knots=knots, degree=degree)
	
	pred <- predictSpline(Spline*c(1.0,   newcoef), newdata) 
	
	return(pred)
}

predictCLT.flexrsurvclt <- function(object, newdata = NULL,
		type=c("clt", "correction"),
		se.fit = FALSE, 
		na.action = na.pass, newcoef = NULL, ...){
	
	
	call <- match.call()
	type <-match.arg(type)
	
	type <- switch(type,
			clt = "clt",
			correction = "clt")  

	
	# il faut reconstruire le BX0 et le Spline_CLT
	paramclt <- get_SplinebasisCLT(object)


	if(!is.null(newcoef)){
		if(length(newcoef) != length(object$list_coefficients$brass0)){
			stop("wrong lenght o newcoef")
		}
	} else {
		newcoef <- object$list_coefficients$brass0
	}
	
	
	if ((missing(newdata) || is.null(newdata)) && !se.fit ){
		pred <- switch(type, clt = predictSpline(paramclt$Spline*c(1.0,   newcoef), object$logit_end) ) 
		if (!is.null(na.action)) {
			pred <- napredict(na.action, pred)
		}
	} else {
		if (missing(newdata) || is.null(newdata)) {
			newdata <- object$logit_end
			if (!is.null(na.action)) {
				newdata <- na.action(newdata)
			}
		} 
		pred <- switch(type, clt = predictSpline(paramclt$Spline*c(1.0,  newcoef ), newdata) ) 
		
		}
		
	
	return(pred)
	
	
}






