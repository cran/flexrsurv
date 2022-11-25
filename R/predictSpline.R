# in AllGeneric.R
# predictSpline <- function(object, x, ...) UseMethod("predictSpline")

setMethod("predictSpline",
		signature(object="character",x="numeric"),
		function(object, x, ...)predictSpline.default(object=object, x=x, ...))


predictSpline.default <- function(object=c("b-spline", "tp-spline"), x, beta = 1, knots, degree, keep.duplicates = FALSE, ...){

	object <-match.arg(object)
	
	knots <- sort(knots)
	
	Spline <-switch(object,
			"b-spline" =  BSplineBasis(knots=knots, 
					degree=degree, 
					keep.duplicates=keep.duplicates, 
					log=FALSE),
			"tp-spline" =  TPSplineBasis(knots=knots[c(-1, -length(knots))], 
					degree=degree, 
					min=min(knots),
					max=max(knots),
					log=FALSE,
					type="standard"))
		

	
	pred <- predictSpline(Spline*beta, x )
	
	return(pred)
	
}

