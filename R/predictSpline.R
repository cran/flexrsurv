# method to evaluate a spline function

setMethod("predictSpline",
		signature(object="character",x="numeric", beta="numeric"),
		function(object, x, beta, ...)predictSpline.type(type=object, x=x,  coef = beta, ...))

predictSpline.type <- function(type, x, knots, degree, keep.duplicates = FALSE, coef = 1){
	
	knots <- sort(knots)
	
	Spline <-switch(type,
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
		

	
	pred <- predictSpline(Spline*coef, x )
	
	return(pred)
	
}

