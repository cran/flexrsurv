PivotConstraintSplinceCLT.R2bBSplineBasis <- function(object){
	knots <- unique(getKnots(object))
	nk <- length(knots)
	ndf <- getNBases(object)
	ndeg <- getDegree(object)
	
	if(nk == ndf){
		pivot <- knots
	} else if(nk == (1+ndf)){
		# middle of interknots interval
		pivot<-sort(c(knots[1], (knots[1:(nk-1)]+knots[2:nk]/2), knots[nk]))
	} else if( ndf == 2*nk -1){
		# knots + midle of intervale
		pivot <- sort(c(knots, (knots[1:(nk-1)]+knots[2:nk]/2), knots[nk]))
	} else if( ndf == 2*nk){
		# boundary knots + 2 pivots per interknots interval
		pivot <- sort(c(knots[1], 
						1/3*(knots[2:nk]-knots[1:(nk-1)])+knots[1:(nk-1)], 
						2/3*(knots[2:nk]-knots[1:(nk-1)])+knots[1:(nk-1)], 
						knots[nk]))
	} else {
		stop("PivotConstraintSplinceCLT.R2bBSplineBasis: combination of ndf & knots not yet implemented") 
	}
	
	return(pivot)
}





# pivot for Sturm's constraints
# evaluate(B, pivot) >= 0
# if degree > 3 SturmSpline(B) = 0 
# if degree == 3 Chan, V., et al. (2021). "Efficient estimation of smoothing spline with exact shape constraints." Statistical Theory and Related Fields 5(1): 55-69

PivotConstraintSplinceMonotoneSturm.BSplineBasis <- function(object){
	knots <- unique(getKnots(object))
	nk <- length(knots)
	ndf <- getNBases(object)
	degree <- getDegree(object)
	
	
	if( degree == 3){
		d2objatknots <- evaluate(deriv(deriv(object)), knots)
		# note that testval[i] = +Inf if d2objatknots[i]  == d2objatknots[i+1]
		testval <- d2objatknots[-nk]/(d2objatknots[-nk]-d2objatknots[-1])
		stationary_points <- na.omit(ifelse( 0 < testval & testval < 1, knots[-nk] + testval , NA ))
		pivot <-  c(knots, stationary_points)
	}
	else {
		pivot <- PivotConstraintSplinceCLT.R2bBSplineBasis(object)
		
	}
	
	return(pivot)
}


setGeneric("PivotConstraintSplinceCLT",function(object)standardGeneric("PivotConstraintSplinceCLT"))
setMethod("PivotConstraintSplinceCLT",signature("R2bBSplineBasis"),function(object)PivotConstraintSplinceCLT.R2bBSplineBasis(object=object))

# 
# output : evaluate(B, pivot) >= 0
# if degree > 3 SturmSpline(B) = 0 
# if degree == 3 Chan, V., et al. (2021). "Efficient estimation of smoothing spline with exact shape constraints." Statistical Theory and Related Fields 5(1): 55-69

ConstraintSplineNonDecreassing.BSplineBasis <- function(object){
	knots <- unique(getKnots(object))
	nk <- length(knots)
	ndf <- getNBases(object)
	degree <- getDegree(object)
	
	
	if( degree == 3){
		d2objatknots <- evaluate(deriv(deriv(object)), knots)
		# note that testval[i] = +Inf if d2objatknots[i]  == d2objatknots[i+1]
		testval <- d2objatknots[-nk]/(d2objatknots[-nk]-d2objatknots[-1])
		stationary_points <- na.omit(ifelse( 0 < testval & testval < 1, knots[-nk] + testval , NA ))
		support <- cbind(knots[-nk], knots[-1], stationary_points)
		Constraint <- apply(X=support, MARGIN = 2, FUN = function(x, spline){ min(predictSpline(deriv(spline), x[!is.na(x)])) }, spline  = object) 
	}
	else {
		Constraint <-  c(predictSpline(deriv(PivotConstraintSplinceCLT.R2bBSplineBasis(object))))	
	}
	
	return(Constraint)
}

setGeneric("ConstraintSplineNonDecreassing",function(object)standardGeneric("ConstraintSplineNonDecreassing"))
setMethod("ConstraintSplineNonDecreassing",signature("BSplineBasis"),function(object) ConstraintSplineNonDecreassing.BSplineBasis(object=object))



