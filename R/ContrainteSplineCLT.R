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

setGeneric("PivotConstraintSplinceCLT",function(object)standardGeneric("PivotConstraintSplinceCLT"))
setMethod("PivotConstraintSplinceCLT",signature("R2bBSplineBasis"),function(object)PivotConstraintSplinceCLT.R2bBSplineBasis(object=object))