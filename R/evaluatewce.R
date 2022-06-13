# methods for computing  Weighted Cummulative Exposure at t of one exposure profile (Wi, fromt) :
#
# WCE(t) = sum_{fromt[i]<t} W_i (beta %*% spline(t - fromt[i]) 
#
# t : vector of time at which the wce is computed : it is assumed that min(t) >= max(frmot)
# object : INTEGRATED spline bases object, assumed define on [0, Tmax], wdht max(t) <= Tmax
# beta : coefficient of the basis
# W : vector of exposure INCREMENTS
# fromT : vector of beginnig of time interval
# tId : index such that each t is in ] fromT[tId], toT[tId]]
#        the Id of each t is Id[tId]
#  it is assumed that for each Id, the profile is such that the last interval is right-open [fromT[last] + infty[
#             (ie that t <= toT[last])  
# returned value



# M-spline
PredictWCEMBasis<-function(object, t, Increment, fromT, tId, FirstId, LastId, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(t))
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	wce <- .Call("predict_wce_spline_basis", as.double(knots), as.integer(order), M,
			as.integer(intercept),
			as.double(Increment), as.double(fromT), as.integer(FirstId), as.integer(LastId), 
			as.double(t),  as.integer(tId),  as.integer(outer.ok))
	
	if (object@log) {
		stop("option \"log=TRUE\" unavailable for Weighted Cumulative Exposure effects")
	}
	else {
		return(wce)
	}
	
}

# EM-spline
PredictWCEEMBasis<-function(object, t, Increment, fromT, tId, FirstId, LastId, intercept=TRUE, ...) {
	stopifnot(is.numeric(t))
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	wce <- .Call("predict_wce_espline_basis", as.double(knots), as.integer(order), M,
			as.integer(intercept),
			as.double(Increment), as.double(fromT), as.integer(FirstId), as.integer(LastId), 
			as.double(t),  as.integer(tId), PACKAGE="flexrsurv")
	
	if (object@log) {
		stop("option \"log=TRUE\" unavailable for Weighted Cumulative Exposure effects")
	}
	else {
		return(wce)
	}
	
}


PredictWCETPBasis<-function(object, t, Increment, fromT, tId, FirstId, LastId, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(t))
	
	min<-object@min
	max<-object@max
	allknots<-object@knots
	order<-object@degree+1
	knots <- unique(allknots)
	replicates <- table(allknots)
	degrees <- object@degrees
	coef <- object@coef
	
	if(object@type == "increasing"){
		stop("Weighted Cumulative Exposure effects are not defined for increasing truncated power splines")
	}
	else {
		wce <- .Call("predict_wce_trunc_power_basis", as.double(knots), as.double(replicates),
				as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
				as.double(Increment), as.double(fromT), as.integer(FirstId), as.integer(LastId), 
				as.double(t),  as.integer(tId), as.integer(outer.ok), PACKAGE="flexrsurv")
	}
	
	if (object@log) {
		stop("option \"log=TRUE\" unavailable for Weighted Cumulative Exposure effects")
	}
	else {
		return(wce)
	}
}





setMethod("predictwce",signature("BSplineBasis", "numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...) PredictWCEMBasis(object=object, t=t, Increment=Increment, fromT=fromT, ...))
setMethod("predictwce",signature("MSplineBasis", "numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...) PredictWCEMBasis(object=object, t=t, Increment=Increment, fromT=fromT, ...))
setMethod("predictwce",signature("EMSplineBasis", "numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...)PredictWCEEMBasis(object=object, t=t, Increment=Increment, fromT=fromT, ...))
setMethod("predictwce",signature("TPSplineBasis","numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...)PredictWCETPBasis(object=object, t=t, Increment=Increment, fromT=fromT, ...))



#S4 method
# computes the cumulative WCE effect Cum_WCE(0, t) = int_0^t (WCE(Wn, fromT,toT, ISW, t) dt)
# assuming ISW is the integrated basis scaled by the coefficients   (beta_i * b_i)
#    Cum_WCE(0, t[m]) = sum_(i=firstId[tId[m]])^(fromT[i]<=t[m])Increment_i [W(t[m]-FromT[i] )]
# Cum_WCE(0, t) = predicWCE(integrate(object), t)
PredictCumWCE <- function(object=object, t, Increment, fromT, tId, FirstId, LastId, intercept=TRUE, outer.ok=TRUE, ...){
# object : scaled integrated spline basis of the WCE
# t      : the the times at which cumWCE is computed
# tId    : the Id to which corresponds t
# Icrement, fromT, toT, First Id, :
#          the exposure profile
# intercept : whether intercep is in the spline function
# outer.ok  : if true, wce(t) = 0 is t not in range(boundary.knots)
	
	# returned value: a vector of the same lenght than t 
	stopifnot(is.numeric(t))
	
	iobject <- integrate(object)
	cumwce <- predictwce(object=iobject,  t=t, Increment=Increment, fromT=fromT,
			tId=tId, FirstId=FirstId, LastId=LastId, intercept=intercept, outer.ok=outer.ok)
	return(cumwce)
	
}



setMethod("predictcumwce",signature("MSplineBasis", "numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...) PredictCumWCE(object=object, t=t, Increment=Increment, fromT=fromT, ...))
setMethod("predictcumwce",signature("EMSplineBasis", "numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...) PredictCumWCE(object=object, t=t, Increment=Increment, fromT=fromT, ...))
setMethod("predictcumwce",signature("TPSplineBasis", "numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...) PredictCumWCE(object=object, t=t, Increment=Increment, fromT=fromT, ...))


################################################################################

# methods for computing the gradient (with respect to the spline coefs) of the
#       Weighted Cummulative Exposure at t of one exposure profile (Wi, fromt) :
#      d WCE(t)/ d beta with 
# 
# WCE(t) = sum_{fromt[i]<t} W_i (beta %*% spline(t - fromt[i]) 
#
# t : vector of time at which the wce is computed : it is assumed that min(t) >= max(frmot)
# object : INTEGRATED spline bases object, assumed define on [0, Tmax], wdht max(t) <= Tmax
# beta : coefficient of the basis
# W : vector of exposure INCREMENTS
# fromT : vector of beginnig of time interval  
# returned value



# M-spline
GradientWCEMBasis<-function(object, t, Increment, fromT, tId, FirstId, LastId, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(t))
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	
	grdwce <- .Call("grad_wce_spline_basis", as.double(knots), as.integer(order), M,
			as.integer(intercept),
			as.double(Increment), as.double(fromT), as.integer(FirstId), as.integer(LastId), 
			as.double(t),  as.integer(tId),  as.integer(outer.ok))
	
	if (object@log) {
		stop("option \"log=TRUE\" unavailable for Weighted Cumulative Exposure effects")
	}
	else {
		return(grdwce)
	}
	
}

# EM-spline
GradientWCEEMBasis<-function(object, t, Increment, fromT, tId, FirstId, LastId, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(t))
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	
	grdwce <- .Call("grad_wce_espline_basis", as.double(knots), as.integer(order), M,
			as.integer(intercept),
			as.double(Increment), as.double(fromT), as.integer(FirstId), as.integer(LastId), 
			as.double(t),  as.integer(tId), PACKAGE="flexrsurv")
	
	if (object@log) {
		stop("option \"log=TRUE\" unavailable for Weighted Cumulative Exposure effects")
	}
	else {
		return(grdwce)
	}
	
}


GradientWCETPBasis<-function(object, t, Increment, fromT, tId, FirstId, LastId, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(t))
	
	min<-object@min
	max<-object@max
	allknots<-object@knots
	order<-object@degree+1
	knots <- unique(allknots)
	replicates <- table(allknots)
	
	
	degrees <- object@degrees
	coef <- object@coef
	
	if(object@type == "increasing"){
		stop("Weighted Cumulative Exposure effects are not defined for increasing truncated power splines")
	}
	else {
		grdwce <- .Call("grad_wce_trunc_power_basis", as.double(knots), as.double(replicates),
				as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
				as.double(Increment), as.double(fromT), as.integer(FirstId), as.integer(LastId), 
				as.double(t),  as.integer(tId), as.integer(outer.ok), PACKAGE="flexrsurv")
	}
	
	if (object@log) {
		stop("option \"log=TRUE\" unavailable for Weighted Cumulative Exposure effects")
	}
	else {
		return(grdwce)
	}
}





setMethod("gradientwce",signature("MSplineBasis", "numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...) GradientWCEMBasis(object=object, t=t, Increment=Increment, fromT=fromT, ...))
setMethod("gradientwce",signature("EMSplineBasis", "numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...)GradientWCEEMBasis(object=object, t=t, Increment=Increment, fromT=fromT, ...))
setMethod("gradientwce",signature("TPSplineBasis","numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...)GradientWCETPBasis(object=object, t=t, Increment=Increment, fromT=fromT, ...))

setMethod("gradientwce",signature("BSplineBasis", "numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...) GradientWCEMBasis(object=object, t=t, Increment=Increment, fromT=fromT, ...))
