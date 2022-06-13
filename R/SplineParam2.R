# BDSplineBasis = (BSplineBasis,  Dirac delta function , Log)
# last basis is the primitive of order "dirac" of the Dirac delta function at min(knots)
# if dirac == 0, last basis is the Dirac delta function

BDSplineBasis<-function(knots, degree=3, keep.duplicates=FALSE, dirac=0L, cdirac=1.0, log=FALSE, clog=ifelse(log,1.0, 0.0)) {
	# for M-splines,
	#
	# the bases are not defined outside of [kmin, kMax]
	#
	# knots are 1 boundary min knots, interior knots, 1 boundary max knots
	# boundary knots are not tested to be duplicated
	# code from orthogonalsplinebasis
	degree <- as.integer(degree)
	order<-degree+1L
	n<-length(knots)
	if( n ==2 ){
# no interior knots
		knots<-c(rep(knots[1], order), rep(knots[n], order))
	}
	else {
		if(any(table(knots[2:(n-1)])>1)&&!keep.duplicates){
			warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
			# modif MGk 06/06/2011 to keep the (order-1) first and last knots 
			oknots <- knots
			iknots<-unique(oknots[2:(n-1)])
			knots<-c(rep(oknots[1], order), iknots, rep(oknots[n], order))
		}
		else {
			# just duplicates first and last knots
			oknots <- knots
			iknots<-oknots[2:(n-1)]
			knots<-c(rep(oknots[1], order), iknots, rep(oknots[n], order))
		}            
	}
	# recompute n number of knots
	n <- length(knots)
	q <- n-order
	
	SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=order, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	new("BDSplineBasis", knots=knots, min=knots[1], max=knots[length(knots)],
			degree=as.integer(degree), nbases=as.integer(q+1), Matrices=M, SplineBasis=SB, log=log, dirac=as.integer(dirac), cdirac=cdirac, clog=clog)
}


MDSplineBasis<-function(knots, degree=3, keep.duplicates=FALSE, dirac=0L, cdirac=1.0, log=FALSE, clog=ifelse(log,1.0, 0.0) ) {
	# for M-splines basis int_boundinf^^boudsup b_i(t)=1
	# just scaled BSplineBasis
	# knots are interior knots; uplicated knots are allowed for discontinuity conditions
	# boundary knots are not tested to be duplicated
	# code using thogonalsplinebasis
	tmpobj <- BSplineBasis(knots=knots, degree=degree, keep.duplicates=keep.duplicates)
	xscale<-evaluate(integrate(tmpobj), max(knots), intercept=TRUE)
	tmpobj<- tmpobj * (1/as.numeric(xscale))
	tmpobj <- as(tmpobj, "MDSplineBasis")
	tmpobj@nbases <- tmpobj@nbases + 1L
	tmpobj@dirac <- as.integer(dirac)
	tmpobj@cdirac <- cdirac
	tmpobj@log <- log
	tmpobj@clog <- clog
	tmpobj
#	new("MDSplineBasis", knots=tmpobj@knots, min=tmpobj@min, max=tmpobj@max,
#			degree=tmpobj@degree, nbases=tmpobj@nbases, Matrices=tmpobj@Matrices, SplineBasis=tmpobj@SplineBasis, log=log, clog=clog)
}




# getteurs
print.BDSplineBasis<-function(object) { 
	cat(class(object),"\n")
	cat("BSpline with Dirac at min(knos) in first basis")
	cat("Order: ",object@degree+1,"\n",
			"Degree: ",object@degree,"\n",
			"Knots: ", paste(object@knots,collapse=" "),"\n",
			"Number of bases: ", object@nbases+object@log,"\n",
			"Order of the Dirac component: ", object@dirac, "\n",
			"Range: ", paste(c(object@min, object@max),collapse=" ; "),"\n",
			sep="")
	invisible(object)
}



setMethod("show","BDSplineBasis",  print.BDSplineBasis)


setGeneric("getDirac",function(object)standardGeneric("getDirac"))
setMethod("getDirac",signature("BDSplineBasis"),function(object)object@dirac)



FEvaluateBDBasis<-function(object, x, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(x))
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	basis <- .Call(C_eval_spline_basis, as.double(knots), as.integer(order), M, 
			as.integer(intercept), as.double(x), as.integer(outer.ok))
	if(getDirac(object)>0){
		bdirac <- ifelse(x > object@min, object@cdirac*(x-object@min)^(object@dirac-1L), 0)
	}
	else {
		bdirac <- ifelse(x == object@min, object@cdirac, 0)
	}
	
	if (object@log) {
		return(cbind(basis, bdirac , log(x)))
	}
	else {
		return(cbind(basis, bdirac))
	}
}

EvaluateBDBasis<-function(object,x,intercept=TRUE, xname=NULL, outer.ok=TRUE, namespline= "B", ...) {
	
	basis <- FEvaluateBDBasis(object=object, x=x,
			intercept=intercept, 
			xname=xname,
			outer.ok= outer.ok, ...) 
	
	Aknots<-object@knots
	degree<-object@degree
	nbases<-object@nbases
	
	
	#add dimnames
	if(!is.null(xname)){
		if (intercept) {
			dbs <- c(paste(namespline, "-", xname, ":", 1:(getNBases(object)-1), sep=""), "B_Dirac_min" )
		}
		else {
			dbs <- c(paste(namespline, "-", xname, ":", 2:(getNBases(object)-1), sep=""), "B_Dirac_min" )
		}
		if( getLog(object) ){
			dbs <- c(dbs, "B_log" )
		}
		dimnames(basis)[[2]] <- dbs
	}
	
	
#  dimnames(basis) <- list(nx, 1L:n.col)
	a <- list(degree = degree, knots =  Aknots, 
			intercept = intercept, log=getLog(object), dirac=getDirac(object))
	attributes(basis) <- c(attributes(basis), a)
	class(basis) <- c("bs", "basis", "matrix")
	basis
}




setMethod("evaluate", signature("BDSplineBasis","numeric"),function(object, x, ...) EvaluateBDBasis(object=object, x=x, ...))
setMethod("fevaluate",signature("BDSplineBasis","numeric"),function(object, x, ...)FEvaluateBDBasis(object=object, x=x, ...))


PredictBDBasisBeta <- function(object=object, x=x, beta=beta, intercept=TRUE, outer.ok=TRUE, ...){
	if(intercept){
		predict(object * beta, x, intercept=TRUE, ...)
	}
	else {
		predict(object * c(0, beta), x,  intercept=FALSE, ...)
	}
}

predict.BDSplineBasis <- function(object, x, intercept=TRUE, outer.ok=TRUE, ...){
	stopifnot(is.numeric(x))
	
	cl <- predict(object=as(object, "BSplineBasis"),  x=x, intercept=intercept, outer.ok=outer.ok, ...)
	if(getDirac(object)>0){
		cl <- cl + ifelse(x > object@min, object@cdirac*(x-object@min)^(object@dirac-1L), 0)
	}
	else {
		cl <- cl + ifelse(x == object@min, object@cdirac, 0)
	}
	if (object@log) {
		return(cl + log(x) * object@clog)
	}
	else {
		return(cl)
	}
	
}



setMethod("predictSpline",signature(object="BDSplineBasis",x="numeric", beta="missing"),function(object, x, ...)predict.BDSplineBasis(object=object, x=x,  ...))



integrate.BDSplineBasis<-function(object){
	SB <- orthogonalsplinebasis::integrate(object@SplineBasis)
	new("BDSplineBasis", knots=SB@knots, min=object@min, max=object@max,
			degree=as.integer(SB@order-1), nbases=as.integer(dim(SB@Matrices)[2]+1), Matrices=SB@Matrices,
			SplineBasis=SB, log=FALSE, dirac=object@dirac+1L,cdirac=object@cdirac/(object@dirac+1) )
}

integrate.MDSplineBasis<-function(object){
	SB <- orthogonalsplinebasis::integrate(object@SplineBasis)
	new("MDSplineBasis", knots=SB@knots, min=object@min, max=object@max,
			degree=as.integer(SB@order-1), nbases=as.integer(dim(SB@Matrices)[2]+1), Matrices=SB@Matrices,
			SplineBasis=SB, log=FALSE, dirac=object@dirac+1L,cdirac=object@cdirac/(object@dirac+1) )
}

# compute values of integrated basis
IntegrateBDBasis<-function(object,x,intercept=TRUE, xname=NULL, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(x))
	if (object@log) {
		stop("no method 'integrate' for BDSplineBasis with additional log basis" )
	}
	else {
		evaluate(integrate(object), x=x, intercept=intercept, outer.ok=outer.ok, ...)
	}
}


setMethod("integrate",signature("BDSplineBasis", "missing"),integrate.BDSplineBasis)
setMethod("integrate",signature("MDSplineBasis", "missing"),integrate.MDSplineBasis, valueClass = "MDSplineBasis")
setMethod("integrate",signature("BDSplineBasis", "numeric"),function(object, x, ...)IntegrateBDBasis(object=object, x=x, ...))
setMethod("integrate",signature("MDSplineBasis", "numeric"),function(object, x, ...)IntegrateBDBasis(object=object, x=x, ...))


######################################################################
# deriv method
# if dirac == 0, no derivative
deriv.BDSplineBasis<-function(expr){
	if(expr@dirac > 0){
		SB <- orthogonalsplinebasis::deriv(expr@SplineBasis)
		
		no <- new("BDSplineBasis", knots=SB@knots, min=expr@min, max=expr@max,
				degree=as.integer(SB@order-1), nbases=as.integer(dim(SB@Matrices)[2]+1), Matrices=SB@Matrices, SplineBasis=SB,
				dirac = expr@dirac-1L, cdirac =  expr@cdirac *  expr@dirac, log=FALSE)
		
		no
	}
	else {
		warning("tentative to derive a Dirac delta function")
		return(NULL)
	}
}


setMethod("deriv",signature("BDSplineBasis"),deriv.BDSplineBasis)

InitCoefBDBasis<-function(object,ncol, init=1, intercept=TRUE, xname=NULL, ...) { 
	# output matrix with all "init" 
	# knots are all knots (first and last replicated
	stopifnot(is.integer(ncol))
	nb <- object@nbases + 1 - intercept  + object@log
	
	matrix(init, ncol=ncol , nrow=nb)
	
}


setMethod("initcoef",
		signature("BDSplineBasis","integer"),
		function(object, ncol, ...)InitCoefBDBasis(object=object, ncol=ncol, ...)
)



######################################################################
# operators

Prod.BDS.n <- function(e1, e2) {
	if(length(e2) == 1){
		e1@SplineBasis <- e1@SplineBasis * e2
		e1@Matrices <- e1@SplineBasis@Matrices
		e1@cdirac <- e1@dirac * e2
	}
	else {
		e1@SplineBasis <- e1@SplineBasis * e2[-length(e2)]
		e1@Matrices <- e1@SplineBasis@Matrices
		e1@cdirac <- e1@cdirac * e2[length(e2)]
	}
	e1
}
setMethod("*", signature(e1 = "BDSplineBasis", e2 = "numeric"), Prod.BDS.n)
setMethod("*", signature(e1 = "MDSplineBasis", e2 = "numeric"), Prod.BDS.n)

Prod.n.BDS <- function(e1, e2) {
	e2 * e1
}
setMethod("*", signature(e1 ="numeric" , e2 = "BDSplineBasis"), Prod.n.BDS)
setMethod("*", signature(e1 ="numeric" , e2 = "MDSplineBasis"), Prod.n.BDS)


Sum.BDS.BDS <- function(e1, e2) {
	if(e1@dirac != e2@dirac){
		stop("tentative to add BDSplineBasis of different order of integration")
	}
	else {
		e1@SplineBasis <- e1@SplineBasis + e2@SplineBasis 
		e1@Matrices <- e1@SplineBasis@Matrices
		e1@cdirac <- e1@cdirac + e2@cdirac 
		
		e1
	}
}
setMethod("+", signature(e1 = "BDSplineBasis", e2 = "BDSplineBasis"), Sum.BDS.BDS)
setMethod("+", signature(e1 = "MDSplineBasis", e2 = "MDSplineBasis"), Sum.BDS.BDS, valueClass = "MDSplineBasis")


Sum.BDS.n <- function(e1, e2) {
	e1 + ( BDSplineBasis(knots=e1@knots, degree=e1@degree, keep.duplicates=TRUE, dirac=e1@dirac, cdirac=0) * e2)
}
setMethod("+", signature(e1 = "BDSplineBasis", e2 = "numeric"), Sum.BDS.n)
setMethod("+", signature(e1 = "MDSplineBasis", e2 = "numeric"), Sum.BDS.n, valueClass = "MDSplineBasis")

Sum.n.BDS <- function(e1, e2) {
	e2 + e1
}
setMethod("+", signature(e1 ="numeric" , e2 = "BDSplineBasis"), Sum.n.BDS)
setMethod("+", signature(e1 ="numeric" , e2 = "MDSplineBasis"), Sum.n.BDS, valueClass = "MDSplineBasis")



Dif.BDS.BDS <- function(e1, e2) {
	e1@SplineBasis <- e1@SplineBasis - e2@SplineBasis 
	e1@Matrices <- e1@SplineBasis@Matrices
	e1@cdirac <- e1@cdirac - e2@cdirac 
	e1
}
setMethod("-", signature(e1 = "BDSplineBasis", e2 = "BDSplineBasis"), Dif.BDS.BDS)
setMethod("-", signature(e1 = "MDSplineBasis", e2 = "MDSplineBasis"), Dif.BDS.BDS, valueClass = "MDSplineBasis")

Dif.BDS.n <- function(e1, e2) {
	e1 + ((-1) * e2)
}
setMethod("-", signature(e1 = "BDSplineBasis", e2 = "numeric"), Dif.BDS.n)
setMethod("-", signature(e1 = "MDSplineBasis", e2 = "numeric"), Dif.BDS.n, valueClass = "MDSplineBasis")

Dif.n.BDS <- function(e1, e2) {
	( e2 * (-1)) + e1
}
setMethod("-", signature(e1 = "numeric", e2 = "BDSplineBasis"), Dif.n.BDS)
setMethod("-", signature(e1 = "numeric", e2 = "MDSplineBasis"), Dif.n.BDS, valueClass = "MDSplineBasis")


# no product matrix %*% BDSplineBasis 

#WCE_dirac(t) = sum_{fromt[i]<t} W_i beta basedirac(t - fromt[i]) 
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
PredictWCEMDBasis<-function(object, t, Increment, fromT, tId, FirstId, LastId, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(t))
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	wce <- .Call("predict_wce_spline_basis", as.double(knots), as.integer(order), M,
			as.integer(intercept),
			as.double(Increment), as.double(fromT), as.integer(FirstId), as.integer(LastId), 
			as.double(t),  as.integer(tId),  as.integer(outer.ok))
	wcedirac <- rep(0, length(t))
	for(i in 1:length(t)){
		wcedirac[i] <- t(Increment[FirstId[tId[i]]:LastId[tId[i]]]) %*%
				ifelse((t[i] > fromT[FirstId[tId[i]]:LastId[tId[i]]]),
						(t[i] - fromT[FirstId[tId[i]]:LastId[tId[i]]])^(object@dirac-1L),
						0.0)
	}
	wcedirac <- object@cdirac * wcedirac
	
	if (object@log) {
		stop("option \"log=TRUE\" unavailable for Weighted Cumulative Exposure effects")
	}
	else {
		return(wce + wcedirac)
	}
	
}


setMethod("predictwce",signature("MDSplineBasis", "numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...) PredictWCEMDBasis(object=object, t=t, Increment=Increment, fromT=fromT, ...))


setMethod("predictwce",signature("BDSplineBasis", "numeric","numeric","numeric"),
		function(object, t, Increment, fromT, ...) PredictWCEMDBasis(object=object, t=t, Increment=Increment, fromT=fromT, ...))

# no specific function for predcumwce because it calls predictwce() 
setMethod("predictcumwce",signature("MDSplineBasis", "numeric","numeric","numeric"),
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



# MD-spline
GradientWCEMDBasis<-function(object, t, Increment, fromT, tId, FirstId, LastId, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(t))
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	
	grdwce <- .Call("grad_wce_spline_basis", as.double(knots), as.integer(order), M,
			as.integer(intercept),
			as.double(Increment), as.double(fromT), as.integer(FirstId), as.integer(LastId), 
			as.double(t),  as.integer(tId),  as.integer(outer.ok))
	wcedirac <- rep(0, length(t))
	for(i in 1:length(t)){
		wcedirac[i] <- t(Increment[FirstId[tId[i]]:LastId[tId[i]]]) %*%
				ifelse((t[i] > fromT[FirstId[tId[i]]:LastId[tId[i]]]) ,
						(t[i] - fromT[FirstId[tId[i]]:LastId[tId[i]]])^(object@dirac-1L),
						0.0)
	}
	
	if (object@log) {
		stop("option \"log=TRUE\" unavailable for Weighted Cumulative Exposure effects")
	}
	else {
		return(cbind(grdwce, wcedirac))
	}
	
}

setMethod("gradientwce",
		signature(object="MDSplineBasis", t="numeric", Increment="numeric", fromT="numeric"),
		function(object, t, Increment, fromT, ...) GradientWCEMDBasis(object=object, t=t, Increment=Increment, fromT=fromT, ...))

