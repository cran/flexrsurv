predictHazard <- function(object, ...) UseMethod("predictHazard")

predictHazard.flexrsurv <- function(object, newdata = NULL,
		type=c("lp", "link", "terms", "risk", "hazard", "hazardrate", "rate", "loghazard", "log", "lograte"),
		se.fit = FALSE, 
		na.action = na.pass, ...){
	
	
	call <- match.call()
	type <-match.arg(type)
	
	type <- switch(type,
			risk = "risk",
			rate = "risk",
			hazard = "risk",
			hazardrate = "risk",
			lp = "link",
			log = "link",
			loghazard = "link",
			lograte = "link",
			link = "link",
			terms = "terms")  
	
	if (type == "terms") {
		stop("not yet implemented")
	}
	
	if ((missing(newdata) || is.null(newdata)) && !se.fit &&  type != "terms" ){
		pred <- switch(type, link = object$linear.predictors, 
				risk = object$fitted.values)
		if (!is.null(na.action)) {
			pred <- napredict(na.action, pred)
		}
	} else {
		if (missing(newdata) || is.null(newdata)) {
			newdata <- object$data
			if (!is.null(na.action)) {
				newdata <- na.action(newdata)
			}
		} else {
			newdata <- as.data.frame(newdata)
		}
		
		Terms <- object$terms
		
		# only the end times are needed, add beginning and stuts vars if needed
		survvars <- all.vars(lhs(formula(Terms)))
		if(length(survvars)==3){
# original survival object (tfrom , tto)
# if tfrom missing in newdata, add new columns with tfrom=0
			if(!(survvars[1] %in% dimnames(newdata)[[2]])){
				newdata[[survvars[1]]] <- rep(0, dim(newdata)[1])
			}
		}
		if(!(survvars[length(survvars)] %in% dimnames(newdata)[[2]])){
			newdata[[survvars[length(survvars)]]] <- rep(1, dim(newdata)[1])
		}
		
		m <- model.frame(Terms, data=newdata)
		
		
		Y <- model.extract(m, "response")
		# get T-splines
		Max_T <- eval(as.expression(object$call$Max_T))
		Min_T <- eval(as.expression(object$call$Min_T))
		Spline <- eval(as.expression(object$call$Spline))
		if(is.null(Max_T)){
			Max_T <- max(Y[,1:(ncol(Y)-1)])
		}
		if(is.null(Min_T)){
			Min_T <- 0
		}
		
		
		if(!is.null(object$baselinehazard)){
			baselinehazard <- object$baselinehazard
		} else if (!is.null(eval(as.expression(object$call$baselinehazard)))){
			baselinehazard <- eval(as.expression(object$call$baselinehazard))
		} else {
			baselinehazard <- TRUE
		}
		if(baselinehazard == TRUE){
			bhlink <- object$bhlink
			knots.Bh <- eval(as.expression(object$call$knots.Bh)) 
			degree.Bh <- eval(as.expression(object$call$degree.Bh)) 
			intercept.Bh <- eval(as.expression(object$call$intercept.Bh)) 
			if(!is.null(intercept.Bh)){
				Intercept_t0 <- intercept.Bh				
			} else {
				Intercept_t0 <- TRUE				
			}
			if(Spline=="b-spline"){
				Spline_t0 <- BSplineBasis(knots=c(Min_T, knots.Bh, Max_T),
						degree=degree.Bh,
						keep.duplicates=TRUE)
				nT0basis <- getNBases(Spline_t0) - 1 +  Intercept_t0
				ngamma0 <- nT0basis
				
				Spline_t<- BSplineBasis(knots=c(Min_T, knots.Bh, Max_T),
						degree=degree.Bh,
						keep.duplicates=TRUE)
			} else if(Spline=="tp-spline") {
				Spline_t0 <-TPSplineBasis(knots=knots.Bh,
						degree=degree.Bh,
						min=Min_T,
						max=Max_T)

				nT0basis <- getNBases(Spline_t0) - 1 +  Intercept_t0
				ngamma0 <- nT0basis
				Spline_t <-TPSplineBasis(knots=knots.Bh,
						degree=degree.Bh,
						min=Min_T,
						max=Max_T)
			}
		}
		else {
			Spline_t0 <- NULL
			Intercept_t0 <- TRUE
			nT0basis <- 0
			ngamma0 <- 0
			knots.Bh <- NULL 
			degree.Bh <- NULL 
			bhlink <- "log"
		}
		
		# extract components
		des <- object$des
		newdes <- ReadDesignFlexrsurv(Terms=Terms, modframe=m, data=newdata, rate=NULL, Spline_t0=Spline_t0, intercept.Bh=Intercept_t0)
		
		# test if all time spline are identical
		if(!is.null(des$X) | !is.null(des$Z)){
			allSpline_T <- c(des$Spline_XT, des$Spline_ZT)
			if(length(allSpline_T) > 1){
				for(sb in allSpline_T[-1]){
					if(!identical(sb, allSpline_T[[1]])){
						stop("flexrsurv cannot handle different spline basis for time dependent effects")
					}
				}
			}
			Spline_t <- allSpline_T[[1]] 
		}
		
		
		
		# X0 linear and non linear effects
		X0<-des$X0
		if(!is.null(X0)){
			is.PH <- TRUE
			if(is.matrix(X0)) {
				nX0 <- dim(X0)[2]
			} else if(is.vector(X0)) {
				nX0 <- 1L
			} else {
				stop("error flexrsurv_LL.fit(): wrong type of X0") 
			}
			nalpha0<-nX0
			ialpha0<-1:nX0 + ngamma0
			Ialpha0<-1:nX0 
			First.alpha0<-ngamma0+1
			first.alpha0<-1
			X0<-newdes$X0
		} else {
			is.PH <- FALSE
			nX0 <- 0L
			nalpha0 <- 0L
			alpha0 <- NULL
			ialpha0 <- NULL
			Ialpha0 <- NULL
			First.alpha0 <- NULL
			first.alpha0 <- NULL
		}
		
		
		
# X (NPH effects)
		if(!is.null(des$X)){
			X<-as.matrix(des$X)
			#  Intercept.t of NPH() effets are set in Flersurv:fix.flexrsurv.formula
			#  now build the vector of Intercept_t_NPH for each X var
			# {linear or non linear} and non prop effect are possible; set intercept for NPH effects
			Intercept_t_NPH <- rep(TRUE, length(des$XVars))
			for (i in attr(des$TermsX, "specials")[["NPH"]]){
				thecall <-  match.call(NPH, attr(des$TermsX,"variables")[[i+1]])
				# the variable name is the second argument of the special function
				Intercept_t_NPH[i] <- ifelse( length(thecall[["Intercept.t"]]) == 0 ,
						formals(NPH)[["Intercept.t"]],
						thecall[["Intercept.t"]] )
			}
			nTbasis <- getNBases(Spline_t)
			nTbasis_NPH <- getNBases(Spline_t) - 1 + Intercept_t_NPH
			if(is.matrix(X)) {
				nX <- dim(X)[2]
			} else if(is.vector(X)) {
				nX <- 1L
			}
			
			nbeta0 <- sum(nTbasis_NPH)
			ibeta0 <- 1:nbeta0 + ngamma0 + nalpha0
			Ibeta0 <- 1:nbeta0 
			First.beta0 <- ngamma0 + nalpha0 + 1
			first.beta0 <- 1
			X<-as.matrix(newdes$X)
		} else  {
			X <- NULL
			Intercept_t_NPH <- NULL
			is.NPHLIN <- FALSE
			nTbasis <- 0L
			nTbasis_NPH <- 0L
			nX <- 0L
			nbeta0 <- 0L
			beta0 <- NULL
			ibeta0 <- NULL
			Ibeta0 <- NULL
			First.beta0 <- NULL
			first.beta0 <- NULL
		}
		
		
		# Z (NPHNLL effects
		# get splines for each NPHNLL variable
		
		if( !is.null(des$Z) ){
			Spline_Z <- des$Spline_Z
			Z<-DesignMatrixNPHNLL(Z=des$Z, listsplinebasis=Spline_Z, timesplinebasis=Spline_t)
			
			if( getNvar(Z)>0){
				is.NPHNLIN <- TRUE
				nTbasis <- getNBases(Spline_t)
				nZ<-getNvar(Z)
				nalpha <- getNparam(Z)
				ialpha <- 1:nalpha + ngamma0 + nalpha0 + nbeta0  
				Ialpha <- 1:nalpha  + nalpha0   
				First.alpha <- ngamma0 + nalpha0 + nbeta0 + 1  
				first.alpha <- nalpha0  + 1  
				
				# as first beta is constraints, nTbasis -1 beta per Z variable
				nbeta <- nZ * (nTbasis-1)
				ibeta <- 1:nbeta + ngamma0 + nalpha0 + nbeta0 + nalpha
				Ibeta <- 1:nbeta  + nbeta0 
				First.beta <- ngamma0 + nalpha0 + nbeta0 + nalpha + 1
				first.beta <- 1 + nbeta0
				
				Z <- DesignMatrixNPHNLL(Z=newdes$Z, listsplinebasis=Spline_Z, timesplinebasis=Spline_t)
				
			}
		} else  {
			is.NPHNLIN <- FALSE
#    nTbasis <- 0 already done with beta0
#      Z <- new("DesignMatrixNPHNLL")
			Z <- NULL
			nZ <- 0L
			nalpha <- 0L
			alpha <- NULL
			Ialpha <- NULL
			ialpha <- NULL
			First.alpha <- NULL
			first.alpha <- NULL
			nbeta <- 0L
			beta <- NULL
			Ibeta <- NULL
			ibeta <- NULL
			First.beta <- NULL
			first.beta <- NULL
		}
		
# get coefficients in order baseline hazard | linear & non-linear effects | non Prop efefct | nonlin_nonprop effects 
		# coefficient <- param[des$param2coef]
		# param[des$param2coef] <- coefficients
		alltheparameters <- rep(0, length(object$coefficients))
		alltheparameters[des$param2coef] <- object$coefficients
		# computes linear predictors
		
		if ( type != "terms" ){
			linpred <- .computeLinearPredictor_GA0B0AB(allparam=alltheparameters,
					Y=Y, X0=X0, X=X, Z=Z,
					nT0basis=nT0basis,
					Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
					ialpha0=ialpha0, nX0=nX0,
					ibeta0= ibeta0, nX=nX, 
					ialpha=ialpha, 
					ibeta= ibeta, 
					nTbasis=nTbasis,
					Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
					bhlink=bhlink,
					debug=10000)    
			if (se.fit) {
				stderr <- .computeStdErrorLinearPredictor_GA0B0AB(allparam=alltheparameters,
						var=object$var,
						Y=Y, X0=X0, X=X, Z=Z,
						nT0basis=nT0basis,
						Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0= ibeta0, nX=nX, 
						ialpha=ialpha, 
						ibeta= ibeta,
						listZSplineBasis = des$Spline_Z,
						nTbasis=nTbasis,
						Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
						bhlink=bhlink,
						debug=FALSE)
			}
		} else {
#      linearPredictors <- .computeTermsLinearPredictor_GA0B0AB(allparam=alltheparameters,
#                                                                            var=object$var,
#                                                               Y=Y, X0=X0, X=X, Z=Z,
#                                                               nT0basis=nT0basis,
#                                                               Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
#                                                               ialpha0=ialpha0, nX0=nX0,
#                                                               ibeta0= ibeta0, nX=nX, 
#                                                               ialpha=ialpha, 
#                                                               ibeta= ibeta, 
#                                                               nTbasis=nTbasis,
#                                                               Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
#                                                               terms = terms,
#                                                               debug=FALSE)
#      
#      if (!se.fit) {
#        stdErrorLinearPredictors <- .computeStdErrorTermsLinearPredictor_GA0B0AB(allparam=alltheparameters,
#                                                                                 Y=Y, X0=X0, X=X, Z=Z,
#                                                                                 nT0basis=nT0basis,
#                                                                                 Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
#                                                                                 ialpha0=ialpha0, nX0=nX0,
#                                                                                 ibeta0= ibeta0, nX=nX, 
#                                                                                 ialpha=ialpha, 
#                                                                                 ibeta= ibeta, 
#                                                                                 nTbasis=nTbasis,
#                                                                                 Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
#                                                                                 terms = terms,
#                                                                                 debug=FALSE)
#      }
		} 
		
		if (type=="risk"){
			if(bhlink == "log"){
				linearPredictors <- pmax(exp(linpred ), .Machine$double.eps)
				if (se.fit) {
					stdErrorLinearPredictors <- stderr * linearPredictors
				}
			} else {
				linearPredictors <- pmax(linpred[,1] * exp(linpred[,2] ), .Machine$double.eps)
				if (se.fit) {
					varerr <- attr(stderr, "varerr")
					stdErrorLinearPredictors <- linearPredictors * (varerr[,1]/linpred[,1]^2 + 2*varerr[,2]/linpred[,1] + varerr[,3])
				}
			}
		} else {
			linearPredictors <- linpred
			if (se.fit) {
				stdErrorLinearPredictors <- stderr
			}
		}
		
		if (se.fit) {
			# structure similar to predic.lm() 
			pred <- list(fit=linearPredictors , se.fit=stdErrorLinearPredictors)
		} else {
			pred <- linearPredictors
		}
		
	}
	
	pred
	
	
}




