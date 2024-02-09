predictCumulativeHazard <- function(object, ...) UseMethod("predictCumulativeHazard")

predictCumulativeHazard.flexrsurv <- function(object, newdata = NULL,
		type=c("cumulative.rate", "cumulative.hazard", "cumulative", "cum", "survival", "surv", "netsurv"),
		se.fit = FALSE, 
		ci.fit = FALSE,
		level = .95,
		na.action = na.pass, ...){
# predict cumulative hazard from 0 to time  
	
	call <- match.call()
	type <-match.arg(type)
	
	type <- switch(type,
			cumulative.rate="cum",
			cumulative.hazard="cum",
			cumulative="cum",
			cum="cum",
			surv="surv",
			survival="surv",
			netsurv="surv")  
	
	if ((missing(newdata) || is.null(newdata)) && !se.fit ){
		pred <- object$cumulative.hazard
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
			hasfromto <- TRUE
			if(!(survvars[1] %in% dimnames(newdata)[[2]])){
				newdata[[survvars[1]]] <- rep(0, dim(newdata)[1])
			}
		} else {
			hasfromto <- FALSE
		}
		if(!(survvars[length(survvars)] %in% dimnames(newdata)[[2]])){
			newdata[[survvars[length(survvars)]]] <- rep(1, dim(newdata)[1])
		}
		
		
		m <- model.frame(Terms, data=newdata)
		
		Y <- model.extract(m, "response")
		
# numerical_integration_method 
		method <- object$numerical_integration_method 
		
		# numerical integration method
		# computes steps for time integtration
		if( !hasfromto){
			if(method$int_meth == "CAV_SIM"){
				int_meth <- "NC"
				intTD <- intTD_NC
				intTD_debug <- intTD_NC_debug
				intTD_base <- intTD_base_NC
				intTD_base_debug <- intTD_base_NC_debug
				intweightsfunc <- intweights_CAV_SIM
				step <- method$step
				mult <- 2
			} else if(method$int_meth == "SIM_3_8"){
				int_meth <- "NC"
				intTD <- intTD_NC
				intTD_base <- intTD_base_NC
				intweightsfunc <- intweights_SIM_3_8
				step <- method$step
				mult <- 3
			} else if(method$int_meth == "BOOLE"){
				int_meth <- "NC"
				intTD <- intTD_NC
				intTD_base <- intTD_base_NC
				intweightsfunc <- intweights_BOOLE
				step <- method$step
				mult <- 4      
			} else if(method$int_meth == "Gauss-Legendre"){
				int_meth <- "GL"
				intTD <- intTD_GL
				intTD_base<- intTD_base_GL
				intweightsfunc <-NULL
				gq <- gauss.quad(method$npoints, kind="legendre")
				step <- gq$nodes
				Nstep <- gq$weights
			} else if(method$int_meth == "GLM"){
				int_meth <- "GLM"
				intTD <- intTD_GLM
				intTD_base <- fastintTD_base_GLM
				intweightsfunc <- NULL
				step <- GLMStepParam(cuts=method$bands)
				Nstep <- WhichBand(Y[,1], step)-1L
			}
			
			if( int_meth == "NC"){
				STEPS <- cutT(Y[,1], step=method$step, mult=mult)
				Nstep <- STEPS$NstepT
				step <- STEPS$stepT
			}
		} else {
			if(method$int_meth == "CAV_SIM"){
				int_meth <- "NC"
				intTD <- intTDft_NC
				intTD_debug<- intTDft_NC_debug
				intTD_base<- intTDft_base_NC
				intTD_base_debug<- intTDft_base_NC_debug
				intweightsfunc <-intweights_CAV_SIM
				step <-method$step
				mult <- 2
			} else if(method$int_meth == "SIM_3_8"){
				int_meth <- "NC"
				intTD <- intTDft_NC
				intTD_base<- intTDft_base_NC
				intweightsfunc <-intweights_SIM_3_8
				step <-method$step
				mult <- 3
			} else if(method$int_meth == "BOOLE"){
				int_meth <- "NC"
				intTD <- intTDft_NC
				intTD_base<- intTDft_base_NC
				intweightsfunc <-intweights_BOOLE
				step <-method$step
				mult <- 4      
			} else if(method$int_meth == "Gauss-Legendre"){
				int_meth <- "GL"
				intTD <- intTDft_GL
				intTD_base<- intTDft_base_GL
				intweightsfunc <-NULL
				gq <- gauss.quad(method$npoints, kind="legendre")
				step <- gq$nodes
				Nstep <- gq$weights
			} else if(method$int_meth == "GLM"){
				int_meth <- "GLM"
				intTD <- intTDft_GLM
				intTD_base<- fastintTDft_base_GLM
				intweightsfunc <- NULL
				step <-GLMStepParam(cuts=method$bands)
				Firststep <- WhichBandInf(Y[,1], step) + 1L
				Laststep  <- WhichBand(Y[,2], step)- 1L
				Nstep <- cbind(Firststep, Laststep)
			}
			
			if( int_meth == "NC"){
				STEPS<-cutTfromto(Y[,1], Y[,2], step=method$step, mult=mult)
				Nstep<-STEPS$Nstep
				step<-STEPS$step
			}
		}
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
		
		do.se.fit <- se.fit | ci.fit
		if(do.se.fit){
			if(is.null(var)){
				warning("the variance matrix of the parameters in NULL, unable to compute standard error and confidence interval of cumulative hazards are not computed")
				stdErrorPredictors <- NULL
				do.se.fit <- FALSE
			}
		}
		
		
		
		if(bhlink=="log"){
			if(dim(Y)[2] == 2){
				Predictors <- .computeCumulativeHazard_GA0B0AB(allparam=alltheparameters,
						Y=Y, X0=X0, X=X, Z=Z,
						step=step, Nstep=Nstep, 
						intTD=intTD, intweightsfunc=intweightsfunc,
						nT0basis=nT0basis,
						Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0= ibeta0, nX=nX, 
						ialpha=ialpha, 
						ibeta= ibeta, 
						nTbasis=nTbasis,
						Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH)
			} else {
				Predictors <- .computeCumulativeHazard_fromto_GA0B0AB(allparam=alltheparameters,
						Y=Y, X0=X0, X=X, Z=Z,
						step=step, Nstep=Nstep, 
						intTD=intTD, intweightsfunc=intweightsfunc,
						nT0basis=nT0basis,
						Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0= ibeta0, nX=nX, 
						ialpha=ialpha, 
						ibeta= ibeta, 
						nTbasis=nTbasis,
						Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH)
			}
			
			if (do.se.fit) {
				# dim(Y) managed in .computeStdErrorCumulativeHazard_GA0B0AB
				stdErrorPredictors <- .computeStdErrorCumulativeHazard_GA0B0AB(allparam=alltheparameters,
						var=object$var,
						Y=Y, X0=X0, X=X, Z=Z, 
						step=step, Nstep=Nstep, 
						intTD=intTD, intweightsfunc=intweightsfunc,
						intTD_base=intTD_base,
						nT0basis=nT0basis,
						Spline_t0=Spline_t0,
						Intercept_t0=Intercept_t0,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0=ibeta0, nX=nX,
						ialpha=ialpha, ibeta=ibeta,                             
						nTbasis=nTbasis,
						Spline_t=Spline_t,
						Intercept_t_NPH=Intercept_t_NPH,
						bhlink=bhlink)
			}
		} else {
			# bhlink == "identity"
			if(dim(Y)[2] == 2){
				Predictors <- .computeCumulativeHazard_GA0B0AB_bh(allparam=alltheparameters,
						Y=Y, X0=X0, X=X, Z=Z,
						step=step, Nstep=Nstep, 
						intTD=intTD, intweightsfunc=intweightsfunc,
						nT0basis=nT0basis,
						Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0= ibeta0, nX=nX, 
						ialpha=ialpha, 
						ibeta= ibeta, 
						nTbasis=nTbasis,
						Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH)
			} else {
				Predictors <- .computeCumulativeHazard_fromto_GA0B0AB_bh(allparam=alltheparameters,
						Y=Y, X0=X0, X=X, Z=Z,
						step=step, Nstep=Nstep, 
						intTD=intTD, intweightsfunc=intweightsfunc,
						nT0basis=nT0basis,
						Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0= ibeta0, nX=nX, 
						ialpha=ialpha, 
						ibeta= ibeta, 
						nTbasis=nTbasis,
						Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH)
			}
			
			if (do.se.fit) {
				# dim(Y) managed in .computeStdErrorCumulativeHazard_GA0B0AB_bh
				stdErrorPredictors <- .computeStdErrorCumulativeHazard_GA0B0AB_bh(allparam=alltheparameters,
						var=object$var,
						Y=Y, X0=X0, X=X, Z=Z, 
						step=step, Nstep=Nstep, 
						intTD=intTD, intweightsfunc=intweightsfunc,
						intTD_base=intTD_base,
						nT0basis=nT0basis,
						Spline_t0=Spline_t0,
						Intercept_t0=Intercept_t0,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0=ibeta0, nX=nX,
						ialpha=ialpha, ibeta=ibeta,                             
						nTbasis=nTbasis,
						Spline_t=Spline_t,
						Intercept_t_NPH=Intercept_t_NPH,
						bhlink=bhlink)
			}
		}
		
		if(ci.fit){
			qtnorm <- stats::qnorm(1 - (1 - level)/2)
			stdErrorLogPredictors <- stdErrorPredictors/Predictors 
			Predictors <- cbind(Predictors, Predictors * exp(qtnorm * stdErrorLogPredictors %o% c(-1, 1))) 
			colnames(Predictors) <- c("fit", "lwr", "upr")
			attr(Predictors, "level") <- level 
		}
		
		
		if (type=="surv"){
			Predictors <- pmax(exp(-Predictors ), .Machine$double.eps)
			if (se.fit) {
				stdErrorPredictors <- stdErrorPredictors * Predictors
			}
		}
		
		# build results object
		if(se.fit){
			pred <- list(fit=Predictors, se.fit=stdErrorPredictors)
		} else {
			pred <- Predictors
		}
		
		
	}
	
	pred
	
	
}




