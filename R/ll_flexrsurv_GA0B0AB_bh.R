ll_flexrsurv_GA0B0AB_bh<-function(allparam, Y, X0, X, Z, 
		expected_rate,  
		weights=NULL,
		step, Nstep,
		intTD=intTD_NC, intweightsfunc=intweights_CAV_SIM,
		nT0basis,
		Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
		ialpha0, nX0,
		ibeta0, nX,
		ialpha, ibeta,                             
		nTbasis,
		Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
		Intercept_t_NPH=rep(TRUE, nX), 
		debug=FALSE,  ...){
	# compute log likelihood of the relative survival model
	# rate = (f(t)%*%gamma )*exp( X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
	#################################################################################################################
	#################################################################################################################
	#  the coef of the first t-basis is constraint to 1 for nat-spline, and n-sum(other beta) if bs using expand() method
	#################################################################################################################
	#################################################################################################################
	#################################################################################################################
	# allparam ; vector of all coefs
	# gamma0 = allparam[1:nY0basis]
	# alpha0= allparam[ialpha0]
	# beta0= matrix(allparam[ibeta0], ncol=nX, nrow=nTbasis)
	# alpha= diag(allparam[ialpha])
	# beta= expand(matrix(allparam[ibeta], ncol=Z@nZ, nrow=nTbasis-1))
	# beta does not contains coef for the first t-basis
	#################################################################################################################
	# Y : object of class Surv
	# X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
	# X : log lineair but time dependante variable 
	# Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
	# expected_rate : expected rate at event time T
	# weights : vector of weights  : LL = sum_i w_i ll_i
	# step : object of class "NCLagParam" or "GLMLagParam"
	# intTD : function to perform numerical integration 
	# intweightfunc : function to compute weightsfor numerical integration
	# nT0basis : number of spline basis 
	#  Spline_t0, spline object for baseline hazard, with evaluate() method
	#  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
	# nTbasis : number of time spline basis for NPH or NLL effects
	# nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
	# nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
	#  Spline_t, spline object for time dependant effects,  with evaluate() method
	# Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
	#  ... not used args
	# the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
	# returned value : the log liikelihood of the model
	
	if ( debug) cat("  # omputinfg the log likelihood: ll_flexrsurv_GOA0B0AB\n")
	
	if(is.null(Z)){
		nZ <- 0
	} else {
		nZ <- Z@nZ
	}
	if(is.null(Spline_t0)){
		YT0Gamma0 <- 1.0
		Spt0g <- NULL
		igamma0 <- NULL
		
	}
	else {
		igamma0 <- 1:nT0basis
		if(Intercept_t0){
			tmpgamma0 <- allparam[igamma0]
		}
		else {
			tmpgamma0 <- c(0, allparam[igamma0])
		}
		# baseline hazard at the end of the interval
		
		Spt0g <- Spline_t0*tmpgamma0
		YT0Gamma0 <- predictSpline(Spt0g, Y[,2])
	}
	
	# contribution of non time dependant variables
	if( nX0){
		PHterm <-exp(X0 %*% allparam[ialpha0])
	} else {
		PHterm <- 1
	}
	# contribution of time d?pendant effect
	# parenthesis are important for efficiency
	if(nZ) {
		# add a row of one for the first T-basis 
		Beta <- t(ExpandAllCoefBasis(allparam[ibeta], ncol=nZ,  value=1))
		# parenthesis important for speed ?
		Zalphabeta <- Z@DM %*%( diag(allparam[ialpha]) %*% Z@signature  %*% Beta )
		if(nX) {
			# add a row of 0 for the first T-basis when !Intercept_T_NPH
			Zalphabeta <- Zalphabeta + X %*% t(ExpandCoefBasis(allparam[ibeta0],
							ncol=nX,
							splinebasis=Spline_t,
							expand=!Intercept_t_NPH,
							value=0))
		}
	} else {
		if(nX) {
			Zalphabeta <- X %*% t(ExpandCoefBasis(allparam[ibeta0],
							ncol=nX,
							splinebasis=Spline_t,
							expand=!Intercept_t_NPH,
							value=0))
		}
	}
	
	if(nX + nZ) {
		NPHterm <- intTD(rateTD_bh_alphabeta, intTo=Y[,1], intToStatus=Y[,2],
				step=step, Nstep=Nstep,
				intweightsfunc=intweightsfunc, 
				gamma0=allparam[igamma0], Zalphabeta=Zalphabeta, 
				Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
				Spline_t = Spline_t, Intercept_t=TRUE)
	} else {
#    NPHterm <- intTD(rateTD_gamma0_bh, intTo=Y[,1], intToStatus=Y[,2], 
#                     step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
#                     gamma0=allparam[igamma0],
#                     Spline_t0=Spt0g, Intercept_t0=Intercept_t0)
		
#   NPHterm <- integrate(Spline_t0, Y[,1], intercep=Intercept_t0) %*% allparam[igamma0]
		if(is.null(Spline_t0)){
			NPHterm <- 1.0
		}
		else {
			NPHterm <- predict(integrate(Spt0g), Y[,1], intercep=Intercept_t0)
		}
	}
	# spline bases for each TD effect
	if(nX + nZ){
		# spline bases for each TD effect
		YT <- evaluate(Spline_t, Y[,1], intercept=TRUE)
		eventterm <- ifelse(Y[,2] ,
				log( PHterm * YT0Gamma0 * exp(apply(YT * Zalphabeta, 1, sum)) + expected_rate ),
				0)
	} else {
		eventterm <- ifelse(Y[,2] , 
				log( PHterm * YT0Gamma0 + expected_rate ), 
				0)
	}
	
	if (!is.null(weights)) {
		if( nX0){
			ret <- crossprod(eventterm - PHterm * NPHterm , weights)
		} else {
			ret <- crossprod(eventterm - NPHterm , weights)
		}
	}
	else {
		if( nX0){
			ret <- sum( eventterm - PHterm * NPHterm )
		} else {
			ret <- sum( eventterm - NPHterm )
		}
	}
	
	if ( debug) {
		attr(ret, "eventterm") <- eventterm
		attr(ret, "PHterm") <- PHterm
		attr(ret, "NPHterm") <- NPHterm
	}
	ret
}






