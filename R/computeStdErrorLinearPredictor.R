.computeStdErrorLinearPredictor_GA0B0AB<-function(allparam,
		var,
		Y, X0, X, Z, 
		nT0basis,
		Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
		ialpha0, nX0,
		ibeta0, nX,
		ialpha, ibeta,
		nTbasis,
		Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
		listZSplineBasis,
		Intercept_t_NPH=rep(TRUE, nX), 
		bhlink=c("log", "identity"),
		debug=FALSE,  ...){
	# compute jacobian matrix of the excess hazard of the relative survival model
	# rate = invlink(f(t)%*%gamma) exp(X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
	# if bhlink = log : return the gradient of f(t)%*%gamma + X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) )
	# if bhlink = identity : return the jacobian of the function F( allparam) = c(baseline = f(t)%*%gamma, linpred = X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) )
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
	# Y : object of class Surv (with ncol=2 or more)
	#                the time at which the predictors are computed is Y[,1] if ncol=2, Y[,2] if ncol>2
	#
	# X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
	# X : log lineair but time dependante variable 
	# Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
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
	
	if ( debug) cat("  # computinf stderr of the linear predictor: .computeStdErrorLlinearPredictor\n")
	
	bhlink  <- match.arg(bhlink)       # type baseline hazard
	
	if(dim(Y)[2] == 2){
		gr <- gr_link_flexrsurv_GA0B0AB(allparam=allparam,
				Y=Y, X0=X0, X=X, Z=Z, 
				nT0basis=nT0basis,
				Spline_t0=Spline_t0,
				Intercept_t0=Intercept_t0,
				ialpha0=ialpha0, nX0=nX0,
				ibeta0=ibeta0, nX=nX,
				ialpha=ialpha, ibeta=ibeta,                             
				nTbasis=nTbasis,
				Spline_t=Spline_t,
				Intercept_t_NPH=Intercept_t_NPH,
				debug=debug)
	}
	else {
		gr <- gr_link_flexrsurv_fromto_GA0B0AB(allparam=allparam,
				Y=Y, X0=X0, X=X, Z=Z, 
				nT0basis=nT0basis,
				Spline_t0=Spline_t0,
				Intercept_t0=Intercept_t0,
				ialpha0=ialpha0, nX0=nX0,
				ibeta0=ibeta0, nX=nX,
				ialpha=ialpha, ibeta=ibeta,                             
				nTbasis=nTbasis,
				Spline_t=Spline_t,
				Intercept_t_NPH=Intercept_t_NPH,
				debug=debug)
	}
	
	if(bhlink == "log"){
		varerr <- apply(gr * tcrossprod(gr, var), 1, sum)
		stderr <- sqrt(varerr)
	} else {
		# varerr is a nobs X 3 matrix with for each obs,
		#         var(bh(obs), linpred(obs) ) =  varerr[,1], varerr[,2]
		#                                        varerr[,2], varerr[,3]
		if(is.null(Spline_t0)){
			ngamma0 <- 0
			ibh <- NULL
			YT0Gamma0 <- 0.0
			varerr <- apply(gr * tcrossprod(gr, var), 1, sum)
		}
		else {
			ngamma0 <- getNBases(Spline_t0)-(1-Intercept_t0)
			ibh <-1:ngamma0
			ilinpred <- (ngamma0+1):length(allparam)
			varerr <- apply(gr[,ibh] * tcrossprod(gr[,ibh], var), 1, sum)
			varerr <- cbind(varerr, apply(gr[,ibh] * tcrossprod(gr[,ilinpred], var), 1, sum))
			varerr <- cbind(varerr, apply(gr[,ilinpred] * tcrossprod(gr[,ilinpred], var), 1, sum))
			names(varerr) <- c("baseline", "coverr", "linpred")
		}
		stderr <- sqrt(varerr[,c(1,3)])
	}
	
	attr(stderr, "varerr") <- varerr
	
	return(stderr)
}




