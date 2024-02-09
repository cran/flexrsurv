.computeLinearPredictor_fromto_1wceadd<-function(allparam,
		Y, X0, X, Z, W, 
		Id, FirstId, LastId, 
		ialpha0, nX0,
		ibeta0, nX,
		ialpha, ibeta,
		ieta0,
		nTbasis,
		Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
		Intercept_t_NPH=rep(TRUE, nX),
		wcelink=c("log", "identity"),
		ISpline_W =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
		Intercept_W=TRUE,
		debug=FALSE,  ...){
	# compute linearpredictor (log rate) if the model
	# rate = link(wce(W,t))exp(X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
	#################################################################################################################
	#################################################################################################################
	#  the coef of the first t-basis is constraint to 1 for nat-spline, and n-sum(other beta) if bs using expand() method
	#################################################################################################################
	#################################################################################################################
	#################################################################################################################
	# coef=
	# allparam c(eta0, alpha0, beta0, beta, alpha, brass0, balpha0  )
	#      ; vector of all coefs
	# with
	# eta0  : vector of all the coef for the WCE effects
	# alpha0 ; vector of all coefs for non time dependant variables (may contain non-loglinear terms such as spline)
	# beta0 ; matrix of all coefs for log-linear but  time dependant variables  X%*%beta0(t)
	# beta  : matrix of coefs for beta(t) nTbasis * nTDvars for NLG and NPH
	# alpha : vector of coef for alpha(z) for NLG and NPH
	# eta = allparam[1:neta0]
	# alpha0= allparam[ialpha0]
	# beta0= matrix(allparam[ibeta0], ncol=nX, nrow=nTbasis)
	# alpha= diag(allparam[ialpha])
	# beta= expand(matrix(allparam[ibeta], ncol=Z@nZ, nrow=nTbasis-1))
	# beta does not contains coef for the first t-basis
	
	#################################################################################################################
	# Y : object of class Surv (with ncol=2 or 3) Y[,ncol-1] is the time at which the predictors are computed
	# X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
	# X : log lineair but time dependante variable 
	# Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
	# nTbasis : number of time spline basis for NPH or NLL effects
	# nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
	# nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
	#  Spline_t, spline object for time dependant effects,  with evaluate() method
	# Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
	#  ... not used args
	# the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
	# returned value : the log liikelihood of the model
	
	wcelink  <- match.arg(wcelink)       # type baseline hazard
	
	if(is.null(Z)){
		nZ <- 0
	} else {
		nZ <- Z@nZ
	}
	
	
	IS_W<- ISpline_W
	if(Intercept_W){
		eta0 <- allparam[ieta0]
	}
	else {
		eta0 <- c(0, allparam[ieta0])
	}
	IS_W <- ISpline_W * eta0
	
	
	# contribution of non time dependant variables
	if( nX0){
		PHterm <-X0 %*% allparam[ialpha0]
	} else {
		PHterm <- 0.0
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
	
	# spline bases for baseline hazard
	colEndTime <- ifelse(ncol(Y)==2, 1, 2)
	# WCE at end of interval
	# eta0 = NULL because IS_W = ISpline_W * eta0
	WCEevent <- predictwce(object=IS_W, t=Y[,2], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
			FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
	
	if(nX + nZ){
		# spline bases for each TD effect
		YT <- evaluate(Spline_t, Y[,colEndTime], intercept=TRUE)
		linpred <- PHterm + apply(YT * Zalphabeta, 1, sum)
	} else {
		linpred <- PHterm 
	}
	if(wcelink == "log"){
		linpred <- as.vector(linpred + WCEevent) 
	} else {
		linpred <- cbind(linpred ,  WCEevent)
		dimnames(linpred)[[2]] <- c("linpred", "WCE")
	}
	
	linpred
}







