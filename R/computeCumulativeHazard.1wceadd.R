# when dim(Y)==3, Y[,1] is beginTime, Y[,2] is endTime, Y[,3] is event
.computeCumulativeHazard_fromto_1wceadd<-function(allparam,
		Y, X0, X, Z, W,
		Id, FirstId, LastId, 
		step, Nstep,
		intTD=intTDft_NC, intweightsfunc=intweights_CAV_SIM,
		ialpha0, nX0,
		ibeta0, nX,
		ialpha, ibeta,
		ieta0,                             
		nTbasis,
		Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
		Intercept_t_NPH=rep(TRUE, nX), 
		ISpline_W =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
		Intercept_W=TRUE,
		debug=FALSE,  ...){
	# compute the cumulative hazard frm Y[,1] to Y[,2]
	# rate = wce(W,t)exp(X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
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
	# Y : object of class Surv with beginning and end of interval
	#
	# X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
	# X : log lineair but time dependante variable 
	# Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
	# expected_rate : expected rate at event time T
	# step : object of class "NCLagParam" or "GLMLagParam"
	# intTD : function to perform numerical integration 
	# intweightfunc : function to compute weightsfor numerical integration
	# nTbasis : number of time spline basis for NPH or NLL effects
	# nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
	# nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
	#  Spline_t, spline object for time dependant effects,  with evaluate() méthod
	# Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
	#  ... not used args
	# the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
	# returned value : the log liikelihood of the model
	
	
	if(is.null(Z)){
		nZ <- 0
	} else {
		nZ <- Z@nZ
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
	
	IS_W<- ISpline_W
	if(Intercept_W){
		eta0 <- allparam[ieta0]
	}
	else {
		eta0 <- c(0, allparam[ieta0])
	}
	IS_W <- ISpline_W * eta0
	if(nX + nZ) {
		stop("NPH effect not yet implemented", call.=TRUE)
		NPHterm <- intTD(rateTD_alphabeta_1addwce, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],
				step=step, Nstep=Nstep,
				intweightsfunc=intweightsfunc, 
				fromT=Y[,1], toT=Y[,2], fail=Y[,3], FirstId=FirstId,
				Zalphabeta=Zalphabeta,
				W = W, 
				Spline_t = Spline_t, Intercept_t=TRUE,
				ISpline_W = IS_W, Intercept_W=Intercept_W)
	} else {
		# no time dependent terms in the exp()
		# NPHTERM is the cumulative WCE effect between Tfrom and Tto
		# algebric formula
		
		wce2 <-  predictwce(object=integrate(IS_W), t=Y[,2], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
				FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE) 
		wce1 <- predictwce(object=integrate(IS_W), t=Y[,1], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
				FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
		NPHterm <- wce2 - wce1
	}
	
	# contribution of non time dependant variables
	if( nX0){
		ret <- exp(X0 %*% allparam[ialpha0]) * NPHterm 
	} else {
		ret <- NPHterm 
	}
	ret
}


