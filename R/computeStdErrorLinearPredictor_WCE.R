.computeStdErrorLinearPredictor_GA0B0ABE0<-function(allparam,
		var,
		Y, X0, X, Z, W, 
		Id, FirstId,
		nT0basis,
		Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
		ialpha0, nX0,
		ibeta0, nX,
		ialpha, ibeta,
		nTbasis,
		ieta0, iWbeg, iWend, nW, 
		Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
		listZSplineBasis,
		Intercept_t_NPH=rep(TRUE, nX), 
		ISpline_W =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
		Intercept_W=TRUE,
		debug=FALSE,  ...){
	# compute std error of a linearpredictor (log rate) if the model
	# rate = f(t)%*%gamma + X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
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
	# Y : object of class Surv but the matrix has 4 columns :
	# Y[,1] beginning(1) , fromT
	# Y[,2] end(2), toT,
	# Y[,3] status(3) fail
	# Y[,4] end of followup(4) 
	#     end of followup is assumed constant by Id
	# X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
	# X : log lineair but time dependante variable 
	# Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
	# W : Exposure variables used in Weighted Cumulative Exposure Models
	# Id : varibale indicating individuals Id, lines with the same Id are considered to be from the same individual
	# FirstId : all lines in FirstId[iT]:iT in the data comes from the same individual 
	# nT0basis : number of spline basis 
	#  Spline_t0, spline object for baseline hazard, with evaluate() méthod
	#  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
	# nTbasis : number of time spline basis for NPH or NLL effects
	# nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
	# nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
	#  Spline_t, spline object for time dependant effects,  with evaluate() méthod
	# Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
	# nW    : nb of WCE variables dim(W)=c(nobs, nW)
	# iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
	#  ISpline_W, list of nW spline object for WCE effects,  with evaluate() méthod
	#  ... not used args
	# the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
	# returned value : the log liikelihood of the model
	
	if ( debug) cat("  # computinfg stderr of the linear predictor: .computeStdErrorLlinearPredictor\n")
	
	nbeta0 <- length(ibeta0)
	nalpha0 <- length(ialpha0)
	
	# spline bases for baseline hazard
	colEndTime <- ifelse(ncol(Y)==2, 1, 2)
	if(is.null(Spline_t0)){
		DesignMatrix <-     X0
	}
	else {
		# baseline hazard at the end of the interval
		YT0 <- evaluate(Spline_t0, Y[,colEndTime], intercept=Intercept_t0)
		DesignMatrix <- cbind(YT0,
				X0)
	}
	
	
	# spline bases for each TD effect
	if(nX){
		# spline bases for each TD effect
		YT <- evaluate(Spline_t, Y[,colEndTime], intercept=TRUE)
		if(Intercept_t_NPH[1]){
			cbindDiagT <- YT
		}
		else {
			cbindDiagT <- YT[, -1]
		}
		
		cbindDiagX <- as.matrix(rep(1, nTbasis +  Intercept_t_NPH[1] -1 ))
		if( nX > 1 ){
			for( i in 2:nX){
				if(Intercept_t_NPH[i]){
					cbindDiagT <- cbind(cbindDiagT, YT)
				}
				else {
					cbindDiagT <- cbind(cbindDiagT, YT[, -1])
				}
				cbindDiagX <- cbindDiagX %sd% as.matrix(rep(1, nTbasis + Intercept_t_NPH[i] -1))
			}
		}
		# remove first basis if
		
		DesignMatrix <- cbind(DesignMatrix,
				(X %*%  t(cbindDiagX)) * cbindDiagT)
	}
	
	if(is.null(Z) & is.null(W)){
		nZ <- 0
		return(sqrt(diag(DesignMatrix %*% var %*% t(DesignMatrix))))
	} else {
		nlincoef <- nT0basis + nalpha0 + nbeta0
		dGdBeta <- diag(nlincoef)
		if(!is.null(Z)){
			#number of linear coef
			# spline bases for each TD effect
			YT <- evaluate(Spline_t, Y[,colEndTime], intercept=TRUE)
			nZ <- Z@nZ
			alpha <- allparam[ialpha]
			beta <- allparam[ibeta]
			for(i in 1:nZ){
				thenalpha <- getNBases(listZSplineBasis[[i]])-1
				cbindDiag1 <- duplicMat(diag(thenalpha), nTbasis)
				cbindDiag2 <- duplicMat(diag(nTbasis), thenalpha)
				thealpha <- alpha[Z@index[i,1]:Z@index[i,2]]
				thebeta <- beta[(i-1)*(nTbasis-1)+1:(nTbasis-1)]
				BdiagA <- duplicSumDirect(c(1, thebeta), thenalpha)
				BdiagB <- StackDiag( thealpha, nTbasis)[,-1]
				if(i == 1){
					GA <- BdiagA
					GB <- BdiagB
				}
				else {
					GA <- GA %sd% BdiagA
					GA <- GB %sd% BdiagB
				}
				DesignMatrix <- cbind(DesignMatrix,
						((Z@DM)[,Z@index[i,1]:Z@index[i,2]] %*% cbindDiag1) * (YT %*% cbindDiag2)
				)
				dGdBeta <- dGdBeta %sd% cbind(GA , GB)
			}
		}
		if(!is.null(W)){
			#number of linear coef
			# spline bases for each WCE effect
			
			################################################################################
			#### A COMPLETER
			################################################################################
		}
		
		DesignMatrix <- DesignMatrix %*% dGdBeta
		
		return(sqrt(diag(DesignMatrix %*% var %*% t(DesignMatrix))))
		
	}
	
}


