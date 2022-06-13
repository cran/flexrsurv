# idem rateTD.R but baseline hasard is gamma0(t), not exp(gamma0(t))
rateTD_bh_beta0alphabeta<- function(T, iT, gamma0, Zbeta0, Zalphabeta, 
		Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
		Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE, ...){
	# compute the contribution of the time dependent variables to the rate 
	# of relative survival model for patient iT with Zalphabeta[iT, ]
	# at a vector of T (useful to compute numerical integration 
	# spline bases for baseline hazard
	# Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
	#           : thus no nead to multiply each spline coordinate by its coef
	
	# spline bases for baseline hazard
	if(is.null(Spline_t0)){
		YT0Gamma0 <- 1.0
	}
	else {
		YT0Gamma0 <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
	}
	# spline bases for each TD effect
#  YT  <- bs(T, knots=Knots_t, intercept=Intercept_t, degree=degree_t, Boundary.knots =  Boundary.knots_t)
	if(!is.null(Zbeta0)) {
		if(!is.null(Zalphabeta)) {
			YT0Gamma0 * 
					exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] +
									fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
		}
		else {
			YT0Gamma0 *
					exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] )
		}
	}
	else if(!is.null(Zalphabeta)) {
		YT0Gamma0 *
				exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
	}
}

rateTD_bh_alphabeta<- function(T, iT, gamma0, Zalphabeta, 
		Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
		Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE, ...){
	# compute the contribution of the time dependent variables to the rate 
	# of relative survival model for patient iT with Zalphabeta[iT, ]
	# at a vector of T (useful to compute numerical integration 
	# Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
	#           : thus no nead to multiply each spline coordinate by its coef
	
	# spline bases for baseline hazard
	if(is.null(Spline_t0)){
		YT0Gamma0 <- 1.0
	}
	else {
		YT0Gamma0 <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
	}
	# spline bases for each TD effect
	YT  <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)
	
	# returned value 
	YT0Gamma0 * exp( YT %*% Zalphabeta[iT,])
	
}

rateTD_gamma0_bh<- function(T, iT, gamma0, 
		Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE, ...){
	# compute the contribution of the baseline hazard rate to the rate 
	# of relative survival model for patient iT with Zalphabeta[iT, ]
	# at a vector of T (useful to compute numerical integration 
	# Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
	#           : thus no nead to multiply each spline coordinate by its coef
	
	# spline bases for baseline hazard
	if(is.null(Spline_t0)){
		YT0Gamma0 <- rep(1.0, length(T))
	}
	else {
		YT0Gamma0 <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
	}
	# returned value
	
	YT0Gamma0
	
}

rateTD_bh<- function(T, iT, ... ){
	# compute the contribution of the baseline hazard rate when there are no time dependent effect
	# and when the baseline hazard is null
	
	# returned value
	
	rep(0.0, length(T))
	
}


######################################################################
###
###  WCE

# computes the contribution of time dependent termes in the rate (baseline, NPH, NPHNLL and WCE effects 
rateTD_gamma0alphabetaeta0_bh<- function(T, iT,
		fromT, toT, fail, FirstId,
		gamma0, Zalphabeta,
		nW, W, eta0, iWbeg, iWend,
		Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
		Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE,
		ISpline_W, Intercept_W=rep(TRUE,nW), ...){
	# compute the contribution of the time dependent variables to the rate 
	# of relative survival model for line iT with Zalphabeta[iT, ]
	# at a vector of T (useful to compute numerical integration
	# fromT : begining of the time intervals
	# toT   : end of the time intervals
	# fail  : fail = 1 if event, 0 if censored
	# all T basis for the NPH/td effects are the same (Spline_t)
	# Zalphabeta = X %*% beta0 + f(Z,alpha) %*% beta 
	# FirstId : all lines in FirstId[iT]:iT of fromT, toT, fail and Zalphabeta comes from the same individual 
	# nW number of cols in W (number of WCE effects
	# W matrix of exposure INCREMENT variables W[FirstId[iT]:iT, k] is the vectore of exposure increment x_il - x_(i-1)l for patient l, expo variable k    
	# eta0 : vector all the coef of WCE
	# iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
	# ISpline_W : list of the nW integrated splines parameters for the WCE effects multiplied by eta0
	#           : thus no nead to multiply each spline coordinate by its coef
	# Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
	#           : thus no nead to multiply each spline coordinate by its coef
	
#  print("rateTD_gamma0alphabetaeta0")
	# spline bases for baseline hazard
	if(is.null(Spline_t0)){
		YT0Gamma0 <- 1.0
	}
	else {
		YT0Gamma0 <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
	}
	# spline bases for each TD effect
	YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)
	
	
	WCE <- 0.0
	for(iW in 1:nW){
		for(iId in FirstId[iT]:iT){
			WCE <- WCE + W[iId, iW] * predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)  
		}
	}
	# returned value 
	exp(WCE + YT %*% Zalphabeta[iT,])*YT0Gamma0
	
}


# computes the contribution of time dependent termes in the rate (baseline, NPH, NPHNLL and WCE effects 
rateTD_gamma0alphabetaeta0_bh_bh<- function(T, iT,
		fromT, toT, fail, FirstId,
		gamma0, Zalphabeta,
		nW, W, eta0, iWbeg, iWend,
		Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
		Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE,
		ISpline_W, Intercept_W=rep(TRUE,nW), ...){
	# compute the contribution of the time dependent variables to the rate 
	# of relative survival model for line iT with Zalphabeta[iT, ]
	# at a vector of T (useful to compute numerical integration
	# fromT : begining of the time intervals
	# toT   : end of the time intervals
	# fail  : fail = 1 if event, 0 if censored
	# all T basis for the NPH/td effects are the same (Spline_t)
	# Zalphabeta = X %*% beta0 + f(Z,alpha) %*% beta 
	# FirstId : all lines in FirstId[iT]:iT of fromT, toT, fail and Zalphabeta comes from the same individual 
	# nW number of cols in W (number of WCE effects
	# W matrix of exposure INCREMENT variables W[FirstId[iT]:iT, k] is the vectore of exposure increment x_il - x_(i-1)l for patient l, expo variable k    
	# eta0 : vector all the coef of WCE
	# iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
	# ISpline_W : list of the nW integrated splines parameters for the WCE effects multiplied by eta0
	#           : thus no nead to multiply each spline coordinate by its coef
	# Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
	#           : thus no nead to multiply each spline coordinate by its coef
	
#  print("rateTD_gamma0alphabetaeta0")
	# spline bases for baseline hazard
	if(is.null(Spline_t0)){
		WCE <- 0.0
	}
	else {
		WCE <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
	}
	# spline bases for each TD effect
	YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)
	
	
	for(iW in 1:nW){
		for(iId in FirstId[iT]:iT){
			WCE <- WCE + W[iId, iW] * predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)  
		}
	}
	# returned value 
	WCE*exp(YT %*% Zalphabeta[iT,])
	
}


# WCE additif (link=identity, pas de BH (gamma0)
# computes the contribution of time dependent termes in the rate (NPH, NPHNLL and WCE effects 
rateTD_alphabetaeta0_bh<- function(T, iT,
		fromT, toT, fail, FirstId,
		Zalphabeta,
		nW, W, eta0, iWbeg, iWend,
		Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE,
		ISpline_W, Intercept_W=rep(TRUE,nW), ...){
	# compute the contribution of the time dependent variables to the rate 
	# of relative survival model for line iT with Zalphabeta[iT, ]
	# at a vector of T (useful to compute numerical integration
	# fromT : begining of the time intervals
	# toT   : end of the time intervals
	# fail  : fail = 1 if event, 0 if censored
	# all T basis for the NPH/td effects are the same (Spline_t)
	# Zalphabeta = X %*% beta0 + f(Z,alpha) %*% beta 
	# FirstId : all lines in FirstId[iT]:iT of fromT, toT, fail and Zalphabeta comes from the same individual 
	# nW number of cols in W (number of WCE effects
	# W matrix of exposure INCREMENT variables W[FirstId[iT]:iT, k] is the vectore of exposure increment x_il - x_(i-1)l for patient l, expo variable k    
	# eta0 : vector all the coef of WCE
	# iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
	# ISpline_W : list of the nW integrated splines parameters for the WCE effects multiplied by eta0
	#           : thus no nead to multiply each spline coordinate by its coef
	
#  print("rateTD_gamma0alphabetaeta0")
	# spline bases for baseline hazard
	WCE <- rep(0, length(T))
	# spline bases for each TD effect
	YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)
	
	
	for(iW in 1:nW){
		for(iId in FirstId[iT]:iT){
			WCE <- WCE + W[iId, iW] * predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)  
		}
	}
	# returned value 
	WCE*exp(YT %*% Zalphabeta[iT,])
	
}



# computes the contribution of gamma0 and eta0 (baseline & WCE)
rateTD_gamma0eta0_bh<- function(T, iT,
		fromT, FirstId,
		gamma0, 
		nW, W, eta0, iWbeg, iWend,
		Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
		ISpline_W, Intercept_W=rep(TRUE,nW), ...){
	# compute the contribution of the time dependent variables to the rate 
	# of relative survival model for patient iT 
	# at a vector of T (useful to compute numerical integration
	# fromT : begining of the time intervals
	# toT   : end of the time intervals
	# fail  : fail = 1 if event, 0 if censored
	# FirstId : all lines in FirstId[iT]:iT of fromT, toT, fail and Zalphabeta comes from the same individual 
	# nW number of cols in W (number of WCE effects
	# W matrix of exposure INCREMENT variables W[FirstId[iT]:iT, k] is the vectore of exposure increment x_il - x_(i-1)l for patient l, expo variable k    
	# eta0 : vector all the coef of WCE
	# iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
	# ISpline_W : list of the nW integrated splines parameters for the WCE effects multiplied by eta0
	#           : thus no nead to multiply each spline coordinate by its coef
	# Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
	#           : thus no nead to multiply each spline coordinate by its coef
	
#  print("rateTD_gamma0eta0")
	
	# spline bases for baseline hazard
#  cat("************************************************************************\niT, T, first: ")
#  cat(c(iT, FirstId[iT], T))
#  cat("\n")
	
	if(is.null(Spline_t0)){
		YT0Gamma0 <- 1.0
	}
	else {
		YT0Gamma0 <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
	}
	WCE <- rep(0, length(T))
	for(iW in 1:nW){
		for(iId in FirstId[iT]:iT){
#  cat("iId, W[Iid, iW], fromT[iId]: ")
#  cat(c(iId, W[iId, iW], fromT[iId]))
#  cat("\n")
#  print(cbind(T, predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)))
			WCE <- WCE + W[iId, iW] * predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)  
		}
	}
	# returned value
	exp(WCE)*YT0Gamma0
	
}



# computes the contribution of gamma0 and eta0 (baseline & WCE)
rateTD_gamma0_bh_eta0_bh<- function(T, iT,
		fromT, FirstId,
		gamma0, 
		nW, W, eta0, iWbeg, iWend,
		Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
		ISpline_W, Intercept_W=rep(TRUE,nW), ...){
	# compute the contribution of the time dependent variables to the rate 
	# of relative survival model for patient iT 
	# at a vector of T (useful to compute numerical integration
	# fromT : begining of the time intervals
	# toT   : end of the time intervals
	# fail  : fail = 1 if event, 0 if censored
	# FirstId : all lines in FirstId[iT]:iT of fromT, toT, fail and Zalphabeta comes from the same individual 
	# nW number of cols in W (number of WCE effects
	# W matrix of exposure INCREMENT variables W[FirstId[iT]:iT, k] is the vectore of exposure increment x_il - x_(i-1)l for patient l, expo variable k    
	# eta0 : vector all the coef of WCE
	# iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
	# ISpline_W : list of the nW integrated splines parameters for the WCE effects multiplied by eta0
	#           : thus no nead to multiply each spline coordinate by its coef
	# Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
	#           : thus no nead to multiply each spline coordinate by its coef
	
#  print("rateTD_gamma0eta0")
	
	# spline bases for baseline hazard
#  cat("************************************************************************\niT, T, first: ")
#  cat(c(iT, FirstId[iT], T))
#  cat("\n")
	
	if(is.null(Spline_t0)){
		WCE <- 0.0
	}
	else {
		WCE <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
	}
	for(iW in 1:nW){
		for(iId in FirstId[iT]:iT){
#  cat("iId, W[Iid, iW], fromT[iId]: ")
#  cat(c(iId, W[iId, iW], fromT[iId]))
#  cat("\n")
#  print(cbind(T, predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)))
			WCE <- WCE + W[iId, iW] * predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)  
		}
	}
	# returned value
	WCE
	
}



# computes the contribution eta0 (WCE)
rateTD_eta0_bh<- function(T, iT,
		fromT, FirstId,
		nW, W, eta0, iWbeg, iWend,
		ISpline_W, Intercept_W=rep(TRUE,nW), ...){
	# compute the contribution of the time dependent variables to the rate 
	# of relative survival model for patient iT 
	# at a vector of T (useful to compute numerical integration
	# fromT : begining of the time intervals
	# toT   : end of the time intervals
	# fail  : fail = 1 if event, 0 if censored
	# FirstId : all lines in FirstId[iT]:iT of fromT, toT, fail and Zalphabeta comes from the same individual 
	# nW number of cols in W (number of WCE effects
	# W matrix of exposure INCREMENT variables W[FirstId[iT]:iT, k] is the vectore of exposure increment x_il - x_(i-1)l for patient l, expo variable k    
	# eta0 : vector all the coef of WCE
	# iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
	# ISpline_W : list of the nW integrated splines parameters for the WCE effects multiplied by eta0
	#           : thus no nead to multiply each spline coordinate by its coef
	
#  print("rateTD_eta0")
	
	
	WCE <- rep(0.0, length(T))
	for(iW in 1:nW){
		for(iId in FirstId[iT]:iT){
#  cat("iId, W[Iid, iW], fromT[iId]: ")
#  cat(c(iId, W[iId, iW], fromT[iId]))
#  cat("\n")
#  print(cbind(T, predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)))
			WCE <- WCE + W[iId, iW] * predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)  
		}
	}
	# returned value
	WCE
	
}




######################################################################
# same as rateTD* but witout the baseline


# idem rateTD.R but baseline hasard is gamma0(t), not exp(gamma0(t))
ratioTD_bh_beta0alphabeta<- function(T, iT, Zbeta0, Zalphabeta, 
		Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE, ...){
	# compute the contribution of the time dependent variables to the rate 
	# of relative survival model for patient iT with Zalphabeta[iT, ]
	# at a vector of T (useful to compute numerical integration 
	# spline bases for baseline hazard
	# spline bases for each TD effect
#  YT  <- bs(T, knots=Knots_t, intercept=Intercept_t, degree=degree_t, Boundary.knots =  Boundary.knots_t)
	if(!is.null(Zbeta0)) {
		if(!is.null(Zalphabeta)) {
			exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] +
							fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
		}
		else {
			exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] )
		}
	}
	else if(!is.null(Zalphabeta)) {
		exp(fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
	}
}

ratioTD_bh_alphabeta<- function(T, iT, Zalphabeta, 
		Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE, ...){
	# compute the contribution of the time dependent variables to the rate 
	# of relative survival model for patient iT with Zalphabeta[iT, ]
	# at a vector of T (useful to compute numerical integration 
	
	# spline bases for each TD effect
	YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)
	
	# returned value 
	exp( YT %*% Zalphabeta[iT,])
	
}

ratioTD_bh<- function(T, iT, ...){
	# compute the contribution of the baseline hazard rate to the rate 
	# of relative survival model for patient iT with Zalphabeta[iT, ]
	# at a vector of T (useful to compute numerical integration 
	
	# returned value
	
	return(rep(1, length(T)))
	
}





