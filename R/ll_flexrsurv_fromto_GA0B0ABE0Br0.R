ll_flexrsurv_fromto_GA0B0ABE0Br0<-function(allparam, Y, X0, X, Z, W, BX0,
		Id, FirstId,LastId=NULL,
		isEnter, isEnd, 
		expected_rate, expected_logit_end, expected_logit_enter,
		weights=NULL,
		step, Nstep,
		intTD=intTD_NC, intweightsfunc=intweights_CAV_SIM,
		nT0basis,
		Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
		ialpha0, nX0,
		ibeta0, nX,
		ialpha, ibeta,
		nTbasis,
		ieta0, iWbeg, iWend, nW, 
		Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
		Intercept_t_NPH=rep(TRUE, nX), 
		ISpline_W =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
		Intercept_W=TRUE,
		nBbasis,
		Spline_B, Intercept_B=TRUE,
		ibrass0, nbrass0, 
		ibalpha0, nBX0,
		debug=FALSE,  ...){
	# compute log likelihood of the relative survival model
	# excess rate = exp( f(t)%*%gamma + X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) + sum ( wce(Wi , eta0i)(t)) ))
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
	# eta0 = allparam[ieta0]
	# brass0 = allparam[ibrass0]
	# balpha0 = allparam[ibalpha0]
	
	
	# corection of lifetable according to generalized brass method
	# Cohort-independent generalized Brass model in an age-cohort table
	# stratified brass model according to fixed effects BX0 (one brass function per combination)
	# rate = brass0(expected-rate, expected_logit)*exp(BX0 balpha0) + exp(gamma0(t) + time-independent effect(LL + NLL)(X0) + NPH(X) + NPHNLL(Z) + WCE(W))
	# brass0 : BRASS model wiht parameter Spline_B
	# logit(F) = evaluate(Spline_B, logit(F_pop), brass0) * exp(Balpha %*% BX0)
	# HCum(t_1, t_2) = log(1 + exp(evaluate(Spline_B, logit(F_pop(t_2)), brass0)) -  log(1 + exp(evaluate(Spline_B, logit(F_pop(t_1)), brass0))
	# rate(t_1) = rate_ref * (1 + exp(-logit(F_pop(t)))/(1 + exp(evaluate(Spline_B, logit(F_pop(t)), brass0)))*
	#                        evaluate(deriv(Spline_B), logit(F_pop(t)), brass0)
	# expected_logit_end =  logit(F_pop(t_2))
	# expected_logit_enter =  logit(F_pop(t_1))
	# brass0 = allparam[ibrass0]
	# Spline_B : object of class "AnySplineBasis" (suitable for Brass model) with method deriv() and evaluate()
	# isEnter = 1 ifoinlyf  the row corresponds to entry of followup AND EnterT > 0 (only used to compute modified_cumrate at enter)
	# isEnd   = 1 ifoinlyf  the row corresponds to end of followup
	#
	# parameters for excess rate
	#################################################################################################################
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
	# expected_rate : expected rate at event time T
	# expected_logit_end : logit of the expected survival at the end of the followup
	# expected_logit_enter : logit of the expected survival at the beginning of the followup 
	# weights : vector of weights  : LL = sum_i w_i ll_i
	# step : object of class "NCLagParam" or "GLMLagParam"
	# Nstep : number of lag for each observation
	# intTD : function to perform numerical integration 
	# intweightfunc : function to compute weightsfor numerical integration
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
	#  ISpline_W, list of nW spline object for WCE effects,  with evaluate() method
	#             ISpline is already integreted 
	#  ... not used args
	# the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
	# returned value : the log liikelihood of the model
	
	
	if(is.null(Z)){
		nZ <- 0
		Zalphabeta <- NULL
	} else {
		nZ <- Z@nZ
	}
	
	# LastId
	if(is.null(LastId)){
		first <- unique(FirstId)
		nline <- c(first[-1],length(FirstId)+1)-first
		LastId <- FirstId+rep(nline, nline)-1
	}
	
	
	if(is.null(Spline_t0)){
		YT0Gamma0 <- 0.0
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
							# no log basis for NPH and NPHNLL effects
							expand=!Intercept_t_NPH,
							value=0))
		}
	}
	if(nW){
		IS_W<- ISpline_W
		eta0 <- allparam[ieta0]
		for(iW in 1:nW){
			if(Intercept_W[[iW]]){
				IS_W[[iW]] <- ISpline_W[[iW]] * eta0[iWbeg[iW]:iWend[iW]]
			}
			else {
				IS_W[[iW]]<- ISpline_W[[iW]] * c(0, eta0[iWbeg[iW]:iWend[iW]])
			}
		}
		if(nX + nZ) {
			NPHterm <- intTD(rateTD_gamma0alphabetaeta0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],
					step=step, Nstep=Nstep,
					intweightsfunc=intweightsfunc, 
					fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
					gamma0=allparam[igamma0], Zalphabeta=Zalphabeta,
					nW = nW, W = W, eta0=allparam[ieta0], iWbeg=iWbeg, iWend=iWend,
					Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
					Spline_t = Spline_t, Intercept_t=TRUE,
					ISpline_W = IS_W, Intercept_W=Intercept_W)
		} else {
			NPHterm <- intTD(rateTD_gamma0eta0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
					step=step, Nstep=Nstep,
					intweightsfunc=intweightsfunc, 
					fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
					nW = nW, W = W, eta0=allparam[ieta0], iWbeg=iWbeg, iWend=iWend,
					Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
					ISpline_W = IS_W, Intercept_W=Intercept_W)
		}
	}
	else {
		# no WCE effect, same NPH term than ll_flexrsurv_fromto_GA0B0ABE0Br0
		if(nX + nZ) {
			NPHterm <- intTD(rateTD_gamma0alphabeta, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
					step=step, Nstep=Nstep,
					intweightsfunc=intweightsfunc, 
					gamma0=allparam[igamma0], Zalphabeta=Zalphabeta, 
					Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
					Spline_t = Spline_t, Intercept_t=TRUE)
		} else {
			NPHterm <- intTD(rateTD_gamma0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
					step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
					gamma0=allparam[igamma0],
					Spline_t0=Spt0g, Intercept_t0=Intercept_t0)
		}
	}
	################################################################################
	################################################################################
	################################################################################
	################################################################################
	################################################################################
	################################################################################
	################################################################################
	################################################################################
	################################################################################
	################################################################################
	################################################################################
	#****** eventterm : on peut peut-être faire un vecteur de lenght = nombre d'individu (correspondant au fin de suivi des individus)
	
	if(nW){
		# eta0 = NULL because IS_W = ISpline_W * eta0
		WCEcontrib <- weighted_cummulative_exposure(Increment=W, fromT=Y[,1], toT=, Y[,2], FirstId=FirstId, LastId=LastId,
				theT=Y[,4], tId=LastId,
				eta0=NULL, iWbeg=iWbeg, iWend=iWend, ISpline_W = IS_W, Intercept_W=Intercept_W)
		
	} else {
		WCEcontrib <- NULL
	}
	
# Brass model
	if(is.null(Spline_B)){
		if( nBX0){
			BPHterm <-exp(BX0 %*% allparam[ibalpha0])
			modified_rate <-  expected_rate * BPHterm
#      modified_cumrate <- ifelse(isEnd, log(1 + exp( expected_logit_end)) * BPHterm, 0) -
#        ifelse(isEnter, log(1 + exp(expected_logit_enter)) * BPHterm , 0)
			modified_cumrate <- log((1 + exp( expected_logit_end))/(1 + exp(expected_logit_enter))) * BPHterm  
		}
		else {
			modified_rate <-  expected_rate 
			modified_cumrate <- log((1 + exp( expected_logit_end))/(1 + exp(expected_logit_enter)))   
#      modified_cumrate <- ifelse(isEnd, log(1 + exp( expected_logit_end )), 0) -
#        ifelse(isEnter, log(1 + exp(expected_logit_enter) ), 0) 
		}
	}
	else {
		brass0 <- allparam[ibrass0]
		SB <- Spline_B * brass0
		# contribution of non time dependant variables
		if( nBX0){
			BPHterm <-exp(BX0 %*% allparam[ibalpha0])
			evalbrass <- exp(predict(SB, expected_logit_end)) 
			evalderivbrass <- predict(deriv(SB), expected_logit_end) 
			modified_rate <-  expected_rate * (1 + exp(-expected_logit_end))/(1+ 1/evalbrass) * evalderivbrass * BPHterm
			# modified cumrate is computed once for each individual (from t_enter to t_end of folowup)
			modified_cumrate <- log((1 + exp( evalbrass ))/(1 + exp(predict(SB, expected_logit_enter)))) * BPHterm 
		} else {
			evalbrass <- exp(predict(SB, expected_logit_end))
			evalderivbrass <- predict(deriv(SB), expected_logit_end) 
			modified_rate <-  expected_rate * (1 + exp(-expected_logit_end))/(1+ 1/evalbrass) * evalderivbrass
			# modified cumrate is computed once for each individual (from t_enter to t_end of folowup)
			# isEnter == 1 only when first row AND Tenter > 0
			modified_cumrate <- log((1 + exp( evalbrass ))/(1 + exp(predict(SB, expected_logit_enter))))
		}
	}
	
	# spline bases for each TD effect
	if(nX + nZ){
		# spline bases for each TD effect at the end of the interval
		YT <- evaluate(Spline_t, Y[,2], intercept=TRUE)
		if(nW){
			eventterm <- ifelse(Y[,3] ,
					log( modified_rate + PHterm * exp(YT0Gamma0 + apply(YT * Zalphabeta, 1, sum) + apply(WCEcontrib, 1, sum)) + modified_rate ),
					0)
		}
		else  {
			eventterm <- ifelse(Y[,3] ,
					log( modified_rate + PHterm * exp(YT0Gamma0 + apply(YT * Zalphabeta, 1, sum)) ),
					0)
		}
	} else {
		if(nW){
			eventterm <- ifelse(Y[,3] , 
					log( modified_rate + PHterm * exp(YT0Gamma0 + apply(WCEcontrib, 1, sum)) ), 
					0)
		}
		else {
			eventterm <- ifelse(Y[,3] , 
					log( modified_rate + PHterm * exp(YT0Gamma0) ), 
					0)
		}
	}
	
	if (!is.null(weights)) {
		if( nX0){
			ret <- crossprod(eventterm + modified_cumrate - PHterm * NPHterm , weights)
		} else {
			ret <- crossprod(eventterm + modified_cumrate - NPHterm , weights)
		}
	}
	else {
		if( nX0){
			ret <- sum( eventterm + modified_cumrate - PHterm * NPHterm )
		} else {
			ret <- sum( eventterm + modified_cumrate - NPHterm )
		}
	}
	
	if ( debug) {
		attr(ret, "eventterm") <- eventterm
		attr(ret, "PHterm") <- PHterm
		attr(ret, "NPHterm") <- NPHterm
		attr(ret, "WCEcontrib") <- WCEcontrib
		attr(ret, "modified_rate") <- modified_rate
		attr(ret, "modified_cumrate") <- modified_cumrate
		
		if ( debug > 1000) cat("fin ll_flexrsurv_GA0B0ABE0Br0 **", ret, "++ \n")
	}
	ret
}






