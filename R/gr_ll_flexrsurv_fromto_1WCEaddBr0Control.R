gr_ll_flexrsurv_fromto_1WCEaddBr0Control<-function(allparam,
		Y, X0, X, Z, W,
		BX0,
		Id, FirstId, LastId,
		expected_rate,
		expected_logit_end,
		expected_logit_enter,
		weights=NULL,
		Ycontrol, BX0control, 
		weightscontrol=NULL,
		Idcontrol, FirstIdcontrol, LastIdcontrol,
		expected_ratecontrol,
		expected_logit_endcontrol,
		expected_logit_entercontrol,
		step, Nstep,
		intTD=intTDft_NC, intweightsfunc=intweights_CAV_SIM,
		intTD_base=intTDft_base_NC,
		intTD_WCEbase=intTDft_WCEbase_NC,
		ialpha0, nX0,
		ibeta0, nX,
		ialpha,
		ibeta,
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
		debug.gr=TRUE,  ...){
	# same as ll_flexrsurv_fromto_GA0B0ABE0Br0.R but with a control group
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
	# for control group
	# rate = brass0(expected-ratecontrol, expected_logitcontrol)*exp(BX0control balpha0) 
	# but for exposed
	# rate = brass0(expected-rate, expected_logit)*exp(BX0 balpha0) + exp(gamma0(t) + time-independent effect(LL + NLL)(X0) + NPH(X) + NPHNLL(Z) + WCE(W))
	# brass0 : BRASS model wiht parameter Spline_B
	# logit(F) = evaluate(Spline_B, logit(F_pop), brass0) * exp(Balpha %*% BX0)
	# HCum(t_1, t_2) = log(1 + exp(evaluate(Spline_B, logit(F_pop(t_2)), brass0)) -  log(1 + exp(evaluate(Spline_B, logit(F_pop(t_1)), brass0))
	# rate(t_1) = rate_ref * (1 + exp(-logit(F_pop(t)))/(1 + exp(evaluate(Spline_B, logit(F_pop(t)), brass0)))*
	#                        evaluate(deriv(Spline_B), logit(F_pop(t)), brass0)
	# expected_logit_end =  logit(F_pop(Y[,2]))
	# expected_logit_enter =  logit(F_pop(Y[,1]))
	# brass0 = allparam[ibrass0]
	# Spline_B : object of class "AnySplineBasis" (suitable for Brass model) with method deriv() and evaluate()
	#            IMPORTANT : the coef of the first basis is constraints to one and evaluate(deriv(spline_B), left_boundary_knots) == 1 for Brass transform 
	#
	# parameters for exposed group
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
	# BX0 : non-time dependante variable for the correction of life table (may contain spline bases expended for non-loglinear terms)
	# Id : varibale indicating individuals Id, lines with the same Id are considered to be from the same individual
	# FirstId : all lines in FirstId[iT]:iT in the data comes from the same individual 
	# LastId  : all lines in FirstId[iT]:LastId[iT] in the data comes from the same individual Id[iT] 
	# expected_rate : expected rate at event time T
	# expected_logit_end : logit of the expected survival at the end of the followup
	# expected_logit_enter : logit of the expected survival at the beginning of the followup 
	# weights : vector of weights  : LL = sum_i w_i ll_i
	# parameters for exposd population
	#################################################################################################################
	# parameters for exposed group
	#################################################################################################################
	# Ycontrol : object of class Surv but the matrix has 4 columns :
	# Ycontrol[,1] beginning(1) , fromT
	# Ycontrol[,2] end(2), toT,
	# Ycontrol[,3] status(3) fail
	# Ycontrol[,4] end of followup(4) 
	#     end of followup is assumed constant by Id
	# BX0control : non-time dependante variable for the correction of life table (may contain spline bases expended for non-loglinear terms)
	# Idcontrol : varibale indicating individuals Id, lines with the same Id are considered to be from the same individual
	# FirstIdcontrol : all lines in FirstId[iT]:iT in the data comes from the same individual 
	# LastIdcontrol  : all lines in FirstId[controliT]:LastIdcontrol[iT] in the data comes from the same individual Id[iT] 
	# expected_ratecontrol : expected rate at event time T
	# expected_logit_endcontrol : logit of the expected survival at the end of the followup
	# expected_logit_entercontrol : logit of the expected survival at the beginning of the followup 
	# weightscontrol : vector of weights  : LL = sum_i w_i ll_i
	#################################################################################################################
	# model parameters
	# step : object of class "NCLagParam" or "GLMLagParam"
	# Nstep : number of lag for each observation
	# intTD : function to perform numerical integration 
	# intweightfunc : function to compute weightsfor numerical integration
	# nTbasis : number of time spline basis for NPH or NLL effects
	# nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
	# nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
	#  Spline_t, spline object for time dependant effects,  with evaluate() method
	# Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
	# nW    : nb of WCE variables dim(W)=c(nobs, nW)
	# iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
	#  ISpline_W, list of nW spline object for WCE effects,  with evaluate() method
	#             ISpline is already integreted 
	#  ... not used args
	# the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
	# returned value : the log liikelihood of the model
	
	
#cat("gr ")
#print(allparam, digits=2)
#cat("+=+=+=++++++++++++++++++++++ gr ")
#print(allparam)
	
#  print("+++++++++++++++++++++++++++++++++++++++++++passage par le gradient")
	
	################################################################################
#  excess rate
	if(is.null(Z)){
		nZ <- 0
		Zalphabeta <- NULL
	} else {
		nZ <- Z@nZ
	}
	
	# contribution of non time dependant variables
	if( nX0){
		PHterm <-exp(X0 %*% allparam[ialpha0])
	} else {
		PHterm <- rep(1.0, dim(Y)[1])
	}
	# contribution of time d?pendant effect
	# parenthesis are important for efficiency
	if(nZ) {
		# add a row for the first basis
		tBeta <- t(ExpandAllCoefBasis(allparam[ibeta], ncol=nZ,  value=1))
		# Zalpha est la matrice des alpha(Z)
		# parenthesis important for speed ?
		Zalpha <- Z@DM %*%( diag(allparam[ialpha]) %*% Z@signature )
		Zalphabeta <- Zalpha  %*% tBeta 
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
		else {
			Zalphabeta <- NULL
		}
	}
	
	IS_W <- ISpline_W
	eta0 <- allparam[ieta0]
	if(Intercept_W){
		IS_W <- ISpline_W * eta0
	}
	else {
		IS_W<- ISpline_W * c(0, eta0)
	}
	IIS_W <- integrate(IS_W)
	IISpline_W <- integrate(ISpline_W)
	
	if(nX + nZ) {
#    stop("NPH effect not yet implemented", call.=TRUE)
		NPHterm <- intTD(rateTD_alphabeta_1addwce, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
				step=step, Nstep=Nstep,
				intweightsfunc=intweightsfunc, 
				fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
				Zalphabeta=Zalphabeta,
				W = W, 
				Spline_t = Spline_t, Intercept_t=TRUE,
				ISpline_W = IS_W, Intercept_W=Intercept_W)
		Intb <-  intTD_base(func=rateTD_alphabeta_1addwce,  intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
				Spline=Spline_t,
				step=step, Nstep=Nstep,
				intweightsfunc=intweightsfunc, 
				fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
				Zalphabeta=Zalphabeta, 
				W = W, 
				Spline_t = Spline_t, Intercept_t=TRUE,
				ISpline_W = IS_W, Intercept_W=Intercept_W,
				debug=debug.gr)
		indx_without_intercept <- 2:getNBases(Spline_t)
		
		gradCumWCE    <-  intTD_WCEbase(func=rateTD_alphabeta_1addwce,  intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
				Spline=ISpline_W, intercept=Intercept_W,
				step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
				Zalphabeta=Zalphabeta, 
				theW=W, fromT=Y[,1], toT=Y[,2], FirstId=FirstId,  LastId=LastId,
				W = W, 
				Spline_t = Spline_t, Intercept_t=TRUE,
				ISpline_W = IS_W, Intercept_W=Intercept_W,
				debug=debug.gr)
	} 
	else {
		# no time dependent terms in the exp()
		# NPHTERM is the cumulative WCE effect between Tfrom and Tto
		# algebric formula
		
		wce2 <-  predictwce(object=IIS_W, t=Y[,2], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
				FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE) 
		wce1 <- predictwce(object=IIS_W, t=Y[,1], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
				FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
		NPHterm <- wce2 - wce1
		#d_NPHTerm / d_eta0 = bases of IIS_W = integrated IS_W
		gradCumWCE <-  gradientwce(object=IISpline_W, t=Y[,2], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
				FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
		gr2 <- gradientwce(object=IISpline_W, t=Y[,1], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
				FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
		gradCumWCE <- gradCumWCE -gr2
		
		
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
	#***** 
	
	# WCE at end of interval
	# eta0 = NULL because IS_W = ISpline_W * eta0
	WCEevent <- predictwce(object=IS_W, t=Y[,2], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
			FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
	
	gradWCEevent <- gradientwce(object=ISpline_W, t=Y[,2], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
			FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
	
	
#print(summary(Y[,2]))
#print(summary(Y[unique(LastId),2]))
#print(summary((Y[,2]-Y[,1])[W>0]))
#  print(ISpline_W)
#print(ISpline_W@Matrices)
#print(IS_W)
#print(IS_W@Matrices)
	
#print(cbind(Y[,1], W, Y[,2], gradWCEevent)[Y[LastId,2]>10000,])
#print(cbind(Y[,1], W, Y[,2], gradCumWCE  )[Y[LastId,2]>10000,])
	
	
#  stop("ISPLINEW", call.=TRUE)
	################################################################################
	# control group
# only Brass model
	if(!is.null(Ycontrol)){
		
# computes intermediates
		
		if(is.null(Spline_B)){
			if( nBX0){
				BPHtermcontrol <-exp(BX0control %*% allparam[ibalpha0])
				modified_ratecontrol <-  expected_ratecontrol * BPHtermcontrol
				modified_cumratecontrol <- log((1 + exp( expected_logit_endcontrol))/(1 + exp(expected_logit_entercontrol))) * BPHtermcontrol  
			}
			else {
				BPHtermcontrol <-1.0
				modified_ratecontrol <-  expected_ratecontrol 
				modified_cumratecontrol <- log((1 + exp( expected_logit_endcontrol))/(1 + exp(expected_logit_entercontrol)))
			}
		}
		else {
			# parameter of the first basis is one
			brass0 <- c(1.0, allparam[ibrass0])
			S_B <- Spline_B * brass0
			Y2C <- exp(predictSpline(S_B, expected_logit_endcontrol)) 
			Y1C <- exp(predictSpline(S_B, expected_logit_entercontrol)) 
			evalderivbrasscontrol <- predictSpline(deriv(S_B), expected_logit_endcontrol)
			# E(x2) spline bases of the brass transformation at exit
			E2C <- evaluate(Spline_B, expected_logit_endcontrol)[,-1]
			# E(x1) spline bases of the brass transformation at enter
			E1C <- evaluate(Spline_B, expected_logit_entercontrol)[,-1]
			# E'(x2) derivative of the spline bases of the brass transformation at exit
			DE2C <- evaluate(deriv(Spline_B), expected_logit_endcontrol)[,-1]
			# contribution of non time dependant variables
			modified_ratecontrol <-  expected_ratecontrol * (1 + exp(-expected_logit_endcontrol))/(1+ 1/Y2C) * evalderivbrasscontrol
			# modified cumrate is computed once for each individual (from t_enter to t_end of folowup)
			modified_cumratecontrol <- log((1 + Y2C)/(1 + Y1C))


			
			
			if( nBX0){
				BPHtermcontrol <-exp(BX0control %*% allparam[ibalpha0])
				
				modified_ratecontrol <-  modified_ratecontrol * BPHtermcontrol
				# modified cumrate is computed once for each individual (from t_enter to t_end of folowup)
				modified_cumratecontrol <- modified_cumratecontrol  * BPHtermcontrol 
				
			} else {
				BPHtermcontrol <- 1
			}
			if(sum(is.na(modified_ratecontrol)) | sum(is.na(modified_cumratecontrol))){
				warning(paste0(sum(is.na(modified_ratecontrol)), 
								" NA rate control and ", 
								sum(is.na(modified_cumratecontrol)), 
								" NA cumrate control  with Brass coef", 
								paste(format(brass0), collapse = " ")))
			}
			if(min(modified_ratecontrol, na.rm=TRUE)<0 | min(modified_cumratecontrol, na.rm=TRUE)<0){
				warning(paste0(sum(modified_ratecontrol<0, na.rm=TRUE), 
								" negative rate controle and ", 
								sum(modified_cumratecontrol<0, na.rm=TRUE), 
								" negative cumrate controle with Brass coef", 
								paste(format(brass0), collapse = " ")))
			}
		}
		
		###################
		# compute dL/d brass0
		if(is.null(Spline_B)){
			dLdbrass0 <- NULL
		}
		else {
			if (!is.null(weightscontrol)) {
				dLdbrass0 <- crossprod(DE2C, Ycontrol[,3]*weightscontrol/evalderivbrasscontrol) +
						crossprod(E2C, ( Ycontrol[,3] - Y2C   * BPHtermcontrol)* weightscontrol /(1+ Y2C) ) +
						crossprod(E1C, Y1C * BPHtermcontrol * weightscontrol /(1+ Y1C) )
			} else {
				dLdbrass0 <- crossprod(DE2C, Ycontrol[,3]/evalderivbrasscontrol) +
						crossprod(E2C, ( Ycontrol[,3] - Y2C   * BPHtermcontrol)/(1+ Y2C) ) +
						crossprod(E1C, Y1C * BPHtermcontrol /(1+ Y1C) )
			}
		}
		
		if( nBX0){
			# compute dL/d balpha0
			if (!is.null(weightscontrol)) {
				dLdbalpha0 <-  crossprod(BX0control ,(Ycontrol[,3] - modified_cumratecontrol ) * weightscontrol )
			} else {
				dLdbalpha0 <-  crossprod(BX0control ,(Ycontrol[,3] - modified_cumratecontrol ) )
			}
		}
		else {
			dLdbalpha0 <- NULL
		}
		
		gr_control <- c(rep(0, length(allparam) - nbrass0 - nBX0),
				dLdbrass0,
				dLdbalpha0)
	}
	else {
		modified_ratecontrol <-  NULL
		modified_cumratecontrol <- NULL
		
		gr_control <- 0.0
	}
	
	################################################################################
	# exposed group
# Brass model
	
# computes intermediates
	if(is.null(Spline_B)){
		modified_rate <-  expected_rate 
		modified_cumrate <- log((1 + exp( expected_logit_end))/(1 + exp(expected_logit_enter)))
	}
	else {
		# parameter of the first basis is one
		brass0 <- c(1.0, allparam[ibrass0])
		S_B <- Spline_B * brass0
		Y2E <- exp(predictSpline(S_B, expected_logit_end)) 
		Y1E <- exp(predictSpline(S_B, expected_logit_enter)) 
		evalderivbrass <- predictSpline(deriv(S_B), expected_logit_end)
		# E(x2) spline bases of the brass transformation at exit
		E2E <- evaluate(Spline_B, expected_logit_end)[,-1]
		# E(x1) spline bases of the brass transformation at enter
		E1E <- evaluate(Spline_B, expected_logit_enter)[,-1]
		# E'(x2) derivative of the spline bases of the brass transformation at exit
		DE2E <- evaluate(deriv(Spline_B), expected_logit_end)[,-1]
		# contribution of non time dependant variables
		
		modified_rate <-  expected_rate * (1 + exp(-expected_logit_end))/(1+ 1/Y2E) * evalderivbrass
		modified_cumrate <- log((1 + Y2E)/(1 +  Y1E))
		
	}
	
	if( nBX0){
		BPHterm <-exp(BX0 %*% allparam[ibalpha0])
		modified_rate <-   modified_rate * BPHterm
		modified_cumrate <- modified_cumrate * BPHterm  
	}
	else {
		BPHterm <- 1.0
	}
	
	if(sum(is.na(modified_rate)) | sum(is.na(modified_cumrate))){
		warning(paste0(sum(is.na(modified_rate)), 
						" NA rate and ", 
						sum(is.na(modified_cumrate)), 
						" NA cumrate with Brass coef", 
						paste(format(brass0), collapse = " ")))
	}
	if(min(modified_rate, na.rm=TRUE)<0 | min(modified_cumrate, na.rm=TRUE)<0){
		warning(paste0(sum(modified_rate<0, na.rm=TRUE), 
						" negative rate and ", 
						sum(modified_cumrate<0, na.rm=TRUE), 
						" negative cumrate with Brass coef", 
						paste(format(brass0), collapse = " ")))
	}
	
	
	
	if(nX + nZ){
		# spline bases for each TD effect at the end of the interval
		YT <- evaluate(Spline_t, Y[,2], intercept=TRUE)
		EffectPred <- PHterm  * exp(apply(YT * Zalphabeta, 1, sum))
	} else {
		EffectPred <- PHterm  
	}
	
#  RatePred <- ifelse(Y[,3] ,
#                          EffectPred * WCEevent,
#                          0)
	RatePred <- EffectPred * WCEevent
	
	F <- ifelse(Y[,3] ,
			RatePred/(RatePred + modified_rate ), 
			0)
	FWCE <- ifelse(Y[,3] ,
			EffectPred/(RatePred + modified_rate ), 
			0)
	if(!is.null(Spline_B)){
		Ftable <- ifelse(Y[,3] ,
				modified_rate/(RatePred + modified_rate ), 
				0)
	}
# for each row i of an Id, FId[i] <- F[final_time of the id]
#  first <- unique(FirstId)
#  nline <- c(first[-1],length(FirstId)+1)-first
	
#  LastId <- FirstId+rep(nline, nline)-1
	FId <- F[LastId] 
	
	if(nX + nZ) {
		if(nX0>0) {
			Intb <- Intb * c(PHterm)
		}
		IntbF <- YT*F - Intb
	}
	else {
		IntbF <- NULL
	}
	
	
	
	#####################################################################"
# now computes the mean score and the gradients
	
#^parameters of the  correction of th elife table
	
	if(is.null(Spline_B)){
		dLdbrass0 <- NULL
	}
	else {
		if (!is.null(weights)) {
			# compute dL/d brass0
			dLdbrass0 <- crossprod(DE2E , Ftable *weights/evalderivbrass) +
					crossprod(E2E, ( Ftable - Y2E   * BPHterm)* weights /(1+ Y2E) ) +
					crossprod(E1E,     Y1E * BPHterm * weights /(1+ Y1E) )
		} else {
			# compute dL/d brass0
			dLdbrass0 <- crossprod(DE2E, Ftable / evalderivbrass) +
					crossprod(E2E, ( Ftable - Y2E   * BPHterm)/(1+ Y2E) ) +
					crossprod(E1E, Y1E * BPHterm /(1+ Y1E) )
		}
 	}
	
	if( nBX0){
		# compute dL/d balpha0
		if (!is.null(weights)) {
			dLdbalpha0 <-  crossprod(BX0 ,(Ftable - modified_cumrate ) * weights )
		} else {
			dLdbalpha0 <-  crossprod(BX0 , ( Ftable - modified_cumrate ) )
		}
	}
	else {
		dLdbalpha0 <- NULL
	}
	
	
	if (!is.null(weights)) {
		wIntbF <- IntbF * weights
		if (nX0) {
			dLdalpha0 <- crossprod(X0 , (F - PHterm * NPHterm) * weights )
		}
		else {
			dLdalpha0 <- NULL
		}
		
		if (nX){
#  traiter les Intercept_t_NPH
			IntbF * weights
			dLdbeta0 <- NULL
			for(i in 1:nX){
				if ( Intercept_t_NPH[i] ){
					dLdbeta0 <- c(dLdbeta0,  crossprod(X[,i] ,  wIntbF))
				}
				else {
					dLdbeta0 <- c(dLdbeta0, crossprod(X[,i] ,  wIntbF[,indx_without_intercept] ))
				}
			}
		}
		else {
			dLdbeta0 <- NULL
		}
		
		if (nZ) { 
			baseIntbF <- wIntbF  %*% t(tBeta)
			dLdalpha <- rep(0,getNparam(Z) )
			indZ <- getIndex(Z)
			
			for(iZ in 1:nZ){
				dLdalpha[indZ[iZ,1]:indZ[iZ,2]] <- crossprod(Z@DM[,indZ[iZ,1]:indZ[iZ,2]], baseIntbF[,iZ]  )
			}
			dLdbeta <- c(crossprod((IntbF[,-1, drop=FALSE]),Zalpha * weights))
		}
		else {
			dLdalpha <- NULL
			dLdbeta <- NULL
		}
		
# WCE effect
#       dLdeta0 <- crossprod(weights, FWCE * gradWCEevent - PHterm * gradCumWCE )
		# faster 
		dLdeta0 <- crossprod(weights * FWCE,  gradWCEevent)  - crossprod(weights * PHterm , gradCumWCE )
		
	} # end weights!=NULL
	else {
		if (nX0) {
			dLdalpha0 <- crossprod(X0 , F - PHterm* NPHterm )
		}
		else {
			dLdalpha0 <- NULL
		}
		
		if (nX){
#  traiter les Intercept_t_NPH
			dLdbeta0 <- NULL
			for(i in 1:nX){
				if ( Intercept_t_NPH[i] ){
					dLdbeta0 <- c(dLdbeta0,  crossprod(X[,i] ,  IntbF))
				}
				else {
					dLdbeta0 <- c(dLdbeta0, crossprod(X[,i] ,  IntbF[,indx_without_intercept]))
				}
			}
		}
		else {
			dLdbeta0 <- NULL
		}
		
		if (nZ) { 
			baseIntbF <- IntbF  %*% t(tBeta)
			dLdalpha <- rep(0,getNparam(Z) )
			indZ <- getIndex(Z)
			
			for(iZ in 1:nZ){
				dLdalpha[indZ[iZ,1]:indZ[iZ,2]] <- crossprod(Z@DM[,indZ[iZ,1]:indZ[iZ,2]], baseIntbF[,iZ] )
			}
			dLdbeta <- c(crossprod((IntbF[,-1, drop=FALSE]),Zalpha ))
		}
		else {
			dLdalpha <- NULL
			dLdbeta <- NULL
		}
		
		# WCE effects
# WCE effect
		if (nX0) {
			dLdeta0 <- crossprod(FWCE , gradWCEevent) - crossprod(PHterm , gradCumWCE )
		}
		else {
			dLdeta0 <- crossprod(FWCE , gradWCEevent) -  apply( gradCumWCE, 2, sum)
		}
		
	} # end weights==NULL
	
	
#print(cbind(Y, FWCE, gradWCEevent, PHterm, gradCumWCE)[1:50,])
	
	gr_exposed <- c(dLdeta0,
			dLdalpha0,          
			dLdbeta0,          
			dLdalpha,          
			dLdbeta,
			dLdbrass0,
			dLdbalpha0)
	
	
#      print("*************************************************gr_exposed")
#      print(gr_exposed)
	
	ret <- gr_control + gr_exposed
	
	
#print(cbind(allparam, ret))
#print(table(Y[,3]))
#cat("gC ")
#print(gr_control)
#cat("gE ")
#print(gr_exposed)
	
	
	
	if(debug.gr){
		attr(rep, "F") <- F 
		if(nX+nZ){
			attr(rep, "YT") <- YT
			attr(rep, "intb") <- Intb
			attr(rep, "intbF") <- IntbF
		}
		attr(rep, "RatePred") <-  RatePred
		if(debug.gr > 1000){
			cat("grad value and parameters :", "\n") 
			print(cbind(   allparam, ret))
			
		}
	}
	
	
	if ( debug.gr) {
		attr(ret, "PHterm") <- PHterm
		attr(ret, "NPHterm") <- NPHterm
		attr(ret, "modified_rate") <- modified_rate
		attr(ret, "modified_cumrate") <- modified_cumrate
		attr(ret, "gr_exposed") <- gr_exposed
		attr(ret, "modified_ratecontrol") <- modified_ratecontrol
		attr(ret, "modified_cumratecontrol") <- modified_cumratecontrol
		attr(ret, "gr_control") <- gr_control
		
		if ( debug.gr > 1000) cat("fin gr_flexrsurv_GA0B0ABE0Br0Control **", ret, "++ \n")
	}
	ret
}







