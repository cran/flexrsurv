opg_flexrsurv_fromto_GA0B0ABE0Br0Control<-function(allparam,
		Y, X0, X, Z, W,
		BX0,
		Id, FirstId, LastId=NULL,
		expected_rate,
		expected_logit_end,
		expected_logit_enter,
		weights=NULL,
		Ycontrol, BX0control, 
		weightscontrol=NULL,
		Idcontrol, FirstIdcontrol,
		expected_ratecontrol,
		expected_logit_endcontrol,
		expected_logit_entercontrol,
		step, Nstep,
		intTD=intTDft_NC, intweightsfunc=intweights_CAV_SIM,
		intTD_base=intTDft_base_NC,
		intTD_WCEbase=intTDft_WCEbase_NC,
		nT0basis,
		Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
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
	# nT0basis : number of spline basis 
	#  Spline_t0, spline object for baseline hazard, with evaluate() method
	#  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
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
	
	
#cat("************gr_flexrsurv_fromto_1WCEaddBr0Control ")
#print(format(allparam, scientific = TRUE, digits=12))
	
	
	
	################################################################################
#  excess rate
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
		YT0 <- NULL
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
		YT0 <- fevaluate(Spline_t0, Y[,2], intercept=Intercept_t0)
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
							expand=!Intercept_t_NPH,
							value=0))
		}
		else {
			Zalphabeta <- NULL
		}
	}
	
	if(nW) {
		IS_W <- ISpline_W
		eta0 <- allparam[ieta0]
		for(iW in 1:nW){
			if(Intercept_W[iW]){
				IS_W[[iW]] <- ISpline_W[[iW]] * eta0[iWbeg[iW]:iWend[iW]]
			}
			else {
				IS_W[[iW]]<- ISpline_W[[iW]] * c(0, eta0[iWbeg[iW]:iWend[iW]])
			}
			IntbW <- list()
		}
		
		if(identical(intweightsfunc , intweights_CAV_SIM,  ignore.srcref=TRUE) ){
			degree <- 2L
		}
		else if(identical(intweightsfunc , intweights_SIM_3_8,  ignore.srcref=TRUE) ){
			degree <- 3L
		}
		else if(identical(intweightsfunc , intweights_BOOLE,  ignore.srcref=TRUE) ){
			degree <- 4L
		}
		else {
			degree <- 0L
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
			
			if(is.null(Spline_t0)){
				Intb0 <- rep(0.0, dim(Y)[1])
			} else {
				Intb0 <-  intTD_base(func=rateTD_gamma0alphabetaeta0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
						Spline=Spline_t0,
						step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
						fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
						gamma0=allparam[igamma0], Zalphabeta=Zalphabeta, 
						nW = nW, W = W, eta0=allparam[ieta0], iWbeg=iWbeg, iWend=iWend,
						Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
						Spline_t = Spline_t, Intercept_t=TRUE,
						ISpline_W = IS_W, Intercept_W=Intercept_W, 
						debug=debug.gr)
			}
			
			if( identical(Spline_t0, Spline_t)){
				Intb <- Intb0
			}
			else {
				Intb <-  intTD_base(func=rateTD_gamma0alphabetaeta0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],
						Spline=Spline_t,
						step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
						fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
						gamma0=allparam[igamma0], Zalphabeta=Zalphabeta, 
						nW = nW, W = W, eta0=allparam[ieta0], iWbeg=iWbeg, iWend=iWend,
						Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
						Spline_t = Spline_t, Intercept_t=TRUE,
						ISpline_W = IS_W, Intercept_W=Intercept_W, 
						debug=debug.gr)
			}
			if(!Intercept_t0 & !is.null(Spline_t0)){
				Intb0<- Intb0[,-1]
			}
			indx_without_intercept <- 2:getNBases(Spline_t)
			
			for(iW in 1:nW){
				# in IntbW, the integrated WCE splines (parameter named Spline) are not scaled by eta0
				IntbW[[iW]]    <-  intTD_WCEbase(func=rateTD_gamma0alphabetaeta0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],  
						Spline=ISpline_W[[iW]], intercept=Intercept_W[iW], theW=W[,iW],
						step=step, Nstep=Nstep, degree=degree, intweightsfunc=intweightsfunc, 
						fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
						gamma0=allparam[igamma0], Zalphabeta=Zalphabeta, 
						nW = nW, W = W, eta0=allparam[ieta0], iWbeg=iWbeg, iWend=iWend,
						Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
						Spline_t = Spline_t, Intercept_t=TRUE,
						ISpline_W = IS_W, Intercept_W=Intercept_W, 
						debug=debug.gr)
				
			}
		} 
		else {
			NPHterm <- intTD(rateTD_gamma0eta0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
					step=step, Nstep=Nstep,
					intweightsfunc=intweightsfunc, 
					fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
					gamma0=allparam[igamma0], 
					nW = nW, W = W, eta0=allparam[ieta0], iWbeg=iWbeg, iWend=iWend,
					Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
					ISpline_W = IS_W, Intercept_W=Intercept_W)
			
			if(is.null(Spline_t0)){
				Intb0 <- rep(0.0, dim(Y)[1])
			} else {
				Intb0 <-  intTD_base(func=rateTD_gamma0eta0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
						Spline=Spline_t0,
						step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
						fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
						gamma0=allparam[igamma0], Zalphabeta=Zalphabeta, 
						nW = nW, W = W, eta0=allparam[ieta0], iWbeg=iWbeg, iWend=iWend,
						Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
						Spline_t = Spline_t, Intercept_t=TRUE,
						ISpline_W = IS_W, Intercept_W=Intercept_W, 
						debug=debug.gr)
				
				if(!Intercept_t0){
					Intb0<- Intb0[,-1]
				}
			}
			Intb <- NULL
			for(iW in 1:nW){
				# in IntbW, the integrated WCE splines are not scaled by eta0
				IntbW[[iW]]    <-  intTD_WCEbase(func=rateTD_gamma0eta0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
						Spline=ISpline_W[[iW]], intercept=Intercept_W[iW],  theW=W[,iW],
						step=step, Nstep=Nstep, degree=degree, intweightsfunc=intweightsfunc, 
						fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
						gamma0=allparam[igamma0], 
						nW = nW, W = W, eta0=allparam[ieta0], iWbeg=iWbeg, iWend=iWend,
						Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
						Spline_t = Spline_t, Intercept_t=TRUE,
						ISpline_W = IS_W, Intercept_W=Intercept_W, 
						debug=debug.gr)
			}
		}
	}
	else {
		# no VCE effect, same NPH term than ll_flexrsurv_fromto_G0A0B0AB
		if(nX + nZ) {
			NPHterm <- intTD(rateTD_gamma0alphabeta, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
					step=step, Nstep=Nstep,
					intweightsfunc=intweightsfunc, 
					gamma0=allparam[igamma0], Zalphabeta=Zalphabeta, 
					Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
					Spline_t = Spline_t, Intercept_t=TRUE)
			if(is.null(Spline_t0)){
				Intb0 <- rep(0.0, dim(Y)[1])
			} else {
				Intb0 <-  intTD_base(func=rateTD_gamma0alphabeta, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
						Spline=Spline_t0,
						step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
						gamma0=allparam[igamma0], Zalphabeta=Zalphabeta, 
						Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
						Spline_t = Spline_t, Intercept_t=TRUE,
						debug=debug.gr)
			}
			if( identical(Spline_t0, Spline_t)){
				Intb <- Intb0
			}
			else {
				Intb <-  intTD_base(func=rateTD_gamma0alphabeta, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
						Spline=Spline_t,
						step=step, Nstep=Nstep, intweightsfunc=intweightsfunc,
						gamma0=allparam[igamma0], Zalphabeta=Zalphabeta, 
						Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
						Spline_t = Spline_t, Intercept_t=TRUE)
			}
			if(!Intercept_t0 & !is.null(Spline_t0)){
				Intb0<- Intb0[,-1]
			}
			indx_without_intercept <- 2:getNBases(Spline_t)
		}
		else {
			NPHterm <- intTD(rateTD_gamma0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
					step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
					gamma0=allparam[igamma0],
					Spline_t0=Spt0g, Intercept_t0=Intercept_t0)
			if(is.null(Spline_t0)){
				Intb0 <- rep(0.0, dim(Y)[1])
			} else {
				Intb0 <-  intTD_base(func=rateTD_gamma0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3], 
						Spline=Spline_t0,
						step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
						gamma0=allparam[igamma0], 
						Spline_t0=Spt0g, Intercept_t0=Intercept_t0,
						debug=debug.gr)
				if(!Intercept_t0){
					Intb0<- Intb0[,-1]
				}
			}
			Intb <- NULL
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
	#***** 
# the contribution of the WCE to the excess rate at TFinal is in WCEContrib[LastId, ] 
	if(nW){
		# eta0 = NULL because IS_W = ISpline_W * eta0
#    WCEcontrib <- weighted_cummulative_exposure_old(Increment=W, fromT=Y[,1], finalT=Y[,4], Id=Id,
#                                                eta0=NULL, iWbeg=iWbeg, iWend=iWend, ISpline_W = IS_W)
		
		WCEcontrib <- weighted_cummulative_exposure(Increment=W, fromT=Y[,1], toT=, Y[,2], FirstId=FirstId, LastId=LastId,
				theT=Y[,4], tId=LastId,
				eta0=NULL, iWbeg=iWbeg, iWend=iWend, ISpline_W = IS_W, Intercept_W=Intercept_W)
		
	}
	else {
		WCEcontrib <- NULL
	}
	
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
								" NA cumrate control with Brass coef", 
								paste(format(brass0), collapse = " ")))
			}
			if(min(modified_ratecontrol, na.rm=TRUE)<0 | min(modified_cumratecontrol, na.rm=TRUE)<0){
				warning(paste0(sum(modified_ratecontrol<0, na.rm=TRUE), 
								" negative rate control and ", 
								sum(modified_cumratecontrol<0, na.rm=TRUE), 
								" negative cumrate control with Brass coef", 
								paste(format(brass0), collapse = " ")))
			}
		}
		
		###################
		# compute dL/d brass0
		if(is.null(Spline_B)){
			dLdbrass0 <- NULL
		}
		else {
			dLdbrass0 <- DE2C * (Ycontrol[,3]/evalderivbrasscontrol) +
					E2C * (( Ycontrol[,3] - Y2C   * BPHtermcontrol)/(1+ Y2C) ) +
					E1C * (Y1C * BPHtermcontrol /(1+ Y1C) )
		}
		
		if( nBX0){
			# compute dL/d balpha0
			dLdbalpha0 <-  BX0control * (Ycontrol[,3] - modified_cumratecontrol ) 
		}
		else {
			dLdbalpha0 <- NULL
		}
		
#    gr_control <- cbind(matrix(0, nrow=dim(Ycontrol)[1], ncol=length(allparam) - nbrass0 - nBX0),
		gr_control <- cbind(dLdbrass0,
				dLdbalpha0)
		if (!is.null(weightscontrol)) {
			opg_control<- crossprod(weightscontrol*gr_control , gr_control)
		}
		else {
			opg_control<- crossprod(gr_control)
		}
		# opg_control for other parameters are 0
		opg_control <- rbind(matrix(0, nrow=length(allparam) - nbrass0 - nBX0, ncol=length(allparam)),
				cbind(matrix(0, nrow=nbrass0 + nBX0, ncol=length(allparam) - nbrass0 - nBX0),
						opg_control))
		
		
	}
	else {
		modified_ratecontrol <-  NULL
		modified_cumratecontrol <- NULL
		
		opg_control <- 0.0
	}
#      print("*************************************************gr_control")
#  print(gr_control)
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
	
	
	# spline bases for each TD effect
	if(nX + nZ){
		# spline bases for each TD effect at the end of the interval
		YT <- evaluate(Spline_t, Y[,2], intercept=TRUE)
		if(nW){
			RatePred <- ifelse(Y[,3] ,
					PHterm * exp(YT0Gamma0 + apply(YT * Zalphabeta, 1, sum) + apply(WCEcontrib, 1, sum)),
					0)
		}
		else  {
			RatePred <- ifelse(Y[,3] ,
					PHterm * exp(YT0Gamma0 + apply(YT * Zalphabeta, 1, sum)),
					0)
		}
	} else {
		if(nW){
			RatePred <- ifelse(Y[,3] , 
					PHterm * exp(YT0Gamma0 + apply(WCEcontrib, 1, sum)),
					0)
		}
		else {
			RatePred <- ifelse(Y[,3] , 
					PHterm * exp(YT0Gamma0),
					0)
		}
	}
	
	
	F <- ifelse(Y[,3] ,
			RatePred/(RatePred + modified_rate ), 
			0)
	Ftable <- ifelse(Y[,3] ,
			modified_rate/(RatePred + modified_rate ), 
			0)
# for each row i of an Id, FId[i] <- F[final_time of the id]
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
	Intb0 <- Intb0 * c(PHterm)
	
	WF <- list()
	if(nW){
		for(i in 1:nW){
			if(nX0>0) {
				# rescale IndbW by PHterm
				IntbW[[i]] <- IntbW[[i]] * c(PHterm)
			}
			WF[[i]] <- evaluate(ISpline_W[[i]], Y[,4] - Y[,1], intercept=Intercept_W[i]) * FId
		}
	}
	else {
		WF <- NULL
	}
	
	
	#####################################################################"
# now computes the mean score and the gradients
	
#^parameters of the  correction of the life table
	
	if(is.null(Spline_B)){
		dLdbrass0 <- NULL
	}
	else {
		# compute dL/d brass0
		dLdbrass0 <- DE2E * Ftable /evalderivbrass +
				E2E * (Ftable - Y2E   * BPHterm)/(1+ Y2E)  +
				E1E * Y1E * BPHterm /(1+ Y1E) 
		
 	}
	
	if( nBX0){
		# compute dL/d balpha0
		dLdbalpha0 <-  BX0 * ( Ftable - modified_cumrate ) 
	}
	else {
		dLdbalpha0 <- NULL
	}
	
	
	
# d<dgamma0
	if(is.null(Spline_t0)){
		dLdgamma0 <- NULL
	}
	else {
		dLdgamma0 <- YT0 * F - Intb0 
	}
	if (nX0) {
		dLdalpha0 <- ( F - c(PHterm)* NPHterm ) * X0
	}
	else {
		dLdalpha0 <- NULL
	}
	
	if (nX){
#  traiter les Intercept_t_NPH
		dLdbeta0 <- NULL
		for(i in 1:nX){
			if ( Intercept_t_NPH[i] ){
				dLdbeta0 <- cbind(dLdbeta0,  X[,i] *  IntbF)
			}
			else {
				dLdbeta0 <- cbind(dLdbeta0, X[,i] *  IntbF[,indx_without_intercept])
			}
		}
	}
	else {
		dLdbeta0 <- NULL
	}
	
	if (nZ) { 
		baseIntbF <- IntbF  %*% t(tBeta)
		dLdalpha <- NULL 
		dLdbeta <- NULL 
		indZ <- getIndex(Z)
		
		for(iZ in 1:nZ){
			dLdalpha<- cbind( dLdalpha , Z@DM[,indZ[iZ,1]:indZ[iZ,2]]* baseIntbF[,iZ] )
			dLdbeta <- cbind(dLdbeta, IntbF[,-1, drop=FALSE] * Zalpha[, iZ , drop=TRUE]) 
		}
	}
	else {
		dLdalpha <- NULL
		dLdbeta <- NULL
	}
	
	# WCE effects
	
	if(nW){
		dLdeta0 <- NULL
		for(i in 1:nW){
			dLdeta0 <- cbind(dLdeta0,  W[,i] * WF[[i]] - IntbW[[i]])
		}
	}
	else{
		dLdeta0 <- NULL
	}
	
	
	
	
	gr_exposed <- cbind(dLdgamma0,          
			dLdalpha0,          
			dLdbeta0,          
			dLdalpha,          
			dLdbeta,
			dLdeta0,
			dLdbrass0,
			dLdbalpha0)
	
	if (!is.null(weights)) {
		opg_exposed<- crossprod(weights*gr_exposed , gr_exposed)
	}
	else {
		opg_exposed<- crossprod(gr_exposed)
	}
	
	opg <- opg_control + opg_exposed
	
	
	if(debug.gr){
		attr(opg, "intb0") <- Intb0
		attr(opg, "F") <- F 
		attr(opg, "YT0") <- YT0 
		if(nX+nZ){
			attr(opg, "YT") <- YT
			attr(opg, "intb") <- Intb
			attr(opg, "intbF") <- IntbF
		}
		if(nW){
			attr(opg, "intbW") <- IntbW
		}
		attr(opg, "RatePred") <-  RatePred
		
	}
	
	
	if ( debug.gr) {
		attr(opg, "PHterm") <- PHterm
		attr(opg, "NPHterm") <- NPHterm
		attr(opg, "WCEcontrib") <- WCEcontrib
		attr(opg, "modified_rate") <- modified_rate
		attr(opg, "modified_cumrate") <- modified_cumrate
		attr(opg, "opg_exposed") <- opg_exposed
		attr(opg, "modified_ratecontrol") <- modified_ratecontrol
		attr(opg, "modified_cumratecontrol") <- modified_cumratecontrol
		attr(opg, "opg_control") <- opg_control
		
		if ( debug.gr > 1000) cat("fin opg_flexrsurv_GA0B0ABE0Br0Control ** \n")
	}
#cat("************opg_flexrsurv_fromto_1WCEaddBr0Control ")
#print(cbind(allparam, opg), digits=12)
	opg
}







