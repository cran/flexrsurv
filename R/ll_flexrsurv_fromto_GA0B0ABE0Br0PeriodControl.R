ll_flexrsurv_fromto_GA0B0ABE0Br0PeriodControl<-function(allparam,
		Y, X0, X, Z, W,
		BX0,
		Id, FirstId, LastId=NULL,
		expected_rate,
		expected_logit_end,
		expected_logit_enter,
		expected_logit_end_byperiod, 
		expected_logit_enter_byperiod, 
		weights_byperiod=NULL, 
		Id_byperiod,
		weights=NULL,
		Ycontrol, BX0control, 
		weightscontrol=NULL,
		Idcontrol, FirstIdcontrol,
		expected_ratecontrol,
		expected_logit_endcontrol,
		expected_logit_entercontrol,
		expected_logit_end_byperiodcontrol, 
		expected_logit_enter_byperiodcontrol, 
		weights_byperiodcontrol, 
		Id_byperiodcontrol,
		step, Nstep,
		intTD=intTD_NC, intweightsfunc=intweights_CAV_SIM,
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
		debug=FALSE,  ...){
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
	# but for other
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
# expected_logit_end_byperiod, : expected logit of periode survival at exit of each period (used in the Brass model
# expected_logit_enter_byperiod, : expected logit of periode survival at entry of each period (used in the Brass model
# weights_byperiod,  : weight of each period (used in the Brass model weights_byperiod = weight[Id_byperiod]
# Id_byperiod,    : index in the Y object : XX_byperiod[i] corrsponds to the row Id_byperiod[i] of Y, X, Z, ...
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
# expected_logit_end_byperiodcontrol, : expected logit of periode survival at exit of each period (used in the Brass model
# expected_logit_enter_byperiodcontrol, : expected logit of periode survival at entry of each period (used in the Brass model
# weights_byperiodcontrol,  : weight of each period (used in the Brass model weights_byperiod = weight[Id_byperiod]
# Id_byperiodcontrol,    : index in the Y object : XX_byperiod[i] corrsponds to the row Id_byperiod[i] of Y, X, Z, ...
	#################################################################################################################
	# model parameters
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
	
#cat("ll_flexrsurv_fromto_GA0B0ABE0Br0Control ")
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
		# dim(WCEcondtib)== dim(W), only last line by Id is usable
#    WCEcontrib <- weighted_cummulative_exposure_old(Increment=W, fromT=Y[,1], finalT=Y[,4], fail=Y[,3], Id=Id,
#                                                eta0=NULL, iWbeg=iWbeg, iWend=iWend, ISpline_W = IS_W)
		
		WCEcontrib <- weighted_cummulative_exposure(Increment=W, fromT=Y[,1], toT=, Y[,2], FirstId=FirstId, LastId=LastId,
				theT=Y[,4], tId=LastId,
				eta0=NULL, iWbeg=iWbeg, iWend=iWend, ISpline_W = IS_W, Intercept_W=Intercept_W)
		
	} else {
		WCEcontrib <- NULL
	}
	################################################################################
	# control group
# only Brass model
	if(!is.null(Ycontrol)){
		if(is.null(Spline_B)){
			modified_ratecontrol <-  expected_ratecontrol 
			modified_cumratecontrol <- log((1 + exp( expected_logit_endcontrol))/(1 + exp(expected_logit_entercontrol)))
			modified_cumratebyPcontrol <- log((1 + exp( expected_logit_end_byperiodcontrol))/(1 + exp(expected_logit_enter_byperiodcontrol)))
		}
		else {
			# parameter of the first basis is one
			brass0 <- c(1.0, allparam[ibrass0])
			S_B <- Spline_B * brass0
			Y2C <- exp(predictSpline(S_B, expected_logit_endcontrol)) 
			evalderivbrasscontrol <- predictSpline(deriv(S_B), expected_logit_endcontrol) 
			# contribution of non time dependant variables
			modified_ratecontrol <-  expected_ratecontrol * (1 + exp(-expected_logit_endcontrol))/(1+ 1/Y2C) * evalderivbrasscontrol
			modified_cumratecontrol <- log((1 + Y2C)/(1 + exp(predictSpline(S_B, expected_logit_entercontrol))))
# by period
			Y2CbyP <- exp(predictSpline(S_B, expected_logit_end_byperiodcontrol)) 
			evalderivbrassbyPcontrol <- predictSpline(deriv(S_B), expected_logit_end_byperiodcontrol) 
			# contribution of non time dependant variables
			modified_cumratebyPcontrol <- log((1 + Y2CbyP)/(1 + exp(predictSpline(S_B, expected_logit_enter_byperiodcontrol))))
		}
		if( nBX0){
			BPHtermcontrol <-exp(BX0control %*% allparam[ibalpha0])
			modified_ratecontrol <-  modified_ratecontrol * BPHtermcontrol
			# modified cumrate is computed once for each individual (from t_enter to t_end of folowup)
			modified_cumratecontrol <- modified_cumratecontrol * BPHtermcontrol  
# by period
			modified_cumratebyPcontrol <- modified_cumratebyPcontrol * BPHtermcontrol[Id_byperiodcontrol,]  
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
		
		
		eventtermcontrol <- ifelse(Ycontrol[,3], log( modified_ratecontrol ),  0)
		if (!is.null(weightscontrol)) {
			llcontrol <- crossprod(eventtermcontrol, weightscontrol)  - crossprod(modified_cumratebyPcontrol, weights_byperiodcontrol)
		}
		else {
			llcontrol <- sum( eventtermcontrol) - sum(modified_cumratebyPcontrol )
		}
		
	}
	else {
		modified_ratecontrol <-  NULL
		modified_cumratecontrol <- NULL
		modified_cumratebyPcontrol <- NULL
		llcontrol <- 0.0
	}
	################################################################################
	# exposed group
# Brass model
	
# computes intermediates
	if(is.null(Spline_B)){
		modified_rate <-  expected_rate
		if(!is.null(expected_logit_end)){
			modified_cumrate <-log((1 + exp( expected_logit_end))/(1 + exp(expected_logit_enter)))
			modified_cumratebyP <- log((1 + exp( expected_logit_end_byperiod))/(1 + exp(expected_logit_enter_byperiod)))
		}
		else {
			modified_cumrate <-0.0
			modified_cumratebyP <-0.0
		}
	}
	else {
		brass0 <- c(1.0, allparam[ibrass0])
		S_B <- Spline_B * brass0
		Y2E <- exp(predictSpline(S_B, expected_logit_end)) 
		evalderivbrass <- predictSpline(deriv(S_B), expected_logit_end) 
		
		# contribution of non time dependant variables
		modified_rate <-  expected_rate * (1 + exp(-expected_logit_end))/(1+ 1/Y2E) * evalderivbrass
		modified_cumrate <- log((1 + Y2E)/(1 + exp(predictSpline(S_B, expected_logit_enter))))
		
		# by period
		Y2EbyP <- exp(predictSpline(S_B, expected_logit_end_byperiod)) 
		evalderivbrassbyP <- predictSpline(deriv(S_B), expected_logit_end_byperiod) 
		# contribution of non time dependant variables
		modified_cumratebyP <- log((1 + Y2EbyP)/(1 + exp(predictSpline(S_B, expected_logit_enter_byperiod))))
		
	}
	
	#proportional corrections
	if( nBX0){
		BPHterm <-exp(BX0 %*% allparam[ibalpha0])
		modified_rate <-  modified_rate  * BPHterm
		modified_cumrate <- modified_cumrate * BPHterm 
		# by period
		modified_cumratebyP <- modified_cumratebyP * BPHterm[Id_byperiod,]  
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
	
#    print("WCECONTRIB")
#    print(head(cbind(modified_rate, PHterm, WCEcontrib, NPHterm)))
#    print(head(PHterm))
#    print(summary(modified_cumrate))
	
	
	# spline bases for each TD effect
	if(nX + nZ){
		# spline bases for each TD effect at the end of the interval
		YT <- evaluate(Spline_t, Y[,2], intercept=TRUE)
		if(nW){
			eventterm <- ifelse(Y[,3] ,
					log( modified_rate + PHterm * exp(YT0Gamma0 + apply(YT * Zalphabeta, 1, sum) + apply(WCEcontrib, 1, sum)) ),
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
			llexposed <- crossprod(eventterm - PHterm * NPHterm , weights) - crossprod(modified_cumratebyP, weights_byperiod)
		} else {
			llexposed <- crossprod(eventterm - NPHterm , weights) - crossprod(modified_cumratebyP, weights_byperiod)
		}
	}
	else {
		if( nX0){
			llexposed <- sum( eventterm ) - sum( modified_cumratebyP ) - sum(PHterm * NPHterm )
		} else {
			llexposed <- sum( eventterm) - sum(modified_cumratebyP) - sum(NPHterm )
		}
	}
	
#print(cbind(PHterm, NPHterm, YT0Gamma0))
#print(c(sum(eventtermcontrol), sum(modified_cumratecontrol), sum(eventterm), sum(modified_cumrate), sum(PHterm * NPHterm), sum(PHterm ), sum(NPHterm)))
#print(c(nX0, nX, nZ, nW))      
	
#print(c( sum(eventterm), sum(modified_cumratebyP), sum(PHterm),sum(NPHterm)))
#print(summary(modified_rate))
#print(summary(modified_cumratebyP))
#print(cbind(modified_rate, modified_cumrate, modified_cumratebyP[1:length(modified_rate)]))
#print(allparam)
	
	
	ret <- llcontrol + llexposed
	
#  print(c(ret, llcontrol, llexposed))
	
	if ( debug) {
		attr(ret, "eventterm") <- eventterm
		attr(ret, "PHterm") <- PHterm
		attr(ret, "NPHterm") <- NPHterm
		attr(ret, "WCEcontrib") <- WCEcontrib
		attr(ret, "modified_rate") <- modified_rate
		attr(ret, "modified_cumrate") <- modified_cumrate
		attr(ret, "modified_cumratebyP") <- modified_cumratebyP
		attr(ret, "llexposed") <- llexposed
		attr(ret, "modified_ratecontrol") <- modified_ratecontrol
		attr(ret, "modified_cumratecontrol") <- modified_cumratecontrol
		attr(ret, "modified_cumratebyPcontrol") <- modified_cumratebyPcontrol
		attr(ret, "llcontrol") <- llcontrol
		
		if ( debug > 1000) cat("fin ll_flexrsurv_GA0B0ABE0Br0Control **", ret, "++ \n")
	}
#  cat("ll_flexrsurv_fromto_GA0B0ABE0Br0Control ")
#print(c(ret, allparam), digits=12)
	
	ret
}






