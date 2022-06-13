flexrsurv.ll.fromto.brass0.period.wceadd.fitCoptim<-function (Y, X0, X, Z, W,
		BX0, 
		Id, FirstId, LastId, 
		isEnter,
		isEnd,
		expected_rate, 
		expected_logit_end, 
		expected_logit_enter, 
		expected_logit_end_byperiod, 
		expected_logit_enter_byperiod, 
		weights_byperiod, 
		Id_byperiod,
		weights=NULL,
		Ycontrol=NULL, BX0control, 
		weightscontrol,
		Idcontrol, FirstIdcontrol, LastIdcontrol,
		isEntercontrol,
		isEndcontrol, 
		expected_ratecontrol,
		expected_logit_endcontrol,
		expected_logit_entercontrol,
		expected_logit_end_byperiodcontrol, 
		expected_logit_enter_byperiodcontrol, 
		weights_byperiodcontrol, 
		Id_byperiodcontrol,
		Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t_NPH=TRUE,
		ISpline_W =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_W=TRUE,
		Spline_B =LEBSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_B=TRUE,
		init=list(eta0=NULL, alpha0=NULL, beta0=NULL, alpha=NULL, beta=NULL,brass0=NULL, balpha0=NULL),
		optim.control=list(trace=100, REPORT=1, fnscale=-1, maxit=25), 
		Coptim.control=list(),
		method=list(int_meth="CAV_SIM", step=diff(range(Y[,ncol(Y)-1]))/100, optim_meth="BFGS",
				constOptim_meth="BFGS", lower=-Inf, upper = Inf),
		vartype = c("oim", "opg", "none"),
		varmethod = c("optim", "numDeriv.hessian", "numDeriv.jacobian"),
		numDeriv.method.args=list(eps=5e-7, d=0.001, zero.tol=sqrt(.Machine$double.eps/7e-4), r=4, v=2),
		namebrass="CorTable",
		debug=FALSE, debug.ll=FALSE, debug.gr=FALSE,
		...
) {
# flexible relative survival model using full likelihood and
# one WCE additif proportional with NLL, NPH, NPHNLL  effects
# no baseline hazard, the baseline is the WCE effet
	# the WCE link is identity (ie WCE(t) = sum b_j B_j(t-Tk) X_k
	# corection of lifetable according to generalized brass method
	# Cohort-independent generalized Brass model in an age-cohort table
	# stratified brass model eccording to fixed effects BX0 (one brass function per combination)
	# rate = brass0(expected-rate, expected_logit)*exp(BX0 balpha0) + WCE(W,t)*exp{ time-independent-effect(LL + NLL)(X0) + NPH(X) + NPHNLL(Z)}
	#
	#         input :
	#
	#
	# Y : object of class Surv but the matrix has 4 columns : beginning(1) and end(2) of interval, status(3) and end of followup(4)
	#     end of followup is assumed constant by Id
	#################################################################################################################
	# Y : object of class Surv but the matrix has 4 columns :
	# Y[,1] beginning(1) , fromT, start
	# Y[,2] end(2), toT, time
	# Y[,3] status(3) fail, status
	# Y[,4] end of followup(4) FinalT
	#     end of followup is assumed constant by Id
	#################################################################################################################
	# X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
	# X : log lineair but time dependante variable 
	# Z : object of class time d?pendent variables (spline basis expended)
	# W : Exposure variable used in Weighted Cumulative Exposure Models
	# BX0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms) for the Brass model
	# Id : varibale indicating individuals Id, lines with the same Id are considered to be from the same individual
	# FirstId : all lines in FirstId[iT]:iT in the data comes from the same individual 
	# LastId  : all lines in FirstId[iT]:LastId[iT] in the data comes from the same individual Id[iT] 
	# expected_rate : expected rate at event time T
	# expected_logit_end : expected logit of periode survival at event time T (used in the Brass model
	# expected_logit_enter : expected logit  of periode survival at enter 
	# weights : vector of weights  : LL = sum_i w_i ll_i
	# expected_logit_end_byperiod, : expected logit of periode survival at exit of each period (used in the Brass model
	# expected_logit_enter_byperiod, : expected logit of periode survival at entry of each period (used in the Brass model
	# weights_byperiod,  : weight of each period (used in the Brass model weights_byperiod = weight[Id_byperiod]
	# Id_byperiod,    : index in the Y object : XX_byperiod[i] corrsponds to the row Id_byperiod[i] of Y, X, Z, ...
	#  Spline_t0, spline object for baseline hazard, with evaluate() méthod
	#  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
	#  Spline_t, spline object for time dependant effects,  with evaluate() méthod
	#  Intercept_t=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
	#  ISpline_W, list of one integrated spline objects for the WCE effect,  with evaluate() méthod
	#  Intercept_W=TRUE, option for evaluate the ith wce effect, = TRUE all the basis, =FALSE all but first basis 
	#  Spline_B, spline objects for cohort-independent Brass function,  with evaluate() méthod
	#  Intercept_B=FALSE, option for evaluate the brass function = TRUE all the basis, =FALSE all but first basis 
	# init : list  of initial values
	# optime.control : control parameters/options for optim()
	# method : optimisation method (optim_meth) for optim(), numerical intégration method (int_meth),
	# vartype : type of variance matrix : observed inf. mat (oim inv(-H)), robust/sandwich (robust H inv(S'S) H ),
	#           outer product of the gradients (opg inv(S'S)), wher where S is the matrix of scores
	# namebrass="CorrectionTable" : used to build the names of the parameters of the brass function
	
	#
	# output : coef(with name, structure as attributes), var, conv & LL
	# coef=c(eta0, alpha0, beta0, beta, alpha, brass0, balpha0  )
	# with
	# eta0  : vector of all the coef for the WCE effects
	# alpha0 ; vector of all coefs for non time dependant variables (may contain non-loglinear terms such as spline)
	# beta0 ; matrix of all coefs for log-linear but  time dependant variables  X%*%beta0(t)
	# beta  : matrix of coefs for beta(t) nTbasis * nTDvars for NLG and NPH
	# alpha : vector of coef for alpha(z) for NLG and NPH
	# brass0 : vector of coef for the cohort independent Brass function
	# balpha0 ; vector of all coefs for non time dependant variables (may contain non-loglinear terms such as spline) of the brass model
	#           for correcting life table
	#
	
	
	vartype  <- match.arg(vartype)      # type of var
	varmethod  <- match.arg(varmethod)  # function for computing oim
	dohessian <- ( vartype=="oim") & (varmethod == "optim")
	
	
	
#   if (ncol(Y)==4) {
#     start <- Y[,1]
#     time  <- Y[,2]
#     status <- Y[,3]
#     finaltime <- Y[,4]
#   }
#   else {
#     stop("wrong dim of Y", call.=TRUE) 
#   }
	
	start <- Y[,1]
	time  <- Y[,2]
	status <- Y[,3]
	
	
	n <- nrow(Y)
	
	# structure of parameter vector
	# WCE effects
	if(!is.null(W)){
		is.WCE <- TRUE
		if(is.matrix(W)) {
			nW <- dim(W)[2]
			if(nW>1){
				stop("wrong dim of W, only one WCE effect allowed", call.=TRUE) 
			}
			W <- as.vector(W)
		} else if(is.vector(W)) {
		} else {
			stop("wrong dim of W", call.=TRUE) 
		}
		nW <- 1L
		
		nWbasis <- getNBases(ISpline_W) - 1 + Intercept_W
		neta0 <- sum(nWbasis)
		ieta0 <- 1:neta0 
		Ieta0 <- 1:neta0 
		First.eta0 <- 1
		first.eta0 <- 1
		iWend <- cumsum(nWbasis)
		iWbeg <- c(1, iWend[-nW]-1)
	} else  {
		# no WCE, constant effect WCE(t) = 1, excess model = exp(LIN+NPH+NPHNLL)
		is.WCE <- FALSE
		nWbasis <- 0L
		nW <- 0L
		neta0 <- 0L
		eta0 <- NULL
		ieta0 <- NULL
		Ieta0 <- NULL
		First.eta0 <- NULL
		first.eta0 <- NULL
		iWend <- NULL
		iWbeg <- NULL
	}
	
	
	# PHLIN and PHNLIN effects 
	if(!is.null(X0)){
		is.PH <- TRUE
		if(is.matrix(X0)) {
			nX0 <- dim(X0)[2]
		} else if(is.vector(X0)) {
			nX0 <- 1L
		} else {
			stop("wrong type of X0", call.=TRUE) 
		}
		nalpha0<-nX0
		ialpha0<-1:nX0 + neta0
		Ialpha0<-1:nX0 
		First.alpha0<-neta0+1
		first.alpha0<-1
	} else {
		is.PH <- FALSE
		nX0 <- 0L
		nalpha0 <- 0L
		alpha0 <- NULL
		ialpha0 <- NULL
		Ialpha0 <- NULL
		First.alpha0 <- NULL
		first.alpha0 <- NULL
	}
	
	# NPHLIN effects 
	if(!is.null(X)){
		is.NPHLIN <- TRUE
		nTbasis <- getNBases(Spline_t)
		nTbasis_NPH <- getNBases(Spline_t) - 1 + Intercept_t_NPH
		if(is.matrix(X)) {
			nX <- dim(X)[2]
		} else if(is.vector(X)) {
			nX <- 1L
		} else {
			stop("wrong type of X", call.=TRUE) 
		}
		nbeta0 <- sum(nTbasis_NPH)
		ibeta0 <- 1:nbeta0 + neta0 + nalpha0
		Ibeta0 <- 1:nbeta0 
		First.beta0 <- neta0 + nalpha0 + 1
		first.beta0 <- 1
	} else  {
		is.NPHLIN <- FALSE
		nTbasis <- 0L
		nTbasis_NPH <- 0L
		nX <- 0L
		nbeta0 <- 0L
		beta0 <- NULL
		ibeta0 <- NULL
		Ibeta0 <- NULL
		First.beta0 <- NULL
		first.beta0 <- NULL
	}
	
	
	# NPHNLIN effects 
	if(!is.null(Z)) {
		if( getNvar(Z)>0){
			is.NPHNLIN <- TRUE
			nTbasis <- getNBases(Spline_t)
			nZ<-getNvar(Z)
			
			nTbasis_NPHNLL <- rep(getNBases(Spline_t), nZ)
			nalpha <- getNparam(Z)
			ialpha <- 1:nalpha + neta0 + nalpha0 + nbeta0  
			Ialpha <- 1:nalpha  + nalpha0   
			First.alpha <- neta0 + nalpha0 + nbeta0 + 1  
			first.alpha <- nalpha0  + 1  
			
			# as first beta is constraints, nTbasis -1 beta per Z variable
			nbeta <- nZ * (nTbasis-1)
			ibeta <- 1:nbeta + neta0 + nalpha0 + nbeta0 + nalpha
			Ibeta <- 1:nbeta  + nbeta0 
			First.beta <- neta0 + nalpha0 + nbeta0 + nalpha + 1
			first.beta <- 1 + nbeta0 
		}
	} else  {
		is.NPHNLIN <- FALSE
#    nTbasis <- 0 already done with beta0
		nZ <- 0L
		nalpha <- 0L
		alpha <- NULL
		Ialpha <- NULL
		ialpha <- NULL
		First.alpha <- NULL
		first.alpha <- NULL
		nbeta <- 0L
		beta <- NULL
		Ibeta <- NULL
		ibeta <- NULL
		First.beta <- NULL
		first.beta <- NULL
	}
	
	
	
	
	# Brass model
	# Brass function
	if(!is.null(Spline_B)){
		nBbasis <- dim(evaluate(Spline_B, expected_logit_end, intercept=Intercept_B))[2]
		# first basis is the linear basis, with coef ==1
		nbrass0 <- nBbasis -1
		ibrass0<-1:nbrass0 + neta0 + nalpha0 + nbeta0 + nalpha + nbeta 
		Ibrass0<-1:nbrass0 
		First.brass0<-neta0 + nalpha0 + nbeta0 + nalpha + nbeta  + 1
		first.brass0<-1
	}
	else {
		# no Brass function
		nBbasis <- 0
		nbrass0 <- 0
		ibrass0<-NULL
		Ibrass0<-NULL
		First.brass0<-NULL
		first.brass0<-NULL
	}
	# Proportional corrections
	if(!is.null(BX0)){
		is.BPH <- TRUE
		if(is.matrix(BX0)) {
			nBX0 <- dim(BX0)[2]
		} else if(is.vector(BX0)) {
			nBX0 <- 1L
		} else {
			stop("wrong type of X0", call.=TRUE) 
		}
		nbalpha0<-nBX0
		ibalpha0<-1:nBX0 + neta0 + nalpha0 + nbeta0 + nalpha + nbeta + nbrass0
		Ibalpha0<-1:nBX0 
		First.balpha0<-neta0 + + nalpha0 + nbeta0 + nalpha + nbeta  + nbrass0 + 1
		first.balpha0<-1
	} else {
		is.PH <- FALSE
		nBX0 <- 0L
		nbalpha0 <- 0L
		balpha0 <- NULL
		ibalpha0 <- NULL
		Ibalpha0 <- NULL
		First.balpha0 <- NULL
		first.balpha0 <- NULL
	}
	
	
	
	# number of variables 
	nvar <- nX0 + nX + nZ +nW
	# nb of parameters
	nparam = neta0 + nalpha0 + nbeta0 + nalpha + nbeta  + nbrass0 + nbalpha0
	
	
	LastId <- getLastId(FirstId)
	if(!is.null(Ycontrol)) {
		LastIdcontrol <- getLastId(FirstIdcontrol)
	}
	# default init values for eta0 (same value at each basis)
	
	# contribution to the cumulatice WCE of each line (if eta==1)
	# it is NPHterm of the ll_function when I_Spline in unscaled and no time dependent terms
	wce2 <-  predictwce(object=integrate(ISpline_W), t=Y[,2], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
			FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE) 
	wce1 <- predictwce(object=integrate(ISpline_W), t=Y[,1], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
			FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
	cumwcebasis <- wce2 - wce1
	cumwcebasis2 <-  predictwce(object=integrate(ISpline_W), t=Y[Y[,3]==1,2], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1])[Y[,3]==1],
			FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
	
	
	initeta0 <- rep(sum(status)/sum(cumwcebasis), neta0)
	
	if (!missing(expected_rate) && !is.null(expected_rate)) {
		
		WCEevent0 <- predictwce(object=ISpline_W, t=Y[,2], Increment=W, fromT=Y[,1], tId=(1:dim(Y)[1]),
				FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
		
		
		# keep valaues only at event
		rateevent <-expected_rate[Y[,3]==1]
		WCEevent0 <- WCEevent0[Y[,3]==1]
		
		
#    initeta0 <- initeta0 - mean((rateevent/WCEevent0)[WCEevent0>0])
		initeta0 <- initeta0 - mean(rateevent)
		
	}
	
	if (missing(expected_rate) || is.null(expected_rate)) 
		expected_rate <- rep(0, n)
	if (missing(weights) || is.null(weights)) {
		weights <- NULL
	} else {
		if (any(weights <= 0)) 
			stop("Invalid weights, must be >0", call.=TRUE)
		storage.mode(weights)    <- "double"
	}
	
	# control of init values
	if (!missing(init) && !is.null(init)) {
		if (nW ){
			if (!is.null(init$eta0)) {
				if (length(init$eta0) != neta0){ 
					stop("Wrong length for initial values for eta0", call.=TRUE)
				}
				else {
					initeta0 <- init$eta0
				}
			}
			else {
				# initeta0 is the default value
			}
		}
		else {
			initeta0 <- NULL
		}
		
		if ( nX0 ) {
			if(!is.null(init$alpha0)) {
				if (length(init$alpha0) != nX0 ) {
					stop("Wrong length for initial values for alpha0", call.=TRUE)
				} else {
					initalpha0 <- init$alpha0
				}
			} else {
				initalpha0 <- rep(0, nX0)
			}
		} else {
			initalpha0 <- NULL
		}
		
		if (nX ){
			if (!is.null(init$beta0)) {
				if (length(init$beta0) != nbeta0){ 
					stop("Wrong length for initial values for beta0", call.=TRUE)
				} else {
					initbeta0 <- init$beta0
				}
			} else {
				initbeta0 <- rep(0, nbeta0)
			}
		} else {
			initbeta0 <- NULL
		}
		
		if (nZ ){
			if (!is.null(init$alpha)) {
				if (length(init$alpha) != nalpha ){ 
					stop("Wrong length for initial values for alpha", call.=TRUE)
				} else {
					initalpha <- init$alpha
				}
			} else {
				initalpha <- rep(0, nalpha)
			}        
			if (!is.null(init$beta)) {
				if (length(init$beta) != nbeta)  {
					stop("Wrong length for initial values for beta", call.=TRUE)
				} else {
					initbeta <- init$beta
				}
			} else {
				# else initbeta=0
				initbeta <- rep(0, nbeta)
			}
		}
		else {
			initalpha <- NULL
			initbeta <- NULL
		}
		
		
		if (nbrass0) {
			if (!is.null(init$brass0)) {
				if (length(init$brass0) != nbrass0){ 
					stop("Wrong length for initial values for brass0", call.=TRUE)
				}
				else {
					initbrass0 <- init$brass0
				}
			}
			else {
				initbrass0 <- rep(0, nbrass0)
			}
		} else {
			initbrass0 <- NULL
		}
		
		if (nBX0) {
			if (!is.null(init$balpha0)) {
				if (length(init$balpha0) != nbalpha0){ 
					stop("Wrong length for initial values for balpha0", call.=TRUE)
				}
				else {
					initbalpha0 <- init$balpha0
				}
			}
			else {
				initbalpha0 <- rep(0, nbalpha0)
			}
		} else {
			initbalpha0  <- NULL
		}
		
# end if (!missing(init) && !is.null(init)) 
	} else {
		# user did not provide init values
		cat("default init values used\n")
		# no init values
		if (nW) {
			# initeta0 is the default value (see above)
		} else {
			initeta0 <- NULL
		}
		
		if ( nX0) {
			initalpha0 <- rep(0, nX0)
		} else {
			initalpha0 <- NULL
		}
		if (nX) {
			initbeta0 <- rep(0, nbeta0)
		} else {
			initbeta0 <- NULL
		}
		if (nZ) {
			initalpha <- rep(0, nalpha)
			initbeta <- rep(0 , nbeta)
		} else {
			initalpha <- NULL
			initbeta <- NULL
		}
		
		if (nbrass0) {
			initbrass0 <- rep(0, nbrass0)
		} else {
			initbrass0 <- NULL
		}
		
		if (nBX0) {
			initbalpha0 <- rep(0, nbalpha0)
		} else {
			initbalpha0  <- NULL
		}
		
	}# end else if (!missing(init) && !is.null(init)) 
	
	
	storage.mode(initalpha0) <- "double"
	storage.mode(initbeta0)  <- "double"
	storage.mode(initalpha)  <- "double"
	storage.mode(initbeta)   <- "double"
	storage.mode(initeta0)   <- "double"
	storage.mode(initbrass0)   <- "double"
	storage.mode(initbalpha0)   <- "double"
	
	
	
	# numerical integration method
	# computes steps for time integtration
	if(method$int_meth == "CAV_SIM"){
		int_meth <- "NC"
		intTD <- intTDft_NC
		intTD_debug<- intTDft_NC_debug
		intTD_base<- intTDft_base_NC
		intTD_WCEbase<- intTDft_WCEbase_NC
		intweightsfunc <-intweights_CAV_SIM
		step <-method$step
		mult <- 2
	} else if(method$int_meth == "SIM_3_8"){
		int_meth <- "NC"
		intTD <- intTDft_NC
		intTD_base<- intTDft_base_NC
		intTD_WCEbase<- intTDft_WCEbase_NC
		intweightsfunc <-intweights_SIM_3_8
		step <-method$step
		mult <- 3
	} else if(method$int_meth == "BOOLE"){
		int_meth <- "NC"
		intTD <- intTDft_NC
		intTD_base<- intTDft_base_NC
		intTD_WCEbase<- intTDft_WCEbase_NC
		intweightsfunc <-intweights_BOOLE
		step <-method$step
		mult <- 4      
	} else if(method$int_meth == "Gauss-Legendre"){
		int_meth <- "GL"
		intTD <- intTDft_GL
		intTD_base<- intTDft_base_GL
		intTD_WCEbase<- intTDft_WCEbase_GL
		intweightsfunc <-NULL
		gq <- gauss.quad(method$npoints, kind="legendre")
		step <- gq$nodes
		Nstep <- gq$weights
	} else if(method$int_meth == "GLM"){
		int_meth <- "GLM"
		intTD <- intTDft_GLM
		intTD_base<- fastintTDft_base_GLM
		intTD_WCEbase<- fastintTDft_WCEbase_GLM
		intweightsfunc <- NULL
		step <-GLMStepParam(cuts=method$bands)
		Firststep <- WhichBandInf(start, step) + 1L
		Laststep  <- WhichBand(time, step)- 1L
		Nstep <- cbind(Firststep, Laststep)
	} else {
		stop("unknown method$int_meth", call.=TRUE)
	}
	
	if( int_meth == "NC"){
		STEPS<-cutTfromto(start, time, step=method$step, mult=mult)
		Nstep<-STEPS$Nstep
		step<-STEPS$step
	}
	
	
	LL <- +Inf
	
# objective, gradiant functions
	ll_function    <- ll_flexrsurv_fromto_1WCEaddBr0PeriodControl
	gr_function    <- gr_ll_flexrsurv_fromto_1WCEaddBr0PeriodControl
#  gr_function    <- NULL
	opgFunction    <- opg_flexrsurv_fromto_1WCEaddBr0PeriodControl
	computeLinearPredictor <- .computeLinearPredictor_fromto_1wceadd
	computeCumulativeHazard <- .computeCumulativeHazard_fromto_1wceadd  
	
	alpha0 <- initalpha0 
	alpha  <- initalpha 
	beta0  <- initbeta0  
	beta   <- initbeta
	eta0   <- initeta0
	brass0 <- initbrass0 
	balpha0 <- initbalpha0 
	
	initGA0B0ABE0Br0 <- c(eta0, alpha0, beta0, alpha, beta, brass0, balpha0)
	
	if(debug) cat("********************get init LL values \n")
	LL <- ll_function(allparam=initGA0B0ABE0Br0,
			Y=Y, X0=X0, X=X, Z=Z, W=W,
			BX0=BX0,
			Id=Id, FirstId=FirstId, LastId=LastId,
			isEnter=isEnter,
			isEnd=isEnd,
			expected_rate=expected_rate,
			expected_logit_enter=expected_logit_enter,
			expected_logit_end=expected_logit_end,
			expected_logit_end_byperiod=expected_logit_end_byperiod, 
			expected_logit_enter_byperiod=expected_logit_enter_byperiod, 
			weights_byperiod=weights_byperiod, 
			Id_byperiod=Id_byperiod,
			weights = weights,
			Ycontrol = Ycontrol, BX0control = BX0control, 
			weightscontrol = weightscontrol,
			Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol, LastIdcontrol = LastIdcontrol,
			isEntercontrol = isEntercontrol,
			isEndcontrol = isEndcontrol, 
			expected_ratecontrol = expected_ratecontrol,
			expected_logit_endcontrol = expected_logit_endcontrol,
			expected_logit_entercontrol = expected_logit_entercontrol,
			expected_logit_end_byperiodcontrol=expected_logit_end_byperiodcontrol, 
			expected_logit_enter_byperiodcontrol=expected_logit_enter_byperiodcontrol, 
			weights_byperiodcontrol=weights_byperiodcontrol, 
			Id_byperiodcontrol=Id_byperiodcontrol,
			step=step, Nstep=Nstep, 
			intTD=intTD, intweightsfunc=intweightsfunc,
			intTD_base=intTD_base,
			intTD_WCEbase=intTD_WCEbase,
			ialpha0=ialpha0, nX0=nX0,
			ibeta0=ibeta0, nX=nX,
			ialpha=ialpha, 
			ibeta=ibeta, 
			nTbasis=nTbasis,
			ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW,
			Spline_t = Spline_t,
			Intercept_t_NPH=Intercept_t_NPH,
			ISpline_W = ISpline_W,
			Intercept_W=Intercept_W,
			nBbasis=nBbasis,
			Spline_B=Spline_B, Intercept_B=Intercept_B,
			ibrass0=ibrass0, nbrass0=nbrass0,
			ibalpha0=ibalpha0, nBX0=nBX0,
			debug=debug.ll
	)
	if (debug) {
		print("*************** init values")
		print(initGA0B0ABE0Br0)
		cat("LL at init", LL, "\n")
	}
	
	LLold<- LL
	conv <- TRUE
	iter <- -1
	algo <- "optim"
	
	if(!is.null(gr_function)){
		GRInit <- gr_function(allparam=initGA0B0ABE0Br0,
				Y=Y, X0=X0, X=X, Z=Z, W=W,
				BX0=BX0,
				Id=Id, FirstId=FirstId,  LastId=LastId,
				isEnter=isEnter,
				isEnd=isEnd,
				expected_rate=expected_rate,
				expected_logit_enter=expected_logit_enter,
				expected_logit_end=expected_logit_end,
				expected_logit_end_byperiod=expected_logit_end_byperiod, 
				expected_logit_enter_byperiod=expected_logit_enter_byperiod, 
				weights_byperiod=weights_byperiod, 
				Id_byperiod=Id_byperiod,
				weights = weights,
				Ycontrol = Ycontrol, BX0control = BX0control, 
				weightscontrol = weightscontrol,
				Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol,  LastIdcontrol=LastIdcontrol,
				isEntercontrol = isEntercontrol,
				isEndcontrol = isEndcontrol, 
				expected_ratecontrol = expected_ratecontrol,
				expected_logit_endcontrol = expected_logit_endcontrol,
				expected_logit_entercontrol = expected_logit_entercontrol,
				expected_logit_end_byperiodcontrol=expected_logit_end_byperiodcontrol, 
				expected_logit_enter_byperiodcontrol=expected_logit_enter_byperiodcontrol, 
				weights_byperiodcontrol=weights_byperiodcontrol, 
				Id_byperiodcontrol=Id_byperiodcontrol,
				step=step, Nstep=Nstep, 
				intTD=intTD, intweightsfunc=intweightsfunc,
				intTD_base=intTD_base,
				intTD_WCEbase=intTD_WCEbase,
				ialpha0=ialpha0, nX0=nX0,
				ibeta0=ibeta0, nX=nX,
				ialpha=ialpha, 
				ibeta=ibeta, 
				nTbasis=nTbasis,
				ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW,
				ibrass0=ibrass0, nbrass0=nbrass0,
				ibalpha0=ibalpha0, nBX0=nBX0,
				Spline_t = Spline_t,
				Intercept_t_NPH=Intercept_t_NPH,
				ISpline_W = ISpline_W,
				Intercept_W=Intercept_W,
				nBbasis=nBbasis,
				Spline_B=Spline_B, Intercept_B=Intercept_B,
				debug.gr=debug.ll
		)
#                      method="L-BFGS-B",
#                      method="BFGS",
#                      lower=lowerGA0B0ABE0Br0,
#                      upper=upperGA0B0ABE0Br0,
#                      method="Nelder-Mead",
		
		if (debug) {
			cat("GR at init", "\n")
			print(cbind(initGA0B0ABE0Br0, GRInit))
			cat("GR at init", "\n")
		}
		
	}
	else {
		GRInit <- NULL
	}
	
	if(!is.null(Spline_B)) {
		# for brass model : deriv(Spline_B) > 0
		# first basis is contraints to one
		# constraits is betaclt %*% evaluate(deriv(Spline_B), PivotConstraintSplinceCLT.R2bBSplineBasis(Spline_B)[-1]) > 0
		# constraints on brass parameters : brass0[1] unconstraints (la constante)
		# brass0[-1] >= 0
		# other parameters unconstraines
		# dim(ui) = c( nbrass0-1,  nparam)
		# ui[, 1: First.brass0-1] == 0
		# ui[, First.brass0 -1+(1:nbrass0)] == evaluate(deriv(Spline_B), pivot)[-1, -1]
		# ci = evaluate(deriv(Spline_B), pivot)[-1, 1]
		# ui[, (First.brass0+nbrass0):nparam] = 0
		
		pivot <- PivotConstraintSplinceCLT.R2bBSplineBasis(Spline_B)
		matder <- evaluate(deriv(Spline_B), pivot)
		ui <- matder[-1, -1]
		ci <- - matder[-1, 1]
		if(First.brass0>1){
			ui <- cbind(matrix(0, ncol = First.brass0-1, nrow=nbrass0), ui)
		}
		if(!is.null(BX0)){
			ui <- cbind(ui, matrix(0, ncol = nBX0 , nrow=nbrass0))
		}
#		ui <- diag(nbrass0-1)
#		ci = rep(0, nbrass0-1)
#		if(First.brass0>1){
#			ui <- cbind(matrix(0, ncol = First.brass0, nrow=nbrass0-1), ui)
#		}
#		if(!is.null(BX0)){
#			ui <- cbind(ui, matrix(0, ncol = nBX0 , nrow=nbrass0-1))
#		}
		# unconstraint optimisation with optim
#   fit<-optim(initGA0B0ABE0Br0,
#              fn=ll_function,
#              gr = gr_function,
#              method=method$optim_meth,
		# constraint optimisation with constrOptim
#   fit<-constrOptim(initGA0B0ABE0Br0,
#              f=ll_function,
#              grad = gr_function,
#                    ui=ui,
#                    ci=rep(0, nbrass0-1),
#              method="L-BFGS-B",
		##              method="Nelder-Mead",
#              constraints=NULL,
#              lower=lowerGA0B0ABE0Br0,
#              upper=upperGA0B0ABE0Br0,
		
		Con <- list(outer.iterations = 100, outer.eps = 1e-05, mu = 1e-04)
		nCon <- names(Con)
		Con[(nCoptim.control <- names(Coptim.control))] <- Coptim.control
		Con <- Con[nCon]
		
#print("************************************************** fit")
#print(initGA0B0ABE0Br0)
#print(ui)
#print("************************************************** fit")
#print(optim.control)
		
		
		options(show.error.messages = FALSE)
		fit<-try(constrOptim(initGA0B0ABE0Br0,
						f=ll_function,
						grad = gr_function,
						#                         grad = NULL,
						ui=ui,
						ci=rep(0, nbrass0-1),
						method=method$constOptim_meth,
						lower=method$lower,
						upper=method$upper,
						constraints=NULL,
						outer.iterations = Con$outer.iterations, outer.eps = Con$outer.eps, mu = Con$mu,
						control = optim.control,
						hessian = dohessian,
# ll_flexrsurv_fromto_1WCEaddBr0Control args
						Y=Y, X0=X0, X=X, Z=Z, W=W,
						BX0=BX0,
						Id=Id, FirstId=FirstId,  LastId=LastId,
						isEnter=isEnter,
						isEnd=isEnd,
						expected_rate=expected_rate,
						expected_logit_enter=expected_logit_enter,
						expected_logit_end=expected_logit_end,
						expected_logit_end_byperiod=expected_logit_end_byperiod, 
						expected_logit_enter_byperiod=expected_logit_enter_byperiod, 
						weights_byperiod=weights_byperiod, 
						Id_byperiod=Id_byperiod,
						weights = weights,
						Ycontrol = Ycontrol, BX0control = BX0control, 
						weightscontrol = weightscontrol,
						Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol, LastIdcontrol=LastIdcontrol,
						isEntercontrol = isEntercontrol,
						isEndcontrol = isEndcontrol, 
						expected_ratecontrol = expected_ratecontrol,
						expected_logit_endcontrol = expected_logit_endcontrol,
						expected_logit_entercontrol = expected_logit_entercontrol,
						expected_logit_end_byperiodcontrol=expected_logit_end_byperiodcontrol, 
						expected_logit_enter_byperiodcontrol=expected_logit_enter_byperiodcontrol, 
						weights_byperiodcontrol=weights_byperiodcontrol, 
						Id_byperiodcontrol=Id_byperiodcontrol,
						step=step, Nstep=Nstep, 
						intTD=intTD, intweightsfunc=intweightsfunc,
						intTD_base=intTD_base,
						intTD_WCEbase=intTD_WCEbase,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0= ibeta0, nX=nX, 
						ialpha=ialpha, 
						ibeta= ibeta, 
						nTbasis=nTbasis,
						ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW, 
						ibrass0=ibrass0, nbrass0=nbrass0,
						ibalpha0=ibalpha0, nBX0=nBX0,
						Spline_t = Spline_t,
						Intercept_t_NPH=Intercept_t_NPH,
						ISpline_W = ISpline_W,
						Intercept_W=Intercept_W,
						nBbasis=nBbasis,
						Spline_B=Spline_B, Intercept_B=Intercept_B,
						debug=debug.ll,
						debug.gr=debug.gr),
				silent=TRUE)
		options(show.error.messages = TRUE)
		
#print("fin fit")
	}
	else {
#print("#######################################     # no Brass correction, unconstrait optimisation")
		# no Brass correction, unconstrait optimisation
		options(show.error.messages = FALSE)
		fit<-try(optim(initGA0B0ABE0Br0,
						fn=ll_function,
						gr = gr_function,
#                    gr = NULL,
						method=method$optim_meth,
						constraints=NULL,
# ll_flexrsurv_fromto_GA0B0ABE0Br0 args
						Y=Y, X0=X0, X=X, Z=Z, W=W,
						BX0=BX0,
						Id=Id, FirstId=FirstId, LastId=LastId,
						isEnter=isEnter,
						isEnd=isEnd,
						expected_rate=expected_rate,
						expected_logit_enter=expected_logit_enter,
						expected_logit_end=expected_logit_end,
						expected_logit_end_byperiod=expected_logit_end_byperiod, 
						expected_logit_enter_byperiod=expected_logit_enter_byperiod, 
						weights_byperiod=weights_byperiod, 
						Id_byperiod=Id_byperiod,
						weights = weights,
						Ycontrol = Ycontrol, BX0control = BX0control, 
						weightscontrol = weightscontrol,
						Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol, LastIdcontrol=LastIdcontrol,
						isEntercontrol = isEntercontrol,
						isEndcontrol = isEndcontrol, 
						expected_ratecontrol = expected_ratecontrol,
						expected_logit_endcontrol = expected_logit_endcontrol,
						expected_logit_entercontrol = expected_logit_entercontrol,
						expected_logit_end_byperiodcontrol=expected_logit_end_byperiodcontrol, 
						expected_logit_enter_byperiodcontrol=expected_logit_enter_byperiodcontrol, 
						weights_byperiodcontrol=weights_byperiodcontrol, 
						Id_byperiodcontrol=Id_byperiodcontrol,
						step=step, Nstep=Nstep, 
						intTD=intTD, intweightsfunc=intweightsfunc,
						intTD_base=intTD_base,
						intTD_WCEbase=intTD_WCEbase,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0= ibeta0, nX=nX, 
						ialpha=ialpha, 
						ibeta= ibeta, 
						nTbasis=nTbasis,
						ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW, 
						ibrass0=ibrass0, nbrass0=nbrass0,
						ibalpha0=ibalpha0, nBX0=nBX0,
						Spline_t = Spline_t,
						Intercept_t_NPH=Intercept_t_NPH,
						ISpline_W = ISpline_W,
						Intercept_W=Intercept_W,
						nBbasis=nBbasis,
						Spline_B=Spline_B, Intercept_B=Intercept_B,
						control = optim.control, 
						debug=debug.ll,
						debug.gr=debug.gr,
						hessian = dohessian),
				silent=TRUE)
		options(show.error.messages = TRUE)
	}
	
	if( inherits(fit, "try-error")){
		msg <- fit
		cat(msg)
		
		res <- list(coefficients = rep(NA, ),
				list_coefficients = NULL,
				linear.predictors=NULL,
				fitted.values=NULL,
				cumulative.hazard=NULL,
				var = numeric(0),
				informationMatrix= numeric(0),
				loglik = NA,
				weights = weights,
				optimfunction=algo,
				conv=msg,
				method = "flexrsurv.ll.fromto.brass.wce.fit",
				fit=fit,
				start=list(eta0=initeta0, alpha0=initalpha0, beta0=initbeta0, alpha=initalpha, beta=initbeta, 
						brass0=initbrass0, balpha0 = initbalpha0))
		res
		
	} else {
#print("LLfinal")
		LLfinal <- ll_function(allparam=fit$par,
				Y=Y, X0=X0, X=X, Z=Z, W=W,
				BX0=BX0,
				Id=Id, FirstId=FirstId, LastId=LastId,
				isEnter=isEnter,
				isEnd=isEnd,
				expected_rate=expected_rate,
				expected_logit_enter=expected_logit_enter,
				expected_logit_end=expected_logit_end,
				expected_logit_end_byperiod=expected_logit_end_byperiod, 
				expected_logit_enter_byperiod=expected_logit_enter_byperiod, 
				weights_byperiod=weights_byperiod, 
				Id_byperiod=Id_byperiod,
				weights = weights,
				Ycontrol = Ycontrol, BX0control = BX0control, 
				weightscontrol = weightscontrol,
				Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol, LastIdcontrol=LastIdcontrol,
				isEntercontrol = isEntercontrol,
				isEndcontrol = isEndcontrol, 
				expected_ratecontrol = expected_ratecontrol,
				expected_logit_endcontrol = expected_logit_endcontrol,
				expected_logit_entercontrol = expected_logit_entercontrol,
				expected_logit_end_byperiodcontrol=expected_logit_end_byperiodcontrol, 
				expected_logit_enter_byperiodcontrol=expected_logit_enter_byperiodcontrol, 
				weights_byperiodcontrol=weights_byperiodcontrol, 
				Id_byperiodcontrol=Id_byperiodcontrol,
				step=step, Nstep=Nstep, 
				intTD=intTD, intweightsfunc=intweightsfunc,
				intTD_base=intTD_base,
				intTD_WCEbase=intTD_WCEbase,
				ialpha0=ialpha0, nX0=nX0,
				ibeta0=ibeta0, nX=nX,
				ialpha=ialpha, 
				ibeta=ibeta, 
				nTbasis=nTbasis,
				ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW,
				ibrass0=ibrass0, nbrass0=nbrass0,
				ibalpha0=ibalpha0, nBX0=nBX0,
				Spline_t = Spline_t,
				Intercept_t_NPH=Intercept_t_NPH,
				ISpline_W = ISpline_W,
				Intercept_W=Intercept_W,
				nBbasis=nBbasis,
				Spline_B=Spline_B, Intercept_B=Intercept_B,
				debug=debug.ll
		)
#print("fin LLfinal")
		if(!is.null(gr_function)){
#print("GRfinal")
			GRfinal <- gr_function(allparam=fit$par,
					Y=Y, X0=X0, X=X, Z=Z, W=W,
					BX0=BX0,
					Id=Id, FirstId=FirstId, LastId=LastId,
					isEnter=isEnter,
					isEnd=isEnd,
					expected_rate=expected_rate,
					expected_logit_enter=expected_logit_enter,
					expected_logit_end=expected_logit_end,
					expected_logit_end_byperiod=expected_logit_end_byperiod, 
					expected_logit_enter_byperiod=expected_logit_enter_byperiod, 
					weights_byperiod=weights_byperiod, 
					Id_byperiod=Id_byperiod,
					weights = weights,
					Ycontrol = Ycontrol, BX0control = BX0control, 
					weightscontrol = weightscontrol,
					Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol, LastIdcontrol=LastIdcontrol,
					isEntercontrol = isEntercontrol,
					isEndcontrol = isEndcontrol, 
					expected_ratecontrol = expected_ratecontrol,
					expected_logit_endcontrol = expected_logit_endcontrol,
					expected_logit_entercontrol = expected_logit_entercontrol,
					expected_logit_end_byperiodcontrol=expected_logit_end_byperiodcontrol, 
					expected_logit_enter_byperiodcontrol=expected_logit_enter_byperiodcontrol, 
					weights_byperiodcontrol=weights_byperiodcontrol, 
					Id_byperiodcontrol=Id_byperiodcontrol,
					step=step, Nstep=Nstep, 
					intTD=intTD, intweightsfunc=intweightsfunc,
					intTD_base=intTD_base,
					intTD_WCEbase=intTD_WCEbase,
					ialpha0=ialpha0, nX0=nX0,
					ibeta0=ibeta0, nX=nX,
					ialpha=ialpha, 
					ibeta=ibeta, 
					nTbasis=nTbasis,
					ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW,
					ibrass0=ibrass0, nbrass0=nbrass0,
					ibalpha0=ibalpha0, nBX0=nBX0,
					Spline_t = Spline_t,
					Intercept_t_NPH=Intercept_t_NPH,
					ISpline_W = ISpline_W,
					Intercept_W=Intercept_W,
					nBbasis=nBbasis,
					Spline_B=Spline_B, Intercept_B=Intercept_B,
					debug.gr=debug.ll
			)
#print("finGRfinal")
		}
		else {
			GRfinal <- NULL
		}
		print("*********************************   fin")
		print(c(LLfinal, fit$value, LLfinal - fit$value))
		print(cbind(fit$par, GRfinal, GRInit))
		print("*********************************   fin")
		
		conv <- converged(fit, step="maximisation eta0 alpha0 beta0 alpha and beta")
#print("message conv")
#print(conv)
		if(nW){
			eta0<-fit$par[ieta0]
		}
		if (nX0){
			alpha0 <- fit$par[ialpha0]
		}
		if(nX){
			beta0<-fit$par[ibeta0]
		}
		if (nZ){
			alpha <- fit$par[ialpha]
			beta <- fit$par[ibeta]
		}
		if (nbrass0){
			brass0 <- fit$par[ibrass0]
		}
		if (nBX0){
			balpha0 <- fit$par[ibalpha0]
		}
		LL <- fit$value
		iter <- NA
		
		
		# prepare returned value
		if(nW){
			nameseta0 <- NULL
			for(iw in 1:nW) {
				nameseta0 <- c(nameseta0,
						paste(paste("WCE(",
										dimnames(W)[[2]][iw],
										", ",
										dimnames(Y)[[2]][ncol(Y)],
										"):",
										sep=""),
								(1-Intercept_W[iw])+1:nWbasis[iw],
								sep="")
				)
			}
			names(eta0) <-  nameseta0
		}
		
		if (nX0){
			names(alpha0) <- dimnames(X0)[[2]]
		}
		
		if(nX){
			namesbeta0 <- NULL
			for(ix in 1:nX) {
				namesbeta0 <- c(namesbeta0,
						paste(paste("NPH(",
										dimnames(X)[[2]][ix],
										", ",
										dimnames(Y)[[2]][ncol(Y)-1],
										"):",
										sep=""),
								(2-Intercept_t_NPH[ix]):nTbasis,
								sep="")
				)
			}
			names(beta0) <-  namesbeta0
		}
		
		if (nZ){
			names(alpha) <-paste("NPHNLL:", dimnames(getDesignMatrix(Z))[[2]], sep="")
			namesbeta <- NULL
			for(iz in 1:nZ) {
				namesbeta <- c(namesbeta, 
						paste(paste("NPHNLL(", getNames(Z), sep=""),
								paste(dimnames(Y)[[2]][ncol(Y)-1],
										"):",
										2:nTbasis,
										sep=""),
								sep="_"))
				names(beta) <-  namesbeta
			}
		}
		
		if (nbrass0){
			namesbrass0 <-  paste(namebrass, ":",
					(1-Intercept_B)+1:nbrass0,
					sep="")
			names(brass0) <- namesbrass0 
		}
		if (nBX0){
			names(balpha0) <- paste(namebrass, ":", dimnames(BX0)[[2]], sep="")
		}
		
		
		lcoef <- list(eta0=eta0, alpha0=alpha0, beta0=beta0 , alpha=alpha, beta=beta, brass0=brass0, balpha0=balpha0)
		coef <- c(lcoef$eta0, lcoef$alpha0, lcoef$beta0 , lcoef$alpha, lcoef$beta, lcoef$brass0, lcoef$balpha0)
		
# variance computation
#print("# variance computation")
		var <- NULL
		informationMatrix <- NULL
		cholinformationMatrix <- NULL
		
		if( vartype == "oim"){
# the variance matrix is obtained from the observed information matrix
# numericNHessian in package maxLik
			var <- NULL
			if( varmethod == "optim" ){
#print("varmethod = optim")
				informationMatrix <- - fit$hessian
#print("varmethod = optim 2")
			} else if( varmethod == "numDeriv.hessian" ){
#print("varmethod == numDeriv.hessian") 
				informationMatrix <- - hessian(ll_function,  fit$par,
						method.args=numDeriv.method.args,
						Y=Y, X0=X0, X=X, Z=Z, W=W,
						BX0=BX0,
						Id=Id, FirstId=FirstId, LastId=LastId,
						isEnter=isEnter,
						isEnd=isEnd,
						expected_rate=expected_rate,
						expected_logit_enter=expected_logit_enter,
						expected_logit_end=expected_logit_end,
						expected_logit_end_byperiod=expected_logit_end_byperiod, 
						expected_logit_enter_byperiod=expected_logit_enter_byperiod, 
						weights_byperiod=weights_byperiod, 
						Id_byperiod=Id_byperiod,
						weights = weights,
						Ycontrol = Ycontrol, BX0control = BX0control, 
						weightscontrol = weightscontrol,
						Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol, LastIdcontrol=LastIdcontrol,
						isEntercontrol = isEntercontrol,
						isEndcontrol = isEndcontrol, 
						expected_ratecontrol = expected_ratecontrol,
						expected_logit_endcontrol = expected_logit_endcontrol,
						expected_logit_entercontrol = expected_logit_entercontrol,
						expected_logit_end_byperiodcontrol=expected_logit_end_byperiodcontrol, 
						expected_logit_enter_byperiodcontrol=expected_logit_enter_byperiodcontrol, 
						weights_byperiodcontrol=weights_byperiodcontrol, 
						Id_byperiodcontrol=Id_byperiodcontrol,
						step=step, Nstep=Nstep, 
						intTD=intTD, intweightsfunc=intweightsfunc,
						intTD_base=intTD_base,
						intTD_WCEbase=intTD_WCEbase,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0=ibeta0, nX=nX,
						ialpha=ialpha, 
						ibeta=ibeta, 
						nTbasis=nTbasis,
						ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW,
						ibrass0=ibrass0, nbrass0=nbrass0,
						ibalpha0=ibalpha0, nBX0=nBX0,
						Spline_t = Spline_t,
						Intercept_t_NPH=Intercept_t_NPH,
						ISpline_W = ISpline_W,
						Intercept_W=Intercept_W,
						nBbasis=nBbasis,
						Spline_B=Spline_B, Intercept_B=Intercept_B,
						debug=debug.ll)
				
#print("varmethod == numDeriv.hessian2") 
			} else if( varmethod == "numDeriv.jacobian" ){
#print("varmethod == numDeriv.jacobian") 
				informationMatrix <- -  jacobian(gr_function,  fit$par,
						method.args=numDeriv.method.args,
						Y=Y, X0=X0, X=X, Z=Z, W=W,
						BX0=BX0,
						Id=Id, FirstId=FirstId, LastId=LastId,
						isEnter=isEnter,
						isEnd=isEnd,
						expected_rate=expected_rate,
						expected_logit_enter=expected_logit_enter,
						expected_logit_end=expected_logit_end,
						expected_logit_end_byperiod=expected_logit_end_byperiod, 
						expected_logit_enter_byperiod=expected_logit_enter_byperiod, 
						weights_byperiod=weights_byperiod, 
						Id_byperiod=Id_byperiod,
						weights = weights,
						Ycontrol = Ycontrol, BX0control = BX0control, 
						weightscontrol = weightscontrol,
						Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol, LastIdcontrol=LastIdcontrol,
						isEntercontrol = isEntercontrol,
						isEndcontrol = isEndcontrol, 
						expected_ratecontrol = expected_ratecontrol,
						expected_logit_endcontrol = expected_logit_endcontrol,
						expected_logit_entercontrol = expected_logit_entercontrol,
						expected_logit_end_byperiodcontrol=expected_logit_end_byperiodcontrol, 
						expected_logit_enter_byperiodcontrol=expected_logit_enter_byperiodcontrol, 
						weights_byperiodcontrol=weights_byperiodcontrol, 
						Id_byperiodcontrol=Id_byperiodcontrol,
						step=step, Nstep=Nstep, 
						intTD=intTD, intweightsfunc=intweightsfunc,
						intTD_base=intTD_base,
						intTD_WCEbase=intTD_WCEbase,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0=ibeta0, nX=nX,
						ialpha=ialpha, 
						ibeta=ibeta, 
						nTbasis=nTbasis,
						ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW,
						ibrass0=ibrass0, nbrass0=nbrass0,
						ibalpha0=ibalpha0, nBX0=nBX0,
						Spline_t = Spline_t,
						Intercept_t_NPH=Intercept_t_NPH,
						ISpline_W = ISpline_W,
						Intercept_W=Intercept_W,
						nBbasis=nBbasis,
						Spline_B=Spline_B, Intercept_B=Intercept_B,
						debug=debug.ll)
#print("varmethod == numDeriv.jacobian2") 
				
				informationMatrix <- (informationMatrix + t(informationMatrix))/2
#print("varmethod == numDeriv.jacobian3") 
			}
		} else if( vartype == "opg"){
			
			# the variance matrix is obtained from the outer product of the gradient
			# in the case of Type 1 censoring
			informationMatrix <- opgFunction(allparam=fit$par,
					Y=Y, X0=X0, X=X, Z=Z, W=W,
					BX0=BX0,
					Id=Id, FirstId=FirstId,  LastId=LastId,
					isEnter=isEnter,
					isEnd=isEnd,
					expected_rate=expected_rate,
					expected_logit_enter=expected_logit_enter,
					expected_logit_end=expected_logit_end,
					expected_logit_end_byperiod=expected_logit_end_byperiod, 
					expected_logit_enter_byperiod=expected_logit_enter_byperiod, 
					weights_byperiod=weights_byperiod, 
					Id_byperiod=Id_byperiod,
					weights = weights,
					Ycontrol = Ycontrol, BX0control = BX0control, 
					weightscontrol = weightscontrol,
					Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol,  LastIdcontrol=LastIdcontrol,
					isEntercontrol = isEntercontrol,
					isEndcontrol = isEndcontrol, 
					expected_ratecontrol = expected_ratecontrol,
					expected_logit_endcontrol = expected_logit_endcontrol,
					expected_logit_entercontrol = expected_logit_entercontrol,
					expected_logit_end_byperiodcontrol=expected_logit_end_byperiodcontrol, 
					expected_logit_enter_byperiodcontrol=expected_logit_enter_byperiodcontrol, 
					weights_byperiodcontrol=weights_byperiodcontrol, 
					Id_byperiodcontrol=Id_byperiodcontrol,
					step=step, Nstep=Nstep, 
					intTD=intTD, intweightsfunc=intweightsfunc,
					intTD_base=intTD_base,
					intTD_WCEbase=intTD_WCEbase,
					ialpha0=ialpha0, nX0=nX0,
					ibeta0= ibeta0, nX=nX, 
					ialpha=ialpha, 
					ibeta= ibeta, 
					nTbasis=nTbasis,
					ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW, 
					ibrass0=ibrass0, nbrass0=nbrass0,
					ibalpha0=ibalpha0, nBX0=nBX0,
					Spline_t = Spline_t,
					Intercept_t_NPH=Intercept_t_NPH,
					ISpline_W = ISpline_W,
					Intercept_W=Intercept_W,
					nBbasis=nBbasis,
					Spline_B=Spline_B, Intercept_B=Intercept_B,
					debug.gr=debug.gr
			)
		}
		
		var <- NULL
		if(vartype != "none"){
			# first try inverting from its Choleski decomposition.
#print("      # first try inverting from its Choleski decomposition.")
			cholinformationMatrix <- NULL
			options(show.error.messages = FALSE)
			cholinformationMatrix <- try(chol(informationMatrix), silent=TRUE)
#print("      cholinformationMatrix <- try(chol(informationMatrix), silent=TRUE)")
			options(show.error.messages = TRUE)
			if( !inherits(cholinformationMatrix, "try-error")){
				options(show.error.messages = FALSE)
#print("        var <- try( chol2inv(cholinformationMatrix) , silent=TRUE) 1")
				var <- try( chol2inv(cholinformationMatrix) , silent=TRUE)
#print("        var <- try( chol2inv(cholinformationMatrix) , silent=TRUE) 2")
				options(show.error.messages = TRUE)
				if( inherits(var, "try-error")){
					var <- NULL
					cat(geterrmessage())
				}
				else{
					attr(var, "type") <- vartype
					attr(var, "solve") <- "chol"
				}
			}
			if( is.null(var)){
				# second try inverting with solve()
#print("      # second try inverting with solve()")
				options(show.error.messages = FALSE)
				var <- try( solve(informationMatrix) , silent=TRUE)
#print("      # second try inverting with solve() 2")
				options(show.error.messages = TRUE)
				if( inherits(var, "try-error")){
					var <- numeric(0)
					cat(geterrmessage())
				}
				else{
					attr(var, "type") <- vartype
					attr(var, "solve") <- "solve"
				}
			}
		}
		attr(informationMatrix, "type") <- vartype
		attr(var, "type") <- vartype
		
		
		# computes the linear predictor relative to the excess hazard objfit$linear.predictor
		# Remarque, make sure that .computeLinearPredictor_GA0B0ABE0 and .computeCumulativeHazard_fromto_GA0B0ABE0
		# is able to get the wright parameters from ialpha0, ibeta0, ialpah and ibeta and iteat0,
		# even if ther is parameters for brass model in fit$par
		#print("linpred")
		linearPredictors <- computeLinearPredictor(allparam=fit$par,
				Y=Y, X0=X0, X=X, Z=Z, W=W,
				Id=Id, FirstId=FirstId, LastId=LastId,
				ialpha0=ialpha0, nX0=nX0,
				ibeta0= ibeta0, nX=nX, 
				ialpha=ialpha, 
				ibeta= ibeta, 
				ieta0=ieta0, 
				nTbasis=nTbasis,
				Spline_t = Spline_t,
				Intercept_t_NPH=Intercept_t_NPH,
				ISpline_W = ISpline_W,
				Intercept_W=Intercept_W,
				wcelink="identity",
				debug=debug)
		
		# computes the fited rate
#print("  # computes the fited rate")
		fittedValues <- exp(linearPredictors[,1])*linearPredictors[,2]
		
#  # computes the cumulative exxcess hazard ??  
#print("#  # computes the cumulative exxcess hazard ??  ")
		cumulativeHazard <- computeCumulativeHazard(allparam=fit$par,
				Y=Y, X0=X0, X=X, Z=Z, W=W,
				Id=Id, FirstId=FirstId, LastId=LastId,
				step=step, Nstep=Nstep, 
				intTD=intTD, intweightsfunc=intweightsfunc,
				ialpha0=ialpha0, nX0=nX0,
				ibeta0= ibeta0, nX=nX, 
				ialpha=ialpha, 
				ibeta= ibeta, 
				ieta0=ieta0, 
				nTbasis=nTbasis,
				Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
				ISpline_W = ISpline_W,
				Intercept_W=Intercept_W,
				debug=debug)
		
		# convergence assesment
#print("     # convergence assesment")
		if (conv != TRUE) {
			warning("Ran out of iterations and did not converge")
		}
		
		res <- list(coefficients = coef,
				list_coefficients = lcoef,
				linear.predictors=linearPredictors,
				fitted.values=fittedValues,
				cumulative.hazard=cumulativeHazard, 
				var = var,
				informationMatrix=informationMatrix,
				loglik = LL,
				weights = weights,
				optimfunction=algo,
				conv=conv,
				method = "flexrsurv.ll.fromto.brass.wce.fit",
				numerical_integration_method = method,
				fit=fit,
				start=list(eta0= initeta0, alpha0=initalpha0, beta0=initbeta0, alpha=initalpha, beta=initbeta,
						brass0=initbrass0, balpha0 = initbalpha0))
#print("fin res")    
		res
		
	}
}


