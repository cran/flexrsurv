flexrsurv.ll.fromto.brass0.period.wce.fitCoptim<-function (Y, X0, X, Z, W,
		BX0, 
		Id, FirstId,  
		expected_rate, 
		expected_logit_end, 
		expected_logit_enter, 
		expected_logit_end_byperiod, 
		expected_logit_enter_byperiod, 
		weights_byperiod, 
		Id_byperiod,
		weights=NULL,
		Ycontrol, BX0control, 
		weightscontrol,
		Idcontrol, FirstIdcontrol,
		expected_ratecontrol,
		expected_logit_endcontrol,
		expected_logit_entercontrol,
		expected_logit_end_byperiodcontrol, 
		expected_logit_enter_byperiodcontrol, 
		weights_byperiodcontrol, 
		Id_byperiodcontrol,
		Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
		Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t_NPH=TRUE,
		ISpline_W =MSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_W=TRUE,
		Spline_B =LEBSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_B=TRUE,
		bhlink=c("log", "identity"),
		init=list(gamma0= NULL, alpha0=NULL, beta0=NULL, alpha=NULL, beta=NULL,eta0=NULL, brass0=NULL, balpha0=NULL),
		fastinit=TRUE,
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
) 
{
# flexible relative survival model using full likelihood and 
# non iteratif, paramétrage identifiable
# NLL, NPH, NPHNLL and WCE effects
	# corection of lifetable according to generalized brass method
	# Cohort-independent generalized Brass model in an age-cohort table
	# stratified brass model eccording to fixed effects BX0 (one brass function per combination)
	# rate = brass0(expected-rate, expected_logit)*exp(BX0 balpha0) + exp(gamma0(t) + time-independent effect(LL + NLL)(X0) + NPH(X) + NPHNLL(Z) + WCE(W))
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
	# W : Exposure variables used in Weighted Cumulative Exposure Models
	# BX0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms) for the Brass model
	# Id : varibale indicating individuals Id, lines with the same Id are considered to be from the same individual
	# FirstId : all lines in FirstId[iT]:iT in the data comes from the same individual 
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
	#  ISpline_W, list of nW integrated spline objects for WCE effects,  with evaluate() méthod
	#  Intercept_W[i]=FALSE, option for evaluate the ith wce effect, = TRUE all the basis, =FALSE all but first basis 
	#  Spline_B, spline objects for cohort-independent Brass function,  with evaluate() méthod
	#  Intercept_B=FALSE, option for evaluate the brass function = TRUE all the basis, =FALSE all but first basis 
	# init : list  of initial values
	# fastinit : if init=NULL, when fastinit=TRUE, init=(gamma0=rep(log(sum(status)/sum(time*(status==1)), ngamma0), othercoef=0)
	#                          when fastinit=FALSE, init in 3 steps: initgamma0, initalpha0alpah, initbeta0beta
	# optime.control : control parameters/options for optim()
	# method : optimisation method (optim_meth) for optim(), numerical intégration method (int_meth),
	# vartype : type of variance matrix : observed inf. mat (oim inv(-H)), robust/sandwich (robust H inv(S'S) H ),
	#           outer product of the gradients (opg inv(S'S)), wher where S is the matrix of scores
	# namebrass="CorrectionTable" : used to build the names of the parameters of the brass function
	
	#
	# output : coef(with name, structure as attributes), var, conv & LL
	# coef=c(gamma0, alpha0, beta0, beta, alpha, eta0, brass0, balpha0  )
	# with
	# gamma0 : vector of coef for baseline hazard
	# alpha0 ; vector of all coefs for non time dependant variables (may contain non-loglinear terms such as spline)
	# beta0 ; matrix of all coefs for log-linear but  time dependant variables  X%*%beta0(t)
	# beta  : matrix of coefs for beta(t) nTbasis * nTDvars for NLG and NPH
	# alpha : vector of coef for alpha(z) for NLG and NPH
	# eta0  : vector of all the coef for the WCE effects
	# brass0 : vector of coef for the cohort independent Brass function
	# balpha0 ; vector of all coefs for non time dependant variables (may contain non-loglinear terms such as spline) of the brass model
	#           for correcting life table
	#
	
#  print(head(Y))
#  print(head(W))
#  print(head(X0))
#  print(head(cbind(Id, FirstId, expected_rate)))
#  print(summary(expected_logit_end))
#  print(summary(expected_logit_enter))
	
#  if(!is.null(W))
#    {
#      print( ISpline_W[[1]]@Matrices)
#    }
	
	bhlink  <- match.arg(bhlink)       # type baseline hazard
	
	vartype  <- match.arg(vartype)  # type of var
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
	# baseline hazard 
	if(!is.null(Spline_t0)){
		nT0basis <- dim(evaluate(Spline_t0, time, intercept=Intercept_t0))[2]
		ngamma0 <- nT0basis
	}
	else {
		# no baseline hazard
		# the baseline is in NPH, nPHNLL or WCE effects
		nT0basis <- 0
		ngamma0 <- 0
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
		ialpha0<-1:nX0 + ngamma0
		Ialpha0<-1:nX0 
		First.alpha0<-ngamma0+1
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
		ibeta0 <- 1:nbeta0 + ngamma0 + nalpha0
		Ibeta0 <- 1:nbeta0 
		First.beta0 <- ngamma0 + nalpha0 + 1
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
			ialpha <- 1:nalpha + ngamma0 + nalpha0 + nbeta0  
			Ialpha <- 1:nalpha  + nalpha0   
			First.alpha <- ngamma0 + nalpha0 + nbeta0 + 1  
			first.alpha <- nalpha0  + 1  
			
			# as first beta is constraints, nTbasis -1 beta per Z variable
			nbeta <- nZ * (nTbasis-1)
			ibeta <- 1:nbeta + ngamma0 + nalpha0 + nbeta0 + nalpha
			Ibeta <- 1:nbeta  + nbeta0 
			First.beta <- ngamma0 + nalpha0 + nbeta0 + nalpha + 1
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
	
	
	# WCE effects
	if(!is.null(W)){
		is.WCE <- TRUE
		if(is.matrix(W)) {
			nW <- dim(W)[2]
		} else if(is.vector(W)) {
			nW <- 1L
		} else {
			stop("error flexrsurv.ll.fromto.brass.wce.fit(): wrong dim of W") 
		}
		nWbasis <- rep(0, nW)
		if( nW != length(ISpline_W)){
			stop("error flexrsurv.ll.fromto.brass.wce.fit(): wrong length of ISpline_W") 
		}
		if( nW != length(Intercept_W)){
			stop("error flexrsurv.ll.fromto.brass.wce.fit(): wrong length of Intercept_W") 
		}
		for(iW in 1:nW){
			nWbasis[iW] <- getNBases(ISpline_W[[iW]]) - 1 + Intercept_W[iW]
		}
		neta0 <- sum(nWbasis)
		ieta0 <- 1:neta0 + ngamma0 + nalpha0 + nbeta0 + nalpha + nbeta
		Ieta0 <- 1:neta0 
		First.eta0 <- ngamma0 + nalpha0 + nbeta0 + nalpha + nbeta + 1
		first.eta0 <- 1
		iWend <- cumsum(nWbasis)
		iWbeg <- c(1, iWend[-nW]-1)
	} else  {
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
	
	
	
	# Brass model
	# Brass function
	if(!is.null(Spline_B)){
		nBbasis <- dim(evaluate(Spline_B, expected_logit_end, intercept=Intercept_B))[2]
		# first basis is the linear basis, with coef ==1
		nbrass0 <- nBbasis -1
		ibrass0<-1:nbrass0 + ngamma0 + nalpha0 + nbeta0 + nalpha + nbeta + neta0
		Ibrass0<-1:nbrass0 
		First.brass0<-ngamma0 + nalpha0 + nbeta0 + nalpha + nbeta + neta0 + 1
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
			stop("error flexrsurv.ll.fromto.brass.periode.wce.fit(): wrong type of X0") 
		}
		nbalpha0<-nBX0
		ibalpha0<-1:nBX0 + ngamma0 + nalpha0 + nbeta0 + nalpha + nbeta + neta0 + nbrass0
		Ibalpha0<-1:nBX0 
		First.balpha0<-ngamma0 + + nalpha0 + nbeta0 + nalpha + nbeta + neta0 + nbrass0 + 1
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
	nparam = ngamma0 + nalpha0 + nbeta0 + nalpha + nbeta + neta0 + nbrass0 + nbalpha0
	
	
	if (missing(expected_rate) || is.null(expected_rate)) 
		expected_rate <- rep(0, n)
	if (missing(weights) || is.null(weights)) {
		weights <- NULL
	} else {
		if (any(weights <= 0)) 
			stop("Invalid weights, must be >0", call.=TRUE)
		storage.mode(weights)    <- "double"
	}
	
# compute init values for gamma0
	lambda <- sum(status)/sum(time*(status==1))
	gamma0init <- log(lambda)
	# control of init values
	
	if (!missing(init) && !is.null(init)) {
		if(ngamma0>0){
			if ( !is.null(init$gamma0)) {
				if (length(init$gamma0) != nT0basis){ 
					stop("Wrong length for initial values for gamma0", call.=TRUE)
				} else {
					initgamma0 <- init$gamma0
					do.init.gamma0 <- FALSE
				}
			} else {
				initgamma0 <- initcoef(Spline_t0, ncol=1L, init=gamma0init, intercept=Intercept_t0)
				do.init.gamma0 <- !fastinit
			}
		}
		else {
			initgamma0 <- NULL
			do.init.gamma0 <- FALSE
		}
		if ( nX0 ) {
			if(!is.null(init$alpha0)) {
				if (length(init$alpha0) != nX0 ) {
					stop("Wrong length for initial values for alpha0", call.=TRUE)
				} else {
					initalpha0 <- init$alpha0
					do.init.alpha0 <- FALSE
				}
			} else {
				initalpha0 <- rep(0, nX0)
				do.init.alpha0 <- !fastinit
			}
		} else {
			initalpha0 <- NULL
			do.init.alpha0 <- FALSE
		}
		
		if (nX ){
			if (!is.null(init$beta0)) {
				if (length(init$beta0) != nbeta0){ 
					stop("Wrong length for initial values for beta0", call.=TRUE)
				} else {
					initbeta0 <- init$beta0
					do.init.beta0 <- FALSE
				}
			} else {
				initbeta0 <- rep(0, nbeta0)
				do.init.beta0 <-  !fastinit
			}
		} else {
			initbeta0 <- NULL
			do.init.beta0 <- FALSE
		}
		
		if (nZ ){
			if (!is.null(init$alpha)) {
				if (length(init$alpha) != nalpha ){ 
					stop("Wrong length for initial values for alpha", call.=TRUE)
				} else {
					initalpha <- init$alpha
					do.init.alpha <- FALSE
				}
			} else {
				initalpha <- rep(0, nalpha)
				do.init.alpha <-  !fastinit
			}        
			if (!is.null(init$beta)) {
				if (length(init$beta) != nbeta)  {
					stop("Wrong length for initial values for beta", call.=TRUE)
				} else {
					initbeta <- init$beta
					do.init.beta <- FALSE
				}
			} else {
				# if initalpha given (do.init.alpha=0), initbeta should be such that beta(t)=1
				# else initbeta=0
#        initbeta <- (1-do.init.alpha) * initcoefC(Spline_t, ncol=nZ, intercept=FALSE)
				initbeta <- rep(0, nbeta)
				do.init.beta <-  !fastinit
			}
		}
		else {
			initalpha <- NULL
			initbeta <- NULL
			do.init.alpha <- FALSE
			do.init.beta <- FALSE
		}
		
		if (nW ){
			if (!is.null(init$eta0)) {
				if (length(init$eta0) != neta0){ 
					stop("flexrsurv.ll.fromto.brass.period.wce.fit(): Wrong length for initial values for eta0")
				}
				else {
					initeta0 <- init$eta0
					do.init.eta0 <- FALSE
				}
			}
			else {
				initeta0 <- rep(0, neta0)
				do.init.eta0 <-  !fastinit
			}
		}
		else {
			initeta0 <- NULL
			do.init.eta0 <- FALSE
		}
		
		if (nbrass0) {
			if (!is.null(init$brass0)) {
				if (length(init$brass0) != nbrass0){ 
					stop("flexrsurv.ll.fromto.brass.period.wce.fit(): Wrong length for initial values for brass0")
				}
				else {
					initbrass0 <- init$brass0
					do.init.brass0 <- FALSE
				}
			}
			else {
				# last coef is the slope of the right asymptote
				initbrass0 <- c(rep(0, nbrass0-1), 1.0)
				do.init.brass0 <- !fastinit
			}
		} else {
			initbrass0 <- NULL
			do.init.brass0 <- FALSE
		}
		
		if (nBX0) {
			if (!is.null(init$balpha0)) {
				if (length(init$balpha0) != nbalpha0){ 
					stop("flexrsurv.ll.fromto.brass.period.wce.fit(): Wrong length for initial values for balpha0")
				}
				else {
					initbalpha0 <- init$balpha0
					do.init.balpha0 <- FALSE
				}
			}
			else {
				initbalpha0 <- rep(0, nbalpha0)
				do.init.balpha0 <- !fastinit
			}
		} else {
			initbalpha0  <- NULL
			do.init.alpha0 <- FALSE
		}
		
		
		
# end if (!missing(init) && !is.null(init)) {
	} else {
		cat("no init")
		# no init values
		do.init <- !fastinit
		if(ngamma0>0){
			initgamma0 <- gamma0init * initcoef(Spline_t0, ncol=1L, intercept=Intercept_t0)
			optim.control.gamma0 <- optim.control
			optim.control.gamma0$parscale <- NULL
			optim.control.gamma0$ndeps <- NULL
			
			do.init.gamma0 <- !fastinit
		}else {
			initgamma0 <- NULL
			do.init.gamma0 <- FALSE
		}
		if ( nX0) {
			initalpha0 <- rep(0, nX0)
			do.init.alpha0 <- !fastinit
		} else {
			initalpha0 <- NULL
			do.init.alpha0 <- FALSE
		}
		if (nX) {
			initbeta0 <- rep(0, nbeta0)
			do.init.beta0 <- !fastinit
		} else {
			initbeta0 <- NULL
			do.init.beta0 <- FALSE
		}
		if (nZ) {
			initalpha <- rep(0, nalpha)
			initbeta <- rep(0 , nbeta)
			do.init.alpha <- !fastinit
			do.init.beta <- !fastinit
		} else {
			initalpha <- NULL
			initbeta <- NULL
			do.init.alpha <- FALSE
			do.init.beta <- FALSE
		}
		
		if (nW) {
			initeta0 <- rep(0, neta0)
			do.init.eta0 <- !fastinit
		} else {
			initeta0 <- NULL
			do.init.eta0 <- FALSE
		}
		
		if (nbrass0) {
			initbrass0 <- rep(0, nbrass0)
			do.init.brass0 <- !fastinit
		} else {
			initbrass0 <- NULL
			do.init.brass0 <- FALSE
		}
		
		if (nBX0) {
			initbalpha0 <- rep(0, nbalpha0)
			do.init.balpha0 <- !fastinit
		} else {
			initbalpha0  <- NULL
			do.init.alpha0 <- FALSE
		}
		
	}# end else if (!missing(init) && !is.null(init)) {
	
	optim.control.gamma0 <- optim.control
	optim.control.gamma0$parscale <- NULL
	optim.control.gamma0$ndeps <- NULL
	optim.control.alpha <- optim.control
	optim.control.alpha$parscale <- NULL
	optim.control.alpha$ndeps <- NULL
	optim.control.beta <- optim.control
	optim.control.beta$parscale <- NULL
	optim.control.beta$ndeps <- NULL
	
	do.init <- do.init.gamma0 | do.init.alpha0 | do.init.beta0 | do.init.alpha | do.init.beta  | do.init.eta0 
	
	storage.mode(initgamma0) <- "double"
	storage.mode(initalpha0) <- "double"
	storage.mode(initbeta0)  <- "double"
	storage.mode(initalpha)  <- "double"
	storage.mode(initbeta)   <- "double"
	storage.mode(initeta0)   <- "double"
	storage.mode(initbrass0)   <- "double"
	storage.mode(initbalpha0)   <- "double"
	
	
	
	# numerical integration method
	# computes steps for time integtration
	step <- Firststep <- Laststep  <-  Nstep <- mult <- NULL
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
		stop("flexrsurv.ll.fromto.brass.period.wce.fit(): unknown method$int_meth")
	}
	
	if( int_meth == "NC"){
		STEPS<-cutTfromto(start, time, step=method$step, mult=mult)
		Nstep<-STEPS$Nstep
		step<-STEPS$step
	}
	
	
	LL <- +Inf
	
# objective, gradiant functions
	if( bhlink == "log"){
		ll_function    <- ll_flexrsurv_fromto_GA0B0ABE0Br0PeriodControl
		ll_gamma0      <- ll_flexrsurv_fromto_gamma0
		ll_alpha0alpha <- ll_flexrsurv_fromto_alpha0alpha
		ll_beta0beta   <- ll_flexrsurv_fromto_beta0beta
		gr_function    <- gr_ll_flexrsurv_fromto_GA0B0ABE0Br0PeriodControl
#    gr_function    <- NULL
#	if( vartype == "opg"){
#		stop("flexrsurv.ll.fromto.brass.period.wce.fitCoptim(): vartype == 'opg' not implemented for WCE effect")
#	}
		opgFunction    <- opg_flexrsurv_fromto_GA0B0ABE0Br0PeriodControl
#    opgFunction    <- NULL
	}
	else {
		stop("flexrsurv.ll.fromto.brass0.period.wce.fitCoptim(): bhlink = 'identity' not implemented for WCE effect")
#    ll_function    <- ll_flexrsurv_fromto_GA0B0ABBr0_bh
#    ll_gamma0      <- ll_flexrsurv_fromto_gamma0_bh
		##    ll_alpha0alpha <- ll_flexrsurv_fromto_alpha0alpha_bh
		##    ll_beta0beta   <- ll_flexrsurv_fromto_beta0beta_bh
#    ll_alpha0alpha <- ll_flexrsurv_fromto_alpha0alpha
#    ll_beta0beta   <- ll_flexrsurv_fromto_beta0beta
#    gr_function    <- gr_ll_flexrsurv_fromto_GA0B0ABE0Br0_bh
#    opgFunction    <- opg_flexrsurv_fromto_GA0B0ABE0Br0_bh
	}
	
	#initial value
	if(do.init){
		init_hessian <- FALSE
		# fit model to find gamma0, all aother parameters kept = 0
		if(do.init.gamma0){
			
			
			fit1 <- optim(par=initgamma0, fn=ll_gamma0, gr = NULL, 
					method=method$optim_meth,
					constraints=NULL, 
# ll_flexrsurv_gamma0 args
					alpha0=initalpha0, beta0=initbeta0, alpha=initalpha, beta=initbeta ,
					Y=Y, X0=X0, X=X, Z=Z, 
					expected_rate=expected_rate,
					weights = weights,
					step=step, Nstep=Nstep, 
					intTD=intTD, intweightsfunc=intweightsfunc,
					nT0basis=nT0basis,
					Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
					nX0=nX0,
					nX=nX, 
					nTbasis=nTbasis,
					Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
					control = optim.control.gamma0, 
					hessian = init_hessian, debug=debug.ll
			)
			convgamma0 <- converged(fit1, step="init of gamma0")
			gamma0 <- fit1$par
			LL <- fit1$value
			# end of  if(do.init.gamma0)
		} else {
			gamma0 <- initgamma0
			LL <- +Inf
		}
		
		if(nalpha0+nalpha>0 & (do.init.alpha0 | do.init.alpha)) {
			# fit model to find alpha0 and alpha | previous gamm0 and beta0 and beta
			if (nZ) {
#          initbeta <- initcoefC(Spline_t, ncol=nZ , intercept= Intercept_t)
				initbeta <- rep(0,nbeta)
			}
			initalpha0alpha <- c(initalpha0, initalpha)
			fit2 <- optim(initalpha0alpha, fn=ll_alpha0alpha, gr = NULL, 
					method=method$optim_meth,
					constraints=NULL, 
# ll_flexrsurv_alpha0alpha args
					gamma0=gamma0, beta0=initbeta0, beta=initbeta,
					Y=Y, X0=X0, X=X, Z=Z, 
					expected_rate=expected_rate,
					weights = weights,
					step=step, Nstep=Nstep, 
					intTD=intTD, intweightsfunc=intweightsfunc,
					nT0basis=nT0basis,
					Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
					ialpha0=Ialpha0, nX0=nX0,
					nX=nX, 
					ialpha= Ialpha,
					nTbasis=nTbasis,
					Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
					control = optim.control.alpha,
					debug=debug.ll,
					hessian = init_hessian)
			
			convgalpha0alpha <- converged(fit2, step="init of alpha0 and alpha")
			LL <- fit2$value
			
			if(nX0){
				alpha0 <- fit2$par[Ialpha0]
			}
			if (nZ){
				alpha <- fit2$par[Ialpha]
			}
			# if initvalues for alpha found, update init value for beta
			if(do.init.alpha){
				do.init.beta <- TRUE
			}
		} else {
			alpha0 <- initalpha0
			alpha <- initalpha
			LL <- +Inf
		}
		
		
		if(nbeta0+nbeta>0 & (do.init.beta0 | do.init.beta)){
			# fit model with beta0 beta | gamma0 alpha0 alpha
			initbeta0beta <- c(initbeta0, initbeta)
			fit3<-optim(initbeta0beta, fn=ll_beta0beta, gr = NULL, 
					method=method$optim_meth,
					constraints=NULL, 
# ll_flexrsurv_beta0beta args
					alpha0=alpha0, alpha=alpha, gamma0=gamma0,
					Y=Y, X0=X0, X=X, Z=Z, 
					expected_rate=expected_rate,
					weights = weights,
					step=step, Nstep=Nstep, 
					intTD=intTD, intweightsfunc=intweightsfunc,
					nT0basis=nT0basis,
					Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
					nX0=nX0,
					ibeta0=Ibeta0, 
					nX=nX, 
					ibeta= Ibeta, 
					nTbasis=nTbasis,
					Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
					control = optim.control.beta, 
					debug=debug.ll,
					hessian = init_hessian)
			convbeta0beta<-converged(fit3, step="init of beta0 and bata")
			if(nX){
				beta0<-fit3$par[Ibeta0]
			}
			if (nZ){
				beta <- fit3$par[Ibeta]
			}     
			LL <- fit3$value
		} else {
			beta0 <- initbeta0
			beta <- initbeta
			LL <- +Inf
		}
		
		if(neta0>0){
			# no init procedure for WCE effect
			eta0 <- initeta0
		}
		
		if(nbrass0>0){
			# no init procedure for Brass corrections
			brass0 <- initbrass0
		}
		
		if (nBX0){
			balpha0 <- initbalpha0 
		}    
		
# end   if(do.init)
	} else {
		gamma0 <- initgamma0
		alpha0 <- initalpha0 
		alpha  <- initalpha 
		beta0  <- initbeta0  
		beta   <- initbeta
		eta0   <- initeta0
		brass0 <- initbrass0 
		balpha0 <- initbalpha0 
	}
	
	initGA0B0ABE0Br0 <- c(gamma0, alpha0, beta0, alpha, beta, eta0, brass0, balpha0)
	
	
	
	if(debug) cat("********************get init LL values \n")
	LL <- ll_function(allparam=initGA0B0ABE0Br0,
			Y=Y, X0=X0, X=X, Z=Z, W=W,
			BX0=BX0,
			Id=Id, FirstId=FirstId,
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
			Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol,
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
			nT0basis=nT0basis,
			Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
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
	
	if (debug) {
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
				Id=Id, FirstId=FirstId,
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
				Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol,
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
				nT0basis=nT0basis,
				Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
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
		
	}
	else {
		GRInit <- NULL
	}
	
	if(!is.null(Spline_B)){
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
		
		Con <- list(outer.iterations = 100, outer.eps = 1e-05)
		nCon <- names(Con)
		Con[(nCoptim.control <- names(Coptim.control))] <- Coptim.control
		Con <- Con[nCon]
		
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
						outer.iterations = Con$outer.iterations, outer.eps = Con$outer.eps,
						control = optim.control,
						hessian = dohessian,
						# ll_flexrsurv_fromto_GA0B0ABE0Br0 args
						Y=Y, X0=X0, X=X, Z=Z, W=W,
						BX0=BX0,
						Id=Id, FirstId=FirstId, 
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
						Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol,
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
						nT0basis=nT0basis,
						Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
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
		
		
	}
	else {
		# no Brass correction, unconstrait optimisation 
		options(show.error.messages = FALSE)
		fit<-try(optim(initGA0B0ABE0Br0,
						fn=ll_function,
						gr = gr_function,
#                    gr = NULL,
						method=method$optim_meth,
						lower=method$lower,
						upper=method$upper,
						constraints=NULL,
# ll_flexrsurv_fromto_GA0B0ABE0Br0 args
						Y=Y, X0=X0, X=X, Z=Z, W=W,
						BX0=BX0,
						Id=Id, FirstId=FirstId, 
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
						Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol,
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
						nT0basis=nT0basis,
						Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
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
#print("fin optim")
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
				start=list(gamma0= initgamma0, alpha0=initalpha0, beta0=initbeta0, alpha=initalpha, beta=initbeta, eta0=initeta0,
						brass0=initbrass0, balpha0 = initbalpha0))
		res
		
	} else {
#print("LLfinal")
		LLfinal <- ll_function(allparam=fit$par,
				Y=Y, X0=X0, X=X, Z=Z, W=W,
				BX0=BX0,
				Id=Id, FirstId=FirstId,
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
				Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol,
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
				nT0basis=nT0basis,
				Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
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
#print(LLfinal)
		
		if(!is.null(gr_function)){
#print("GRfinal")
			GRfinal <- gr_function(allparam=fit$par,
					Y=Y, X0=X0, X=X, Z=Z, W=W,
					BX0=BX0,
					Id=Id, FirstId=FirstId,
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
					Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol,
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
					nT0basis=nT0basis,
					Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
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
#print("fin GRfinal")
#print(GRfinal)
		}
		else {
			GRfinal <- NULL
		}
		conv <- converged(fit, step="maximisation gamma0 alpha0 beta0 alpha and beta")
#print("conv msg")
#print(conv)
		if(nT0basis){
			gamma0<-fit$par[1:nT0basis]
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
		if(nW){
			eta0<-fit$par[ieta0]
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
		if(nT0basis){
			names(gamma0) <- names(fit$par[1:nT0basis])
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
		if (nbrass0){
			namesbrass0 <-  paste(namebrass, ":",
					(1-Intercept_B)+1:nbrass0,
					sep="")
			names(brass0) <- namesbrass0 
		}
		if (nBX0){
			names(balpha0) <- paste(namebrass, ":", dimnames(BX0)[[2]], sep="")
		}
		
		
		lcoef <- list(gamma0=gamma0, alpha0=alpha0, beta0=beta0 , alpha=alpha, beta=beta, eta0=eta0, brass0=brass0, balpha0=balpha0)
		coef <- c(lcoef$gamma0, lcoef$alpha0, lcoef$beta0 , lcoef$alpha, lcoef$beta, lcoef$eta0, lcoef$brass0, lcoef$balpha0)
		
		
# variance computation
#print("########## variance computation")
#print(fit$par)
		var <- NULL
		informationMatrix <- NULL
		cholinformationMatrix <- NULL
		if( vartype == "oim"){
# the variance matrix is obtained from the observed information matrix
# numericNHessian in package maxLik
			var <- NULL  
			if( varmethod == "optim" ){ 
				informationMatrix <- - fit$hessian 
			} else if( varmethod == "numDeriv.hessian" ){ 
#print("varmethod == numDeriv.hessian")
				informationMatrix <- - hessian(ll_function, fit$par,
						method.args=numDeriv.method.args,
						Y=Y, X0=X0, X=X, Z=Z, W=W,
						BX0=BX0,
						Id=Id, FirstId=FirstId,
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
						Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol,
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
						nT0basis=nT0basis,
						Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
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
						Spline_B=Spline_B, Intercept_B=Intercept_B)
			} else if( varmethod == "numDeriv.jacobian" ){
#print("varmethod == numDeriv.jacobian" )
				informationMatrix <- -  jacobian(gr_function,  fit$par,
						method.args=numDeriv.method.args,
						Y=Y, X0=X0, X=X, Z=Z, W=W,
						BX0=BX0,
						Id=Id, FirstId=FirstId,
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
						Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol,
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
						nT0basis=nT0basis,
						Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
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
						Spline_B=Spline_B, Intercept_B=Intercept_B)
				
#print("informationMatrix")
#print(informationMatrix)
				informationMatrix <- (informationMatrix + t(informationMatrix))/2.0
			}
#print("informationMatrix")
#print(informationMatrix)
			options(show.error.messages = FALSE)
			cholinformationMatrix <- try(chol(informationMatrix), silent=TRUE)
			options(show.error.messages = TRUE)
			if( inherits(cholinformationMatrix, "try-error")){
				var <- numeric(0)
				cat(geterrmessage())
			} else {
				options(show.error.messages = FALSE)
				var <- try( chol2inv(cholinformationMatrix) , silent=TRUE)
				options(show.error.messages = TRUE)
				if( inherits(var, "try-error")){
					var <- numeric(0)
					cat(geterrmessage())
				}
				else{
					attr(var, "type") <- vartype
				}
			}
		} else if( vartype == "opg"){
			
			# the variance matrix is obtained from the outer product of the gradient
			# in the case of Type 1 censoring
			informationMatrix <- opgFunction(allparam=fit$par,
					Y=Y, X0=X0, X=X, Z=Z, W=W,
					BX0=BX0,
					Id=Id, FirstId=FirstId, 
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
					Idcontrol = Idcontrol, FirstIdcontrol = FirstIdcontrol,  
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
					nT0basis=nT0basis,
					Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
					ialpha0=ialpha0, nX0=nX0,
					ibeta0= ibeta0, nX=nX, 
					ialpha=ialpha, 
					ibeta= ibeta, 
					nTbasis=nTbasis,
					ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW, 
					ibrass0=ibrass0, nbrass0=nbrass0,
					ibalpha0=ibalpha0, nBX0=nBX0,
					Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
					ISpline_W = ISpline_W,
					Intercept_W=Intercept_W,
					nBbasis=nBbasis,
					Spline_B=Spline_B, Intercept_B=Intercept_B,
					debug.gr=debug.gr
			)
			var <- NULL
			cholinformationMatrix <- NULL
			options(show.error.messages = FALSE)
			cholinformationMatrix <- try(chol(informationMatrix), silent=TRUE)
			options(show.error.messages = TRUE)
			if( inherits(cholinformationMatrix, "try-error")){
				cat(geterrmessage())
				var <- numeric(0)
			} else {
				options(show.error.messages = FALSE)
				var <- try( chol2inv(cholinformationMatrix) , silent=TRUE)
				options(show.error.messages = TRUE)
				if( inherits(var, "try-error")){
					cat(geterrmessage())
					var <- numeric(0)
				}
				else {
					attr(var, "type") <- vartype
				}
			}
		} else {
			informationMatrix <- numeric(0)
			var <- numeric(0)
		}
		attr(informationMatrix, "type") <- vartype
		attr(var, "type") <- vartype
		
		
		# computes the linear predictor relative to the excess hazard objfit$linear.predictor
		# Remarque, make sure that .computeLinearPredictor_GA0B0ABE0 and .computeCumulativeHazard_fromto_GA0B0ABE0
		# is able to get the wright parameters from ialpha0, ibeta0, ialpah and ibeta and iteat0,
		# even if ther is parameters for brass model in fit$par 
		linearPredictors <- .computeLinearPredictor_GA0B0ABE0(allparam=fit$par,
				Y=Y, X0=X0, X=X, Z=Z, W=W,
				Id=Id, FirstId=FirstId, LastId=getLastId(FirstId),
				nT0basis=nT0basis,
				Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
				ialpha0=ialpha0, nX0=nX0,
				ibeta0= ibeta0, nX=nX, 
				ialpha=ialpha, 
				ibeta= ibeta, 
				nTbasis=nTbasis,
				ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW, 
				Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
				ISpline_W = ISpline_W,
				Intercept_W=Intercept_W,
				debug=debug)
		
		# computes the fited rate
		fittedValues <- exp(linearPredictors)
		
#  # computes the cumulative exxcess hazard ??  
		cumulativeHazard <- .computeCumulativeHazard_fromto_GA0B0ABE0(allparam=fit$par,
				Y=Y, X0=X0, X=X, Z=Z, W=W,
				Id=Id, FirstId=FirstId, LastId=getLastId(FirstId),
				step=step, Nstep=Nstep, 
				intTD=intTD, intweightsfunc=intweightsfunc,
				nT0basis=nT0basis,
				Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
				ialpha0=ialpha0, nX0=nX0,
				ibeta0= ibeta0, nX=nX, 
				ialpha=ialpha, 
				ibeta= ibeta, 
				nTbasis=nTbasis,
				ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW, 
				Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
				ISpline_W = ISpline_W,
				Intercept_W=Intercept_W,
				debug=debug)
		
		# convergence assesment
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
				bhlink=bhlink,
				optimfunction=algo,
				conv=conv,
				method = "flexrsurv.ll.fromto.brass.wce.fit",
				numerical_integration_method = method,
				fit=fit,
				start=list(gamma0= initgamma0, alpha0=initalpha0, beta0=initbeta0, alpha=initalpha, beta=initbeta, eta0=initeta0,
						brass0=initbrass0, balpha0 = initbalpha0))
		
		res
		
	}
}


