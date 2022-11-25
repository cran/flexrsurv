# fit of a flexible relative survival model with a model of correction of life table 

flexrsurvclt.ll <- function(formula=formula(data),
		formula.table=NULL, 
		data=parent.frame(),
		Id,
		# knots.Bh and degree.Bh allow to define the Baseline Hazard
		baselinehazard=TRUE,
		firstWCEIadditive=FALSE,
		knots.Bh,   
		degree.Bh=3,
		Spline=c("b-spline", "tp-spline", "tpi-spline"), # tp-spline for truncated power basis
		log.Bh=FALSE,
		bhlink=c("log", "identity"),
		intercept.Bh=TRUE,
		Min_T=0,
		Max_T=NULL,
		model=c("additive","multiplicative"),
		rate, 
		logit_start, 
		logit_end,
		logit_start_byperiod=NULL, 
		logit_end_byperiod=NULL, 
		weights_byperiod=NULL, 
		Id_byperiod=NULL,
		# knots.table and degree.table allow to define the correction model of life table
		contrasts.table = NULL,
		knots.table=c(-2.5,0,2.5),   
		degree.table=3,
		Spline.table=c("restricted B-splines"), # restricted M-splines (linear extrapolation + 2nd derivative at boundaries == 0)
		Spline.CLT=R2bBSplineBasis(knots=c(-2.5,0,2.5), degree=3),
		model_correction = c("cohort", "period"),
		weights=NULL,
		na.action=NULL, 
		datacontrol=NULL,
		Idcontrol,
		ratecontrol, 
		logit_startcontrol, 
		logit_endcontrol,
		logit_start_byperiodcontrol=NULL, 
		logit_end_byperiodcontrol=NULL, 
		weights_byperiodcontrol=NULL, 
		Id_byperiodcontrol=NULL,
		weightscontrol=NULL,
		int_meth=c("GL", "CAV_SIM", "SIM_3_8", "BOOLE", "GLM", "BANDS"),
		bands=NULL,
		npoints=20,              
		stept=NULL,              
		init=NULL,
		optim.control=list(trace=100, REPORT=1, fnscale=-1, maxit=25),
		Coptim.control= list(),
		optim_meth=c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "SANN", "Brent"),
		lower = -Inf,
		upper = Inf,
		vartype = c("oim", "opg", "none"),
		varmethod = c("optim", "numDeriv.hessian", "numDeriv.jacobian"),
		numDeriv.method.args=list(eps=5e-7, d=0.001, zero.tol=sqrt(.Machine$double.eps/7e-4), r=4, v=2),
		debug=FALSE
){
	
	
	debug.gr <- debug.ll <- FALSE
	if( debug > 100000){
		debug.gr <- debug%%1000 + debug%/%100000 * 1000 
	}
	if( debug > 10000){
		debug.ll <- debug%%1000 + (debug%/%10000)%%10 * 1000 
	}
	debug <- debug%%1000 + (debug%/%1000)%%10 * 1000
	
# force optim.control$fnscale to -1 (maximum of the objective function in optim)
	optim.control$fnscale <- -1
	
	# choice of spline basis
	if (!missing(Spline))
		Spline <- match.arg(Spline)
	else {
		Spline <- "b-spline"
	}
	
	if (!missing(bhlink))
		bhlink <- match.arg(bhlink) 
	else {
		bhlink <- "log"
	}
	
	# type of model used (additive for remontet's model - multiplicative for mahboubi's model)
	if (!missing(model))
		model  <- match.arg(model)
	else {
		model <- "additive"
	}  
	
	# type of spline in the correction model 
	if (!missing(Spline.table))
		Spline.table <- match.arg(Spline.table)
	else {
		Spline.table <- "restricted B-splines"
	}
	
# type of correction model
	if (!missing(model_correction))
		model_correction <- match.arg(model_correction)
	else {
		model_correction <- "cohort"
	}
	
	if (!missing(int_meth))
		int_meth <- match.arg(int_meth)
	else {
		int_meth <- "GL"
	}
	
	if (!missing(optim_meth))
		optim_meth <- match.arg(optim_meth) 
	else {
		optim_meth <- "BFGS"
	}
	
# type of variance
	if (!missing(vartype))
		vartype <- match.arg(vartype)
	else {
		vartype <- "oim"
	}
	
	# type of variance
	if (!missing(varmethod))
		varmethod <- match.arg(varmethod)
	else {
		varmethod <- "optim"
	}
	
	
	
	if (int_meth == "BANDS" ){
		int_meth = "GLM"
	}
	
	
# manage na.action on all the variables
	Call <- match.call(expand.dots = FALSE)
	
# na.action not yet implemented
# data are supposed to be NA-free in all the used variables
	mna.action <- match("na.action", names(Call), 0L)
	if (mna.action > 0){
		stop(   gettextf("argument 'na.action' = %s' not yet implemented. Data must be NA-free and remove argument 'na.action' in the call",
						dQuote(Call[["na.action"]]),
						domain=NA))
	}
	
# if (missing(na.action)) {
#     if (!is.null(naa <- getOption("na.action"))) {
#       na.action <- naa
#     }
#   }
	
	
	if( !is.null(na.action)){
		oneformula <-  formula(as.Formula(Call[["formula"]],Call[["formula.table"]]),
				collapse=TRUE)
#          oneformula <-  asOneFormula(Call[["formula"]],
#                                Call[["formula.table"]])
		m <- match(c("formula", "data", "Id", "rate", "logit_start", "logit_end", "weights"), names(Call), 0L)
		newmf <- Call[c(1L, m)]
		newmf$drop.unused.levels <- TRUE
		newmf[[1L]] <- quote(stats::get_all_vars)
		newmf$formula <- oneformula
		mf <- eval(newmf,  sys.parent())
		
		newdata <- na.action( mf )
	} else {
		newdata <- data
	}
	
	
# setting up variables of the excess model
	m <- match(c("formula", "data", "Id", "rate", "logit_start", "logit_end", "weights"), names(Call), 0L)
	if (m[1]==0)
		stop ("The formula argument is required.")
	
	mf <- Call[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- quote(stats::model.frame)
	mf$data <- newdata
	special <- c("NPH","NLL", "NPHNLL", "nl", "td", "nltd", "WCEI") 
	Terms <- if (missing(data)){
				terms(formula, specials=special)
			} else {
				terms(formula, specials=special, data = newdata)
			}
	
	mf$formula <- Terms
	
	
	
	mf <- eval(mf,  sys.parent())
	mt <- attr(mf, "terms")
	
	
	
# detection of WCEI effets
	if( baselinehazard == TRUE & firstWCEIadditive == TRUE){
		stop("baselinehazard = TRUE and firstWCEIadditive == TRUE are incompatible, you must choose only one")
	}
	
	list_var_WCEI <- all_specials_vars(Terms, specials="WCEI",
			unique = FALSE,
			order="formula")
	is_wce_model <- length(list_var_WCEI) > 0
	
	
# switch to additive WCEI for the first WCEI
	if( firstWCEIadditive == TRUE & is_wce_model ){
		is_wce1add_model <- TRUE
		var_WCEI.Bh <-  list_var_WCEI[[1]]
	}
	else {
		is_wce1add_model <- FALSE
	}
	
#    if (is_wce_model) {
#      stop("Weighted cumulative exposure effect (WCEI() not yet implementd")
#    }
	
	
	
	rate <- as.vector(model.extract(mf, "rate"))
	if(!is.null(rate) && !is.numeric(rate))
		stop("'rate' must be a numeric vector")
	if (is.null(rate)){
		stop("'rate' must be specified")
	} else {
		newdata$rate <- rate
	}
	
	
	logit_start <- as.vector(model.extract(mf, "logit_start"))
	
	logit_end <- as.vector(model.extract(mf, "logit_end"))
	
	
	weights <- model.extract(mf, "weights")
	
	# pour l'intercept, voir avec l'option utilisée dans gamma0(t)
	intercept <- attr(mt, "intercept")
	
	Y <- model.extract(mf, "response")
	if (!inherits(Y, "Surv")) {
		stop("Response must be a survival object")
	}
	Survtype <- attr(Y, "type")
	
	if ((ncol(Y) ==  2) ) {
		if (is_wce_model) {
			stop(gettextf("With weighted cumulative exposure effect (WCEI()), you must specify a (start, end] survival object in the left hand side of formula", domaine=NA))
		} else {
			if (Survtype != "right"){
				stop(gettextf("flexrsurvclt.ll does not support %s type of censoring with (0, end] survival data", dQuote(Survtype), domaine=NA))
			} else {
				Y <- cbind(rep(0, dim(Y)[1]), Y)
			}
		}
	} else if ((ncol(Y) ==  3) && (Survtype != "counting") ) {
		stop(gettextf("flexrsurvclt.ll does not support %s type of censoring with (start, end] survival data", dQuote(Survtype), domaine=NA))
	}
	
	
	Id <- as.vector(model.extract(mf, "Id"))
	if (is.null(Id)){
		Id <- 1:(dim(Y)[1])
	} else {
		newdata$Id <- Id
	}
	if(is.null(Id)){
		Id <- 1:(dim(Y)[1])
	}
	nobs <-  length(unique(Id)) 
	
#  if (is_wce_model) {
# add finalT in Y
# get FirstId and Lastid
	Index <- getIndexId(Id)
	FirstId <- Index$FirstId
	LastId <- Index$LastId
	Y <- cbind(Y, Y[LastId, 2])
#    } 
	
	rate <- as.vector(model.extract(mf, "rate"))
	if(!is.null(rate) && !is.numeric(rate))
		stop("'rate' must be a numeric vector")
	if (is.null(rate)){
		stop("'rate' must be specified")
	} else {
		newdata$rate <- rate
	}
	
	
	
	
	weights <- as.vector(model.weights(mf))
	if (!is.null(weights) && !is.numeric(weights)) 
		stop("'weights' must be a numeric vector")
	if (!is.null(weights) && any(weights < 0)) 
		stop("negative weights not allowed")
	
	
	
	
# build Design Matrics
	
# get T-splines
	
	
	if (is.null(Min_T)) {
		Min_T <- max(min(Y[,1], min(bands) ))
	}
	if (is.null(Max_T)) {
		Max_T <- max(c(bands, Y[,2]))
	}
	
	
	if(baselinehazard == TRUE){
		if(Spline=="b-spline"){
			Spline_t0 <- BSplineBasis(knots=c(Min_T, knots.Bh, Max_T),
					degree=degree.Bh,
					keep.duplicates=TRUE,
					log=log.Bh)
			# no log basis for NPH and NPHNLL and td and nltd effects
			Spline_t<- BSplineBasis(knots=c(Min_T, knots.Bh, Max_T),
					degree=degree.Bh,
					keep.duplicates=TRUE,
					log=FALSE)
		} else if(Spline=="tp-spline" || Spline=="tpi-spline"){
# because time is >0, tp-splines are naturaly increasing
			Spline_t0 <-TPSplineBasis(knots=knots.Bh,
					degree=degree.Bh,
					min=Min_T,
					max=Max_T,
					log=log.Bh,
					type="standard")
			
			# no log basis for NPH and NPHNLL  and td and nltd effects
			Spline_t <-TPSplineBasis(knots=knots.Bh,
					degree=degree.Bh,
					min=Min_T,
					max=Max_T,
					log=FALSE,
					type="standard")
		}
	} else {
		Spline_t0 <- NULL
		Spline_t <- NULL
		
	}
	
# lists of variables in the model/formula
	list_var_LIN <- all_LIN_vars(Terms)
	list_var_NLL <- all_specials_vars(Terms, specials="NLL",
			unique = FALSE,
			order="formula")
	list_var_NPH <- all_specials_vars(Terms, specials="NPH",
			unique = FALSE,
			order="formula")
	list_var_NPHNLL <- all_specials_vars(Terms, specials="NPHNLL",
			unique = FALSE,
			order="formula")
	
# coherence between LIN, NLL, NPH, NPHNLL terms
	
	var_LIN_NPHNLL <- list_var_LIN %in% c(list_var_NLL, list_var_NPH, list_var_NPHNLL)
	var_NPH_NLL <- list_var_NPH %in% list_var_NLL
	var_NLL_NPHNLL <- list_var_NLL %in% list_var_NPHNLL
	var_NPH_NPHNLL <- list_var_NPH %in% list_var_NPHNLL
	
#    if( any(var_LIN_NPHNLL)){
#      stop("ERROR: some linear variables are also non-linear or non-proportional variables")
#    }
	
#    if( any(var_NLL_NPHNLL)){
#      stop("ERROR: some non-linear variables are also non-linear-non-proportional variables")
#    }
	
#    if( any(var_NPH_NPHNLL)){
#      stop("ERROR: some non-proportional variables are also non-linear-non-proportional variable")
#    }
	
	
	# build designs matrices
	
	if (!is_wce_model) {
		des <- ReadDesignFlexrsurv(Terms=Terms, modframe=mf, data=newdata, rate=rate, Spline_t0=Spline_t0, intercept.Bh=intercept.Bh )
	} else {
		des <- ReadDesignFlexrsurvWCEI(Terms=Terms, modframe=mf, data=newdata, rate=rate, Spline_t0=Spline_t0, intercept.Bh=intercept.Bh)
	}
	
# test if all time spline are identical
	if(!is.null(des$X) | !is.null(des$Z)){
		allSpline_T <- c(des$Spline_XT, des$Spline_ZT)
		if(length(allSpline_T) > 1){
			for(sb in allSpline_T){
				if(!identical(sb, allSpline_T[[1]])){
					stop("flexrsurvclt.ll cannot handle different spline basis for time dependent effects")
				}
			}
		}
		Spline_t <-allSpline_T[[1]] 
	}
	
	# rate
	therate<-des$rate
	
	
# X0 linear and non linear effects
	X0<-des$X0
	
	
	
# X (NPH effects)
	if(!is.null(des$X)){
		X<-as.matrix(des$X)
		#  Intercept.t of NPH() effets are set in Flersurv:fix.flexrsurv.formula
		#  now build the vector of Intercept_t_NPH for each X var
		# {linear or non linear} and non prop effect are possible; set intercept for NPH effects
		Intercept_t_NPH <- rep(TRUE, length(des$XVars))
		for (i in attr(des$TermsX, "specials")[["NPH"]]){
			thecall <-  match.call(NPH, attr(des$TermsX,"variables")[[i+1]])
			# the variable name is the second argument of the special function
			Intercept_t_NPH[i] <- ifelse( length(thecall[["Intercept.t"]]) == 0 ,
					formals(NPH)[["Intercept.t"]],
					thecall[["Intercept.t"]] )
		}
	} else {
		X <- NULL
		Intercept_t_NPH <- NULL
	}
	
# Z (NPHNLL effects
	
	
	# get splines for each NPHNLL variable
	
	Spline_Z <- des$Spline_Z
	if( !is.null(des$Z) ){
		Z<-DesignMatrixNPHNLL(Z=des$Z, listsplinebasis=Spline_Z, timesplinebasis=Spline_t)
	}  else {
		Z <- NULL
	}
	
	# W, WCEI weightefd cumulative exposure index
	if(!is.null(des$W)){
		W<-as.matrix(des$W)
		Spline_W <- des$Spline_WT
		ISpline_W <- des$Spline_WT
		for( i in dim(W)[2]){
			ISpline_W[[i]] <- integrate(Spline_W[[i]])
		}
		#  Intercept.t of NPH() effets are set in Flersurv:fix.flexrsurv.formula
		#  now build the vector of Intercept_t_NPH for each X var
		# {linear or non linear} and non prop effect are possible; set intercept for NPH effects
		Intercept_W <- rep(TRUE, length(des$WVars))
		for (i in attr(des$TermsW, "specials")[["WCEI"]]){
			thecall <-  match.call(WCEI, attr(des$TermsW,"variables")[[i+1]])
			# the variable name is the second argument of the special function
			Intercept_W[i] <- ifelse( length(thecall[["Intercept.t"]]) == 0 ,
					formals(WCEI)[["Intercept.t"]],
					thecall[["Intercept.t"]] )
			
		}
	} else {
		W <- NULL
		ISpline_W <- NULL
		Intercept_W <- NULL
	}	
	
	# get initvalues
	
	listinit_excess<- list(gamma0 = init[des$coef2param$gamma0], 
			alpha0 = init[des$coef2param$alpha0],
			beta0  = init[des$coef2param$beta0],
			alpha  = init[des$coef2param$alpha],
			beta   = init[des$coef2param$beta],
			eta0 = init[des$coef2param$eta0])
	df_excess <- length(unlist(listinit_excess))
	
	
	
	
# setting up the model of correction of the life table
# if formula.table or spline parameters for CLT is given, we need logit_start end logit_end
	
	mfclt <- match.call(expand.dots = FALSE)
	mclt <- match(c("formula.table", "logit_start", "logit_end", "knots.table", "degree.table", "Spline.table", "Spline.CLT"),
			names(mfclt), 0L)
	
	if(sum(mclt[-c(2,3)]) > 0) {
		is_correction_model <- TRUE
		Intercept_B <- TRUE
		
		# analysis with correction model both model_correction =="period" & "cohorte" 	  
		
		# proportional correction
		# logit_start and logit_end
		if ((mclt[2]==0 | mclt[3]==0)){
			stop ("With correction of life table models, 'logit_start' and 'logit_end' are  required.")
		}
		if(!is.null(logit_start) && !is.numeric(logit_start))
			stop("'logit_start' must be a numeric vector")
		if (is.null(logit_start)){
			stop("'logit_start' must be specified")
		} else {
			newdata$logit_start <- logit_start
		}
		
		if(!is.null(logit_end) && !is.numeric(logit_end))
			stop("'logit_end' must be a numeric vector")
		if (is.null(logit_end)){
			stop("'logit_end' must be specified")
		} else {
			newdata$logit_end <- logit_end
		}
		
# non-proportional correction
		if(sum(mclt[4:7])>0){
			if (mclt[7]!=0){
				Spline_B <- Spline.CLT
			} else {
				if(mclt[4]!=0 & mclt[5]!=0) {
					Spline_B <- R2bBSplineBasis(knots=knots.table, degree=degree.table)
				}
				else {
					stop ("With correction of life table models, both 'knots.table' and 'degree.table' are  required.")
				}
			}
			# df for the brass model (fisrt basis has coef equal to one
			nbrass <- getNBases(Spline_B) - 1 - (1 - Intercept_B)
		}
		else {
			Spline_B <- NULL
			nbrass <- 0
		}
		
		if(model_correction =="period" ){
			#logit_enter_byperiod", "logit_end_byperiod", "weights_byperiod", "Id_byperiod" needed
			
			mcltper <- match(c("logit_start_byperiod", "logit_end_byperiod", "weights_byperiod", "Id_byperiod"),
					names(mfclt), 0L)
			if(sum(mcltper) == 0L){
				stop ("With period-type correction of life table models, both 'logit_start_byperiod', 'logit_end_byperiod', 'weights_byperiod', 'Id_byperiod' are  required.")
			}
			
			if (is.null(logit_start_byperiod)){
				stop("'logit_start_byperiod' must be specified")
			} else if( !is.numeric(logit_start_byperiod)){
				stop("'logit_start_byperiod' must be a numeric vector")
			}
			
			
			if (is.null(logit_end_byperiod)){
				stop("'logit_end_byperiod' must be specified")
			} else if(!is.numeric(logit_end_byperiod)){
				stop("'logit_end_byperiod' must be a numeric vector")
			}
		} 
		
#design matrix of the proportional terms of the correction model 
		if (mclt[1]!=0 ){
			mclt2 <- match(c("formula.table", "data"),
					names(mfclt), 0L)
			mfclt <- mfclt[c(1L, mclt2)]
			mfclt$drop.unused.levels <- TRUE
			mfclt[[1L]] <- quote(stats::model.frame)
			names(mfclt)[[2L]] <- "formula"
			mfclt$data <- newdata
			
			# NLL effects are allowed
			# but it is unnecessary to consider NLL effect as "special"
			# because stats::model.frame can manage NLL()
#    special <- c("NLL", "nl") 
#    Terms <- if (missing(data)){
#      terms(formula, specials=special)
#    } else {
#      terms(formula, specials=special, data = data)
#    }
			
			mfclt <- eval(mfclt,  sys.parent())
			mtclt <- attr(mfclt, "terms")
			BX0 <- model.matrix(mtclt, mfclt, contrasts.table)
			nBX0 <- dim(BX0)[2]
		}
		else { 
			# no BX0 now
			#    matrix(numeric(0), ncol=0)
			BX0 <- NULL
			nBX0 <- 0
		}
		
		
		listinit_table<- list(brass0 = if(is.null(init)){
							# last coef is the slope of the right asymptote
							c(rep(0, nbrass-1), 1.0)
						} else {
							init[1:nbrass + df_excess]
						},
				balpha0 = if (is.null(BX0)) {
							NULL
						} else {
							if(is.null(init)){
								rep(0, nBX0)
							} else {
								init[nbrass + df_excess + 1:nBX0]
							}
						}
		)
		
	} else {
# no correction of life table
		is_correction_model <- FALSE
		model_correction <- "cohort"
		Spline_B <- NULL
		Intercept_B <- TRUE
		nbrass <- 0
		BX0 <- NULL
		nBX0 <- 0
		listinit_table <- NULL
		if(is.null(logit_end)){
			logit_end  <- rep(0, dim(Y)[1])
		}
		if(is.null(logit_start)){
			logit_start <- rep(0, dim(Y)[1])
		}
		if(is.null(logit_end)){
			logit_end  <- rep(0, dim(Y)[1])
		}
		if(is.null(logit_end_byperiod)){
			logit_end_byperiod <- rep(0, dim(Y)[1])
			Id_byperiod <- 1:(dim(Y)[1])
		}          
		if(is.null(logit_start_byperiod)){
			logit_start_byperiod <- rep(0, dim(Y)[1])
		}                	  
	}
	
	listinit <- c(listinit_excess, listinit_table)
	
	
	
# control Group
	if(!is_correction_model | is.null(datacontrol)){
		# no control group
		Ycontrol <- NULL
		BX0control <- NULL 
		weightscontrol <- NULL
		Idcontrol <- NULL
		FirstIdcontrol <- NULL
		rate_control <- NULL
		logit_endcontrol <- NULL
		logit_startcontrol <- NULL
		logit_end_byperiodcontrol <- NULL
		logit_start_byperiodcontrol <- NULL
		Id_byperiodcontrol <- NULL
		weights_byperiodcontrol <- NULL
		nobscontrol <- 0
		neventcontrol <- 0
		ncontrol <- 0
	} else {
		# we have a control dataset & it is a correction modle
		
		
		# formulacontrol is formula without rhs
		
		# chck that formula.table is onesided formula
		
		formulacontrol <- formula
		rhs(formulacontrol)<-quote(1)
#  formulacontrol[[3]]<-quote(1)
		if(!is.null(formula.table)){
#  formulacontrol[[3]] <- formula.table[[2]]
			formula.table <- asOneSidedFormula(formula.table)
			rhs(formulacontrol)<-rhs(formula.table)
		}
		
		mfc <- match.call(expand.dots = FALSE)
		indx <- match(c("formula", "datacontrol", "Idcontrol",
						"ratecontrol", "logit_startcontrol", "logit_endcontrol",
						"weightscontrol", "na.action"), names(mfc), 0L)
		mfc <- mfc[c(1, indx)]
		
		mfc$drop.unused.levels <- TRUE
		mfc[[1L]] <- quote(stats::model.frame)
		names(mfc)[3] <- "data"
		Termsc <- if (missing(datacontrol)){
					terms(formulacontrol)
				} else {
					terms(formulacontrol, data = datacontrol)
				}
		
		mfc$formula <- Termsc
		
		mfc <- eval(mfc, parent.frame())
		if (nrow(mfc) == 0) {
			# no control group
			Ycontrol <- NULL
			BX0control <- NULL 
			weightscontrol <- NULL
			Idcontrol <- NULL
			FirstIdcontrol <- NULL
			rate_control <- NULL
			logit_endcontrol <- NULL
			logit_startcontrol <- NULL
			logit_end_byperiodcontrol <- NULL
			logit_start_byperiodcontrol <- NULL
			Id_byperiodcontrol <- NULL
			weights_byperiodcontrol <- NULL
			nobscontrol <- 0
			neventcontrol <- 0
			ncontrol <- 0
			warning("No (non-missing) controls")
		}
		else {
			mtc <- terms(mfc)
			Ycontrol <- model.extract(mfc, "response")
			if (!inherits(Ycontrol, "Surv")) 
				stop("Response must be a survival object")
			type <- attr(Ycontrol, "type")
			if (type != "right" && type != "counting") 
				stop(paste("Cox model doesn't support \"", type, "\" survival data", 
								sep = ""))
			data.n <- nrow(Ycontrol)
			
			#design matrix of the proportional terms of the correction model 
			# BX0control
			if(match(c("formula.table"), names(mfc), 0L)){
				xlevels <- .getXlevels(mtc, mfc)
				contrast.arg <- NULL
				BX0control <- model.matrix(mtc, mfc, contrasts = contrast.arg)
				nBX0control <- dim(BX0)[2]
			}
			else {
				BX0control <- NULL
				nBX0control <- 0
			}
			
			if(indx[3]>0L){
				Idcontrol <- as.vector(model.extract(mfc, "Idcontrol"))
				if (is.null(Idcontrol)){
					Idcontrol <- 1:(dim(Ycontrol)[1])
				}
			} else {
				Idcontrol <- 1:(dim(Ycontrol)[1])
			}
			nobscontrol <- length(unique(Idcontrol))
			neventcontrol <- sum(Ycontrol[,ncol(Ycontrol)])
			ncontrol <- nrow(Ycontrol)
			
			if(indx[4]!=0L){
				rate_control <- as.vector(model.extract(mfc, "ratecontrol"))
				if(!is.null(rate_control) && !is.numeric(rate_control))
					stop("'ratecontrol' must be a numeric vector")
				if (is.null(rate_control)){
					stop("'ratecontrol' must be specified")
				}
			} else {
				stop("'ratecontrol' must be specified")
			}
			
			
			if(indx[5]!=0L){
				logit_startcontrol <- as.vector(model.extract(mfc, "logit_startcontrol"))
				if(!is.null(logit_startcontrol) && !is.numeric(logit_startcontrol))
					stop("'logit_startcontrol' must be a numeric vector")
				if (is.null(logit_startcontrol)){
					stop("'logit_startcontrol' must be specified")
				}
			} else {
				stop("'logit_startcontrol' must be specified")
			}
			
			
			if(indx[6]!=0L){
				logit_endcontrol <- as.vector(model.extract(mfc, "logit_endcontrol"))
				if(!is.null(logit_endcontrol) && !is.numeric(logit_endcontrol))
					stop("'logit_endcontrol' must be a numeric vector")
				if (is.null(logit_endcontrol)){
					stop("'logit_endcontrol' must be specified")
				}
			} else {
				stop("'logit_endcontrol' must be specified")
			}
			
#  if (is_wce_model) {
# add finalT in Y
# get FirstId and Lastid
			Indexcontrol <- getIndexId(Idcontrol)
			FirstIdcontrol <- Indexcontrol$FirstId
			LastIdcontrol <- Indexcontrol$LastId
			Ycontrol <- cbind(Ycontrol, Ycontrol[LastIdcontrol, 2])
			
			weightscontrol <- as.vector(model.weights(mfc))
			if (!is.null(weightscontrol) && !is.numeric(weightscontrol)) 
				stop("'weights' must be a numeric vector")
			if (!is.null(weightscontrol) && any(weightscontrol < 0)) 
				stop("negative weights not allowed")
			
			if(model_correction =="period" ){
				mcltper <- match(c("logit_enter_byperiodcontrol", "logit_end_byperiodcontrol", "weights_byperiodcontrol", "Id_byperiodcontrol"),
						names(mfclt), 0L)
				if(sum(mcltper) == 0L){
					stop ("With period-type correction of life table models, both 'logit_enter_byperiodcontrol', 'logit_end_byperiodcontrol', 'weights_byperiodcontrol', 'Id_byperiodcontrol' are  required.")
				}
				
				if (is.null(logit_start_byperiodcontrol)){
					stop("'logit_start_byperiodcontrol' must be specified")
				} else if( !is.numeric(logit_start_byperiodcontrol)){
					stop("'logit_start_byperiodcontrol' must be a numeric vector")
				}
				
				
				if (is.null(logit_end_byperiodcontrol)){
					stop("'logit_end_byperiodcontrol' must be specified")
				} else if(!is.numeric(logit_end_byperiodcontrol)){
					stop("'logit_end_byperiodcontrol' must be a numeric vector")
				}
			}
			
		}
		
		
		
		
#  # remove NAs
# if( !is.null(na.action)){
#    m <- match(c("formula", "datacontrol", "Idcontrol", "ratecontrol", "logit_startcontrol", "logit_endcontrol", "weightscontrol"), names(Call), 0L)
#    newmf <- Call[c(1L, m)]
#    newmf$drop.unused.levels <- TRUE
#    newmf[[1L]] <- quote(stats::get_all_vars)
#    newmf$formula <- formulacontrol
#    mf <- eval(mf,  sys.parent())
# 
#    newdatacontrol <- na.action( mf )
#  } else {
#    newdatacontrol <- data
#  }
		
	}
	
	
	
# get methods
	
	if( int_meth=="GLM" ){        
		if (is.null(bands)){
			bands <- default_bands(Spline_t)
		}
		method <- list(int_meth=int_meth, bands=bands, optim_meth=optim_meth, constOptim_meth=optim_meth, lower=lower, upper=upper)
	} else if( int_meth == "GL"){
		if (is.null(npoints)){
			npoints <- 20
		}
		method <- list(int_meth="Gauss-Legendre", npoints=npoints, optim_meth=optim_meth, constOptim_meth=optim_meth, lower=lower, upper=upper)
	} else {
		if (is.null(stept)){
			stept <- Max_T /500
		}
		method <- list(int_meth=int_meth, stept=stept, optim_meth=optim_meth, constOptim_meth=optim_meth, lower=lower, upper=upper)
	}
	
	# isEnter == 1 only when first row AND Tenter > 0
	
# fit
	if(is_wce1add_model == TRUE){
		#  case of one time-varying exposure variable use as baseline 
		if(model_correction =="period" ){
			fit<-flexrsurv.ll.fromto.brass0.period.wceadd.fitCoptim(X0=X0, X=X, Z=Z, Y=Y, W=W,
					BX0=BX0,
					Id=Id, FirstId=FirstId, 
					expected_rate = rate,
					expected_logit_end = logit_end,
					expected_logit_enter = logit_start, 
					expected_logit_end_byperiod=logit_end_byperiod, 
					expected_logit_enter_byperiod=logit_start_byperiod, 
					weights_byperiod=weights_byperiod, 
					Id_byperiod=Id_byperiod,
					weights=weights,
					Ycontrol=Ycontrol, BX0control=BX0control, 
					weightscontrol=weightscontrol,
					Idcontrol=Idcontrol, FirstIdcontrol=FirstIdcontrol,
					expected_ratecontrol = rate_control,
					expected_logit_endcontrol = logit_endcontrol,
					expected_logit_entercontrol = logit_startcontrol,
					expected_logit_end_byperiodcontrol=logit_end_byperiodcontrol, 
					expected_logit_enter_byperiodcontrol=logit_start_byperiodcontrol, 
					weights_byperiodcontrol=weights_byperiodcontrol, 
					Id_byperiodcontrol=Id_byperiodcontrol,
					Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
					ISpline_W = ISpline_W[[1]], Intercept_W=Intercept_W,
					Spline_B =Spline_B, Intercept_B=Intercept_B,
					init=listinit,
					optim.control=optim.control,
					Coptim.control=Coptim.control,
					method=method,
					vartype = vartype,
					varmethod = varmethod,
					numDeriv.method.args=numDeriv.method.args,
					namebrass="CorTable",
					debug=debug,
					debug.ll=debug.ll, debug.gr=debug.gr)
		} else if(model_correction =="cohort" ){
			fit<-flexrsurv.ll.fromto.brass0.wceadd.fitCoptim(X0=X0, X=X, Z=Z, Y=Y, W=W,
					BX0=BX0,
					Id=Id, FirstId=FirstId, 
					expected_rate = rate,
					expected_logit_end = logit_end,
					expected_logit_enter = logit_start, 
					weights=weights,
					Ycontrol=Ycontrol, BX0control=BX0control, 
					weightscontrol=weightscontrol,
					Idcontrol=Idcontrol, FirstIdcontrol=FirstIdcontrol,
					expected_ratecontrol = rate_control,
					expected_logit_endcontrol = logit_endcontrol,
					expected_logit_entercontrol = logit_startcontrol,
					Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
					ISpline_W = ISpline_W[[1]], Intercept_W=Intercept_W,
					Spline_B =Spline_B, Intercept_B=Intercept_B,
					init=listinit,
					optim.control=optim.control,
					Coptim.control=Coptim.control,
					method=method,
					vartype = vartype,
					varmethod = varmethod,
					numDeriv.method.args=numDeriv.method.args,
					namebrass="CorTable",
					debug=debug,
					debug.ll=debug.ll, debug.gr=debug.gr)
		}
	}
	else {
		if(model_correction =="period" ){
			fit<-flexrsurv.ll.fromto.brass0.period.wce.fitCoptim(X0=X0, X=X, Z=Z, Y=Y, W=W,
					BX0=BX0,
					Id=Id, FirstId=FirstId, 
					expected_rate = rate,
					expected_logit_end = logit_end,
					expected_logit_enter = logit_start, 
					expected_logit_end_byperiod=logit_end_byperiod, 
					expected_logit_enter_byperiod=logit_start_byperiod, 
					weights_byperiod=weights_byperiod, 
					Id_byperiod=Id_byperiod,
					weights=weights,
					Ycontrol=Ycontrol, BX0control=BX0control, 
					weightscontrol=weightscontrol,
					Idcontrol=Idcontrol, FirstIdcontrol=FirstIdcontrol,
					expected_ratecontrol = rate_control,
					expected_logit_endcontrol = logit_endcontrol,
					expected_logit_entercontrol = logit_startcontrol,
					expected_logit_end_byperiodcontrol=logit_end_byperiodcontrol, 
					expected_logit_enter_byperiodcontrol=logit_start_byperiodcontrol, 
					weights_byperiodcontrol=weights_byperiodcontrol, 
					Id_byperiodcontrol=Id_byperiodcontrol,
					Spline_t0=Spline_t0, Intercept_t0=intercept.Bh,
					Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
					ISpline_W = ISpline_W, Intercept_W=Intercept_W,
					Spline_B =Spline_B, Intercept_B=Intercept_B,
					bhlink=bhlink,
					init=listinit,
					optim.control=optim.control,
					Coptim.control=Coptim.control,
					method=method,
					vartype = vartype,
					varmethod = varmethod,
					numDeriv.method.args=numDeriv.method.args,
					namebrass="CorTable",
					debug=debug,
					debug.ll=debug.ll, debug.gr=debug.gr)
		} else if(model_correction =="cohort" ){
			fit<-flexrsurv.ll.fromto.brass0.wce.fitCoptim(X0=X0, X=X, Z=Z, Y=Y, W=W,
					BX0=BX0,
					Id=Id, FirstId=FirstId, 
					expected_rate = rate,
					expected_logit_end = logit_end,
					expected_logit_enter = logit_start, 
					weights=weights,
					Ycontrol=Ycontrol, BX0control=BX0control, 
					weightscontrol=weightscontrol,
					Idcontrol=Idcontrol, FirstIdcontrol=FirstIdcontrol,
					expected_ratecontrol = rate_control,
					expected_logit_endcontrol = logit_endcontrol,
					expected_logit_entercontrol = logit_startcontrol,
					Spline_t0=Spline_t0, Intercept_t0=intercept.Bh,
					Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
					ISpline_W = ISpline_W, Intercept_W=Intercept_W,
					Spline_B =Spline_B, Intercept_B=Intercept_B,
					bhlink=bhlink,
					init=listinit,
					optim.control=optim.control,
					Coptim.control=Coptim.control,
					method=method,
					vartype = vartype,
					varmethod = varmethod,
					numDeriv.method.args=numDeriv.method.args,
					namebrass="CorTable",
					debug=debug,
					debug.ll=debug.ll, debug.gr=debug.gr)
			
		}
	}
	
	
	
# returned value  
# buid an object with attributes similares to glm object.
	objfit <- fit
	
#reorder coefs
	objfit$coefficients <- fit$coefficients[des$param2coef]
	if(length(des$names_coef_XZ)){
		if(baselinehazard == TRUE){
			names(objfit$coefficients)[(des$df.T0+1):df_excess] <- des$names_coef_XZ
		} else {
			names(objfit$coefficients[1:df_excess]) <- des$names_coef_XZ
		}
	}
	reordercoef <- des$param2coef
	if(is_correction_model){
		objfit$coefficients <- c(objfit$coefficients, fit$coefficients[-(1:df_excess)])
		reordercoef <- c(des$param2coef, (df_excess+1):length(fit$coefficients))
	}
	
	if( !is.null(dim(fit$var))){
#      objfit$var <- (fit$var[des$param2coef, ])[,des$param2coef]  
		objfit$var <- (fit$var[reordercoef, ])[,reordercoef]  
		dimnames(objfit$var)[[1]] <- names(objfit$coefficients)
		dimnames(objfit$var)[[2]] <- names(objfit$coefficients)
		attr(objfit$var, "type") <- vartype
	}
	
	if( !is.null(dim(fit$informationMatrix))){
#      objfit$informationMatrix <- (fit$informationMatrix[des$param2coef, ])[,des$param2coef]  
		objfit$informationMatrix <- (fit$informationMatrix[reordercoef, ])[,reordercoef]  
		dimnames(objfit$informationMatrix)[[1]] <- names(objfit$coefficients)
		dimnames(objfit$informationMatrix)[[2]] <- names(objfit$coefficients)
		attr(objfit$informationMatrix, "type") <- vartype
		attr(objfit$informationMatrix, "method") <- varmethod
	}
	
	objfit$des <- des
	objfit$terms <- Terms
	objfit$mt <- mt
# don't know if it is necessary?
#  attr(objfit$terms, ".Environment") <- parent.frame()
	objfit$assign <- des$assign
	objfit$assignList <- des$assignList
	
	
	if(!is.null(logit_start)){
		objfit$logit_start <- logit_start
	} else {
		objfit$logit_start <- NULL
	}
	if(!is.null(logit_start)){
		objfit$logit_end <- logit_end
	} else {
		objfit$logit_end <- NULL
	}
	
	# Affichage de la formule initiale et de la formule
	#----------------------------------------------------------------------------------------
	objfit[["call"]] <- Call
	objfit[["formula"]] <- formula
	if(dim(Y)[2]>=3){
		objfit[["entertime"]] <-  Y[, 1]
		objfit[["time"]] <-  Y[, 2]
	}
	objfit[["time"]] <-  Y[, 1]
	objfit[["workingformula"]] <- formula
	
	
	
	# Simplified names of coefficients
	#----------------------------------------------------------------------------------------
	
	names(objfit$coefficients) <- make.shortnames.coefficients(names(objfit$coefficients),
			formula=objfit[["workingformula"]] ,
			model=model,
			Spline=Spline,
			baselinehazard=baselinehazard,
			firstWCEIadditive=firstWCEIadditive,
			knots.Bh=knots.Bh,
			degree.Bh=degree.Bh,
			intercept.Bh=intercept.Bh,
			log.Bh=log.Bh)
	
	objfit$baselinehazard <- baselinehazard
	objfit$firstWCEIadditive <- firstWCEIadditive 
	
	objfit$na.action <- attr(m, "na.action")
	
	objfit$optim.control <- optim.control
	
	objfit$converged <- objfit$conv
	objfit$conv <- NULL

#number of events
	objfit$neventexposed <- sum(Y[,ncol(Y)])
	objfit$neventcontrol <- neventcontrol
	objfit$nevent <- objfit$neventexposed + objfit$neventcontrol

# number of lines
	objfit$nexposed <- nrow(Y)
	objfit$ncontrol <- ncontrol
	objfit$n <- objfit$nexposed + objfit$ncontrol

# number of obs
	objfit$nobsexposed <- nobs
	objfit$nobscontrol <- nobscontrol
	objfit$nobs <- objfit$nobsexposed + objfit$nobscontrol

	if(is.null(Spline_t0)){
		class(objfit) <- c("flexrsurvclt2.mle", "flexrsurv.mle")
	}
	else {
		class(objfit) <- c("flexrsurvclt.mle", "flexrsurv.mle")
	}
	return(objfit)
	
}

