# fit of a flexible relative survival model with a model of correction of life table 

flexrsurvclt <- function(formula=formula(data),
		formula.table=NULL, 
		data=parent.frame(),
		Id,
		# knots.Bh and degree.Bh allow to define the Baseline Hazard
		baselinehazard=TRUE,
		firstWCEIadditive=FALSE,
		knots.Bh,   
		degree.Bh=3,
		intercept.Bh=TRUE,
		Spline=c("b-spline", "tp-spline", "tpi-spline"), # tp-spline for truncated power basis
		log.Bh=FALSE,
		bhlink=c("log", "identity"),
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
		contrasts.table = NULL,
		# knots.table and degree.table allow to define the correction model of lifi table
		knots.table=c(-2.5,0,2.5),   
		degree.table=3,
		Spline.table=c("restricted B-splines"), # restricted B-splines (linear extrapolation + 2nd derivative at boundaries == 0)
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
		initbyglm=TRUE,
		initbands=bands,
		optim.control=list(trace=100, REPORT=1, fnscale=-1, maxit=25), 
		optim_meth=c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "SANN", "Brent"),
		Coptim.control=list(),
		lower = -Inf,
		upper = Inf,
		control.glm=list(epsilon = 1e-8, maxit = 100, trace = FALSE, epsilon.glm = 1e-1, maxit.glm = 25 ),
		vartype = c("oim", "opg", "none"),
		varmethod = c("optim", "numDeriv.hessian", "numDeriv.jacobian"),
		numDeriv.method.args=list(eps=5e-7, d=0.001, zero.tol=sqrt(.Machine$double.eps/7e-4), r=4, v=2),
		debug=FALSE
){
	
	
	
	
	debug.ll <- FALSE
	debug.gr <- FALSE
	
	# formula: for example Surv(time,cens)~sex
	# data: the observed data set
	# rate: rate variable in data
	
	
	thecall   <- match.call()
	
	
	## extract x, y, etc from the model formula and frame
	if(missing(data)) data <- environment(formula)
	
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
	
	
	# check options
	#----------------------------------------------------------------------------------------
	
	if (int_meth == "BANDS"){
		if( is.null(bands) ) {
			stop(gettextf("argument 'bands' must be specified if 'int_meth' = %s.", dQuote("BANDS"), domain=NA))
		}
	} else if (int_meth == "GL"){
		if( is.null(npoints) ) {
			stop(gettextf("argument 'npoints' must be specified if 'int_meth' = %s.", dQuote("GL"), domain=NA))
		}
	} else {
		if( is.null(stept) ) {
			stop(gettextf("argument 'step' must be specified if 'int_meth' != %s or != %s.", dQuote("BANDS"), dQuote("GL"), domain=NA))
		}
	}
	
	
	mf <- match.call(expand.dots = FALSE)
	
	m <- match(c("formula", "data", "rate", "weights", "logit_start", "logit_end", "Id"), names(mf), 0L)
	if (m[1]==0)
		stop ("'formula' is required.")
	
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- quote(stats::model.frame)
	specials <- c("NPH","NLL", "NPHNLL", "nl", "td", "nltd", "WCEI") 
	
	Terms <- if (missing(data)){ 
				terms(formula, specials=specials)
			} else {
				terms(formula, data=data, specials=specials)
			}
	
	mf$formula <- Terms
	
	mf <- eval(mf, parent.frame())
	
	if (nrow(mf) == 0) 
		stop("No (non-missing) observations")    
	
# detection of WCEI effets
	list_var_WCEI <- all_specials_vars(Terms, specials="WCEI",
			unique = FALSE,
			order="formula")
	is_wce_model <- length(list_var_WCEI) > 0
	
	
	
	Id <- model.extract(mf, "Id")
	Y <- model.extract(mf, "response")
	if (!inherits(Y, "Surv")) 
		stop("Response must be a survival object.")
	
	Survtype <- attr(Y, "type")
	
	if ((ncol(Y) ==  2) && (Survtype != "right") ) {
		stop(gettextf("flexrsurvclt does not support %s type of censoring with (0, end] survival data", dQuote(Survtype), domaine=NA))
	} else if ((ncol(Y) ==  3) && (Survtype != "counting") ) {
		stop(gettextf("flexrsurvclt does not support %s type of censoring with (start, end] survival data", dQuote(Survtype), domaine=NA))
	}
	
	if ((ncol(Y) ==  2) ) {
		if (is_wce_model) {
			stop(gettextf("With weighted cumulative exposure effect (WCE()), you must specify a (start, end] survival object in the left hand side of formula", domaine=NA))
		} else {
			if (Survtype != "right"){
				stop(gettextf("flexrsurvclt does not support %s type of censoring with (0, end] survival data", dQuote(Survtype), domaine=NA))
			} else {
				Y <- cbind(rep(0, dim(Y)[1]), Y)
			}
		}
	} else if ((ncol(Y) ==  3) && (Survtype != "counting") ) {
		stop(gettextf("flexrsurvclt does not support %s type of censoring with (start, end] survival data", dQuote(Survtype), domaine=NA))
	}
	
	
	entervar <- Y[, 1]   # enter time variable  
	timevar <- Y[, 2]   # time variable  
	fail    <- Y[, 3]   # dummy 1=event, 0=censur
	
# add finalT in Y
	Y <- cbind(Y, Y[Id, 2])
# get FirstId and Lastid
	Index <- getIndexId(Id)
	
	FirstId <- Index$FirstId
	LastId <- Index$LastId
	
	
	
	
	theweights <- as.vector(model.weights(mf))
	if (!is.null(theweights) && !is.numeric(theweights)) 
		stop("'weights' must be a numeric vector")
	if (!is.null(theweights) && any(theweights < 0)) 
		stop("negative weights not allowed")
	
	
	
	if (is.null(Min_T)) {
		Min_T <- max(entervar, min(bands) )
	}
	if (is.null(Max_T)) {
		Max_T <- max(c(bands, timevar))
	}
	
	
	weights <- as.vector(model.weights(mf))
	if (!is.null(weights) && !is.numeric(weights)) 
		stop("'weights' must be a numeric vector")
	if (!is.null(weights) && any(weights < 0)) 
		stop("negative weights not allowed")
	
	rate <- as.vector(model.extract(mf, "rate"))
	if(!is.null(rate) && !is.numeric(rate))
		stop("'rate' must be a numeric vector")
	if (is.null(rate)){
		stop("'rate' must be specified")
	} else {
		data$rate <- rate
	}
	
	
#  logit_start <- as.vector(model.extract(mf, "logit_start"))
#  if(!is.null(logit_start) && !is.numeric(logit_start))
#      stop("'logit_start' must be a numeric vector")
#  if (is.null(logit_start)){
#    stop("'logit_start' must be specified")
#  } else {
#    data$logit_start <- logit_start
#  }
	
#  logit_end <- as.vector(model.extract(mf, "logit_end"))
#  if(!is.null(logit_end) && !is.numeric(logit_end))
#      stop("'logit_end' must be a numeric vector")
#  if (is.null(logit_end)){
#    stop("'logit_end' must be specified")
#  } else {
#    data$logit_end <- logit_end
#  }
	
	
	thedata <- data
	
	
	
	if (baselinehazard == TRUE & is.null(knots.Bh)) {
		stop("The knots for the baseline hazard must be specified in 'knots.Bh'.") 
	}
	
	#----------------------------------------------------------------------------------------
	# working formulas 
	#----------------------------------------------------------------------------------------
	
	# fix the formula
	# pass model and Spline parameter
	# check consistantcy of the formula
	# NLL(x) + NPH(X, t) or nl(x) +td(X,t)
	#----------------------------------------------------------------------------------------
	fixformula <- fix.flexrsurv.formula(formula=as.formula(formula), data=data, Spline=Spline, model=model, debug=debug)
	termsf <- terms(fixformula, specials=c("NLL","NPH","NPHNLL", "nl", "td", "nltd", "WCEI"))
	
	formulamle <- make.mle.formula(formula=fixformula, data= data)
	
	# Init values
	#----------------------------------------------------------------------------------------
	if (initbyglm == TRUE) {
		message("init values are searched with uncorrected life table")
		
		if (int_meth != "BANDS" & is.null(initbands)){
			initbands <- seq(Min_T, Max_T, stept*5)
			message("'initbands' has been set to seq(Min_T, Max_T, stept*5)")
		}
		
		
		name.runningtime <- ".t"
#  # build the formula for the GLM with the splited dataset
#  #----------------------------------------------------------------------------------------
		formulasplit <- make.glm.formula(formula=fixformula, data= data, name.runningtime=name.runningtime,
				Min_T=Min_T,Max_T=Max_T, model=model)
		splitdonnees <- split.data(jeudata=data, bands=initbands, fail=fail, entry=entervar, exit=timevar, name.runningtime=name.runningtime)
		if ((model=="additive") | ((model=="multiplicative") & (length(c(attr(termsf, "specials")$NPHNLL, attr(termsf, "specials")$nltp))==0)))  {
			control.glm <- do.call("control.flexrsurv.additive", control.glm)
			results <- flexrsurv.glm.fit(formula=formulasplit, data=splitdonnees, model=model,
					Spline=Spline, baselinehazard = baselinehazard,
					degree.Bh=degree.Bh, knots.Bh=knots.Bh, intercept.Bh=intercept.Bh, log.Bh=log.Bh,
					control=control.glm,
					Min_T=Min_T, Max_T=Max_T, name.runningtime=name.runningtime, start=init)
			
			
		} else if((model=="multiplicative") &  (length(c(attr(termsf, "specials")$NPHNLL, attr(termsf, "specials")$nltp))!=0)) {
			
			control.glm <- do.call("control.flexrsurv.multiplicative", control.glm)
			results <- flexrsurv.glmiterative.fit(formula=formulasplit, data=splitdonnees, model=model, 
					Spline=Spline, baselinehazard = baselinehazard,
					degree.Bh=degree.Bh, knots.Bh=knots.Bh, intercept.Bh=intercept.Bh, log.Bh=log.Bh, 
					Min_T=Min_T, Max_T=Max_T, name.runningtime=name.runningtime, start=init,
					control=control.glm)
			
		}
		oldinit<-init
		init <- results$coef
		names(init) <- NULL
		if ( ( baselinehazard == TRUE) & ( bhlink == "identity") ){
			init[1:(length(knots.Bh)+degree.Bh+intercept.Bh)] <- exp(init[1:(length(knots.Bh)+degree.Bh+intercept.Bh)]) 
		}
		if( !debug ){
			results$fit <- NULL
		}
		if( !debug  ){
			rm(results, splitdonnees)
		}
		else {
			glmresults <- results
			rm(results)
		}
		
	} 
	# end of initbyglm
	
	# final Estimation,with variance computation
	#----------------------------------------------------------------------------------------
	call.ll <- match.call(definition=flexrsurvclt, expand.dots=FALSE)
	m <- match(names(formals(flexrsurvclt.ll)), names(call.ll), 0L)
	
	# keep flexrsurv.ll() args
	call.ll <- call.ll[c(1L, m)]
	call.ll[[1L]] <- quote(flexrsurvclt.ll)
	call.ll$formula <- formulamle
	
# get the the number of degree of freedom of the excess part and the corection part
	
	call.ndf <- call.ll  
	# ::: because ndf.flexrsurvclt is not exported
#  call.ndf[[1L]] <- quote(flexrsurv:::ndf.flexrsurvclt)
	call.ndf[[1L]] <- quote(ndf.flexrsurvclt)
	
	ndf <- eval(call.ndf, parent.frame())
# update init values when initbyglm == TRUE
	if (initbyglm == TRUE) {
		# rem: there is stil a pb in case of factor variables in the linear effects
		call.ll$init <- init[!is.na(init)]
		if(ndf[["list_df"]][["brass0"]]){
			initbrass <- rep(1, ndf[["list_df"]][["brass0"]])
			initbrass[1] <- 0.0
			call.ll$init <- c(call.ll$init, initbrass)
		}
		if(ndf[["list_df"]][["balpha0"]]){
			call.ll$init <- c(call.ll$init, rep(1, ndf[["list_df"]][["balpha0"]]))
		}
		
	}
	
	
	results <- eval(call.ll, parent.frame()) 
	
	if( !debug ){
		results$fit <- NULL
	}
	
	
	if( debug & initbyglm==TRUE ){
		results$glmfit <- glmresults
		results$splitdonnees <- splitdonnees
	}
	
	
	
	# Affichage de la formule initiale et de la formule
	#----------------------------------------------------------------------------------------
	results[["call"]] <- thecall
	results[["formula"]] <- fixformula
	if(dim(Y)[2]>=3){
		results[["entertime"]] <- entervar
	}
	results[["time"]] <- timevar
	results[["workingformula"]] <- formulamle
	
	
	
	# Simplified names of coefficients
	#----------------------------------------------------------------------------------------
	
#	names(results$coefficients) <- make.shortnames.coefficients(names(results$coefficients),
#			formula=results[["workingformula"]] ,
#			model=model,
#			Spline=Spline,
#			baselinehazard=baselinehazard,
#			firstWCEIadditive=firstWCEIadditive,
#			knots.Bh=knots.Bh,
#			degree.Bh=degree.Bh,
#			intercept.Bh=intercept.Bh,
#			log.Bh=log.Bh)
	
	
	
	if( !is.null(dim(results$var))){
		dimnames(results$var)[[1]] <-  names(results$coefficients)
		dimnames(results$var)[[2]] <-  names(results$coefficients)
	}
	if( !is.null(dim(results$informationMatrix))){
		dimnames(results$informationMatrix)[[1]] <-  names(results$coefficients)
		dimnames(results$informationMatrix)[[2]] <-  names(results$coefficients)
	}
	
	results$baselinehazard <- baselinehazard
	results$firstWCEIadditive <- firstWCEIadditive 
	
	
	results$data <- thedata
	results$rate <- rate
	
	results$ndf <- ndf
	
	results$terms <- terms(formulamle, special=c("NLL", "NPH", "NPHNLL", "td", "nltd"),
			data = data)
	
	if (initbyglm==TRUE) {
		results$init<- oldinit
		results$control.glm <- control.glm
	} else {
		results$init<- init
	}
	
	class(results) <- c("flexrsurvclt", "flexrsurv", class(results))
	
	
	return(results)
	
}



