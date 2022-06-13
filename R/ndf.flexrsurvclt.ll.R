# get the number of degree of freed of the excess part and the corection part

ndf.flexrsurvclt <- function(formula=formula(data),
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
		contrasts.table = NULL,
		# knots.table and degree.table allow to define the correction model of lifi table
		knots.table=c(-2.5,0,2.5),   
		degree.table=3,
		Spline.table=c("restricted B-splines"), # restricted M-splines (linear extrapolation + 2nd derivative at boundaries == 0)
		Spline.CLT=R2bBSplineBasis(knots=c(-2.5,0,2.5), degree=3),
		weights=NULL,
		na.action=NULL, 
		int_meth=c("GL", "CAV_SIM", "SIM_3_8", "BOOLE", "GLM", "BANDS"),
		bands=NULL,
		stept=NULL,              
		init=NULL,
		optim.control=list(trace=100, REPORT=1, fnscale=-1, maxit=25), 
		optim_meth=c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "SANN", "Brent"),
		vartype =  c("oim", "opg", "none"),
		debug=FALSE, 
		...
){
	# choice of spline basis
	if (!missing(Spline))
		Spline <- match.arg(Spline)
	else {
		Spline <- "b-spline"
	}
	
	# type of spline in the correction model 
	if (!missing(Spline.table))
		Spline.table <- match.arg(Spline.table)
	else {
		Spline.table <- "restricted B-splines"
	}
	
	
	
	
# setting up variables
	call <- match.call()
	mf <- match.call(expand.dots = FALSE)
	
	
	
	m <- match(c("formula", "data", "rate", "weights", "logit_start", "logit_end", "na.action", "Id"),
			names(mf), 0L)
	if (m[1]==0)
		stop ("The formula argument is required.")
	
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- quote(stats::model.frame)
	special <- c("NPH","NLL", "NPHNLL", "nl", "td", "nltd", "WCEI") 
	Terms <- if (missing(data)){
				terms(formula, specials=special)
			} else {
				terms(formula, specials=special, data = data)
			}
	
	mf$formula <- Terms
	
	
	
	mf <- eval(mf,  sys.parent())
	mt <- attr(mf, "terms")
	
	# remove data containing NAs
	na.act <- attr(m, "na.action")
	if( !is.null( na.act )) {
		data <- data[na.act,]
	}
	
	
# detection of WCEI effets
	list_var_WCEI <- all_specials_vars(Terms, specials="WCEI",
			unique = FALSE,
			order="formula")
	is_wce_model <- length(list_var_WCEI) > 0
	
	
# switch to additive WCEI for the first WCEI
	if( baselinehazard == FALSE & firstWCEIadditive == TRUE & is_wce_model ){
		is_wce1add_model <- TRUE
		var_WCEI.Bh <-  list_var_WCEI[[1]]
	}
	else {
		is_wce1add_model <- FALSE
	}
	
#    if (is_wce_model) {
#      stop("Weighted cumulative exposure effect (WCEI() not yet implementd")
#    }
	
	
	# pour l'intercept, voir avec l'option utilisée dans gamma0(t)
	intercept <- attr(mt, "intercept")
	
	Y <- model.extract(mf, "response")
	if (!inherits(Y, "Surv")) {
		stop("Response must be a survival object")
	}
	Survtype <- attr(Y, "type")
	
	if ((ncol(Y) ==  2) ) {
		if (is_wce_model) {
			stop(gettextf("With weighted cumulative exposure effect (WCE()), you must specify a (start, end] survival object in the left hand side of formula", domaine=NA))
		} else {
			if (Survtype != "right"){
				stop(gettextf("flexrsurv does not support %s type of censoring with (0, end] survival data", dQuote(Survtype), domaine=NA))
			} else {
				Y <- cbind(rep(0, dim(Y)[1]), Y)
			}
		}
	} else if ((ncol(Y) ==  3) && (Survtype != "counting") ) {
		stop(gettextf("flexrsurv does not support %s type of censoring with (start, end] survival data", dQuote(Survtype), domaine=NA))
	}
	
	
# add finalT in Y
# get FirstId and Lastid
	Id <- model.extract(mf, "Id")
	if (is.null(Id)){
		Id <- 1:(dim(Y)[1])
	} else {
		data$Id <- Id
	}
	if(is.null(Id)){
		Id <- 1:(dim(Y)[1])
	}
	
	Index <- getIndexId(Id)
	FirstId <- Index$FirstId
	LastId <- Index$LastId
	Y <- cbind(Y, Y[LastId, 2])
	
	
	rate <- as.vector(model.extract(mf, "rate"))
	if(!is.null(rate) && !is.numeric(rate))
		stop("'rate' must be a numeric vector")
	if (is.null(rate)){
		stop("'rate' must be specified")
	} else {
		data$rate <- rate
	}
	
	
	
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
	}
	
	
	
	# build designs matrices
	
	if (!is_wce_model) {
		des <- ReadDesignFlexrsurv(Terms=Terms, modframe=mf, data=data, rate=rate, Spline_t0=Spline_t0, intercept.Bh=intercept.Bh)
	} else {
		des <- ReadDesignFlexrsurvWCEI(Terms=Terms, modframe=mf, data=data, rate=rate, Spline_t0=Spline_t0, intercept.Bh=intercept.Bh)
	}
	
	
	
	listdf_excess<- list(gamma0 = length(des$coef2param$gamma0), 
			alpha0 = length(des$coef2param$alpha0),
			beta0  = length(des$coef2param$beta0),
			alpha  = length(des$coef2param$alpha),
			beta   = length(des$coef2param$beta),
			eta0 = length(des$coef2param$eta0))
	df_excess <- sum(unlist(listdf_excess))
	
	
	
# correction model
	mfclt <- match.call(expand.dots = FALSE)
	
	mclt <- match(c("formula.table", "logit_start", "logit_end", "knots.table", "degree.table", "Spline.table", "Spline.CLT"),
			names(mfclt), 0L)
	
	if(sum(mclt[-c(2,3)]) > 0){
		Intercept_B <- TRUE
		
		# analysis with correction model
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
		
#design matrix of the proportional terms of the correction model 
		if (mclt[1]!=0 ){
			mclt2 <- match(c("formula.table", "data"),
					names(mfclt), 0L)
			mfclt <- mfclt[c(1L, mclt2)]
			mfclt$drop.unused.levels <- TRUE
			mfclt[[1L]] <- quote(stats::model.frame)
			names(mfclt)[[2L]] <- "formula"
			
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
		df_correction <- nbrass + nBX0
		
	} else {
		df_correction <- 0
		nbrass <- 0
		nBX0 <- 0
	}
	
	list_df <- c(listdf_excess, list(brass0=nbrass, balpha0=nBX0))
	
	return(list(ndf.excess = df_excess, ndf.correction = df_correction, list_df=list_df))
	
}
