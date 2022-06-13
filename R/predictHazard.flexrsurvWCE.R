predictHazard <- function(object, ...) UseMethod("predictHazard")

predictHazard.flexrsurvWCE <- function(object, newdata = NULL,
		type=c("lp", "link", "terms", "risk", "hazard", "hazardrate", "rate", "loghazard", "log", "lograte"),
		se.fit = FALSE, 
		na.action = na.pass, ...){
	
	
	call <- match.call()
	type <-match.arg(type)
	
	type <- switch(type,
			risk = "risk",
			rate = "risk",
			hazard = "risk",
			hazardrate = "risk",
			lp = "link",
			log = "link",
			loghazard = "link",
			lograte = "link",
			link = "link",
			terms = "terms")  
	
	if (type == "terms") {
		stop("not yet implemented")
	}
	
	if ((missing(newdata) || is.null(newdata)) && !se.fit &&  type != "terms" ){
		pred <- switch(type, link = object$linear.predictors, 
				risk = object$fitted.values)
		if (!is.null(na.action)) {
			pred <- napredict(na.action, pred)
		}
	} else {
		if (missing(newdata) || is.null(newdata)) {
			newdata <- object$data
			if (!is.null(na.action)) {
				newdata <- na.action(newdata)
			}
		} else {
			newdata <- as.data.frame(newdata)
		}
		
		Terms <- object$terms
		
		# only the end times are needed, add beginning and stuts vars if needed
		survvars <- all.vars(lhs(formula(Terms)))
		if(length(survvars)==3){
# original survival object (tfrom , tto)
# if tfrom missing in newdata, add new columns with tfrom=0
			if(!(survvars[1] %in% dimnames(newdata)[[2]])){
				newdata[[survvars[1]]] <- rep(0, dim(newdata)[1])
			}
		}
		if(!(survvars[length(survvars)] %in% dimnames(newdata)[[2]])){
			newdata[[survvars[length(survvars)]]] <- rep(1, dim(newdata)[1])
		}
		
		m <- model.frame(Terms, data=newdata)
		
		
		Y <- model.extract(m, "response")
		
		Id <- as.vector(model.extract(m, "Id"))
		
		if (is.null(Id)){
			Id <- 1:(dim(Y)[1])
		} else {
			newdata$Id <- Id
		}
		if(is.null(Id)){
			Id <- 1:(dim(Y)[1])
		}
		nobs <-  length(unique(Id)) 
		
		if (is_wce_model) {
# add finalT in Y
# get FirstId and Lastid
			Index <- getIndexId(Id)
			FirstId <- Index$FirstId
			LastId <- Index$LastId
			Y <- cbind(Y, Y[LastId, 2])
		} 
		
		# get T-splines
		Max_T <- eval(as.expression(object$call$Max_T))
		Min_T <- eval(as.expression(object$call$Min_T))
		Spline <- eval(as.expression(object$call$Spline))
		if(is.null(Max_T)){
			Max_T <- max(Y[,1:(ncol(Y)-1)])
		}
		if(is.null(Min_T)){
			Min_T <- 0
		}
		
		list_var_WCEI <- all_specials_vars(Terms, specials="WCEI",
				unique = FALSE,
				order="formula")
		is_wce_model <- length(list_var_WCEI) > 0
		
		object$baselinehazard <- FALSE
		firstWCEIadditive <- ifelse(is.null(object$firstWCEIadditive), FALSE, object$firstWCEIadditive)
		
# switch to additive WCEI for the first WCEI
		if( firstWCEIadditive == TRUE & is_wce_model ){
			is_wce1add_model <- TRUE
			var_WCEI.Bh <-  list_var_WCEI[[1]]
		}
		else {
			is_wce1add_model <- FALSE
		}
		
# first set baseline or WCEadd	
		tmp_nparam <- 0
		
		baselinehazard <- object$baselinehazard
		if(baselinehazard == TRUE){
			bhlink <- object$bhlink
			knots.Bh <- eval(as.expression(object$call$knots.Bh)) 
			degree.Bh <- eval(as.expression(object$call$degree.Bh)) 
			intercept.Bh <- eval(as.expression(object$call$intercept.Bh)) 
			if(!is.null(intercept.Bh)){
						Intercept_t0 <- intercept.Bh				
					} else {
						Intercept_t0 <- TRUE				
					}
			if(Spline=="b-spline"){
				Spline_t0 <- BSplineBasis(knots=c(Min_T, knots.Bh, Max_T),
						degree=degree.Bh,
						keep.duplicates=TRUE)
				nT0basis <- getNBases(Spline_t0) - 1 +  Intercept_t0
				ngamma0 <- nT0basis
				
				Spline_t<- BSplineBasis(knots=c(Min_T, knots.Bh, Max_T),
						degree=degree.Bh,
						keep.duplicates=TRUE)
			} else if(Spline=="tp-spline") {
				Spline_t0 <-TPSplineBasis(knots=knots.Bh,
						degree=degree.Bh,
						min=Min_T,
						max=Max_T)
				nT0basis <- getNBases(Spline_t0) - 1 +  Intercept_t0
				ngamma0 <- nT0basis
				Spline_t <-TPSplineBasis(knots=knots.Bh,
						degree=degree.Bh,
						min=Min_T,
						max=Max_T)
			}
		}
		else {
			Spline_t0 <- NULL
			Intercept_t0 <- TRUE
			nT0basis <- 0
			ngamma0 <- 0
			knots.Bh <- NULL 
			degree.Bh <- NULL 
			bhlink <- "log"
		}
		tmp_nparam <- tmp_nparam + ngamma0	
		
		
		
		# extract components
		des <- object$des
		if (!is_wce_model) {
			newdes <- ReadDesignFlexrsurv(Terms=Terms, modframe=m, data=newdata, rate=NULL, Spline_t0=Spline_t0, intercept.Bh=Intercept_t0 )
		} else {
			newdes <- ReadDesignFlexrsurvWCEI(Terms=Terms, modframe=m, data=newdata, rate=NULL, Spline_t0=Spline_t0, intercept.Bh=Intercept_t0)
		}
		
		
		# by component
		
		if(is_wce1add_model == TRUE){
			nW <- 1L
			
			nWbasis <- getNBases(des$ISpline_W) - 1 + des$Intercept_W
			neta0 <- sum(nWbasis)
			ieta0 <- 1:neta0 
			Ieta0 <- 1:neta0 
			First.eta0 <- 1
			first.eta0 <- 1
			iWend <- cumsum(nWbasis)
			iWbeg <- c(1, iWend[-nW]-1)
		}
		else  {
			
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
		tmp_nparam <- 	tmp_nparam +  neta0
		
		# test if all time spline are identical
		if(!is.null(des$X) | !is.null(des$Z)){
			allSpline_T <- c(des$Spline_XT, des$Spline_ZT)
			if(length(allSpline_T) > 1){
				for(sb in allSpline_T[-1]){
					if(!identical(sb, allSpline_T[[1]])){
						stop("flexrsurv cannot handle different spline basis for time dependent effects")
					}
				}
			}
			Spline_t <- allSpline_T[[1]] 
		}
		
		
		
		# X0 linear and non linear effects
		X0<-des$X0
		if(!is.null(X0)){
			is.PH <- TRUE
			if(is.matrix(X0)) {
				nX0 <- dim(X0)[2]
			} else if(is.vector(X0)) {
				nX0 <- 1L
			} else {
				stop("error flexrsurv_LL.fit(): wrong type of X0") 
			}
			nalpha0<-nX0
			ialpha0<-1:nX0 + tmp_nparam
			Ialpha0<-1:nX0 
			First.alpha0<-tmp_nparam+1
			first.alpha0<-1
			X0<-newdes$X0
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
		tmp_nparam <- 	tmp_nparam +  nalpha0
		
		
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
			nTbasis <- getNBases(Spline_t)
			nTbasis_NPH <- getNBases(Spline_t) - 1 + Intercept_t_NPH
			if(is.matrix(X)) {
				nX <- dim(X)[2]
			} else if(is.vector(X)) {
				nX <- 1L
			}
			
			nbeta0 <- sum(nTbasis_NPH)
			ibeta0 <- 1:nbeta0 + tmp_nparam
			Ibeta0 <- 1:nbeta0 
			First.beta0 <- tmp_nparam + 1
			first.beta0 <- 1
			X<-as.matrix(newdes$X)
		} else  {
			X <- NULL
			Intercept_t_NPH <- NULL
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
		tmp_nparam <- 	tmp_nparam +  nbeta0
		
		
		# Z (NPHNLL effects
		# get splines for each NPHNLL variable
		
		if( !is.null(des$Z) ){
			Spline_Z <- des$Spline_Z
			Z<-DesignMatrixNPHNLL(Z=des$Z, listsplinebasis=Spline_Z, timesplinebasis=Spline_t)
			if( getNvar(Z)>0){
				is.NPHNLIN <- TRUE
				nTbasis <- getNBases(Spline_t)
				nZ<-getNvar(Z)
				nalpha <- getNparam(Z)
				ialpha <- 1:nalpha + tmp_nparam  
				Ialpha <- 1:nalpha  + nalpha0   
				First.alpha <- tmp_nparam + 1  
				first.alpha <- nalpha0  + 1  
				
				tmp_nparam <- 	tmp_nparam +  nalpha      
				
				# as first beta is constraints, nTbasis -1 beta per Z variable
				nbeta <- nZ * (nTbasis-1)
				ibeta <- 1:nbeta + tmp_nparam
				Ibeta <- 1:nbeta  + nbeta0 
				First.beta <- tmp_nparam + 1
				first.beta <- 1 + nbeta0
				
				tmp_nparam <- 	tmp_nparam +  nbeta
				Z<-DesignMatrixNPHNLL(Z=newdes$Z, listsplinebasis=Spline_Z, timesplinebasis=Spline_t)
				
			}
		} else  {
			is.NPHNLIN <- FALSE
#    nTbasis <- 0 already done with beta0
#      Z <- new("DesignMatrixNPHNLL")
			Z <- NULL
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
		
		
		# W, WCEI weightefd cumulative exposure index
		if(is_wce1add_model == FALSE & !is.null(des$W)){
			W<-as.matrix(des$W)
			nW <- dim(W)[2]
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
				nWbasis <- rep(0, nW)
				for(iW in 1:nW){
					nWbasis[iW] <- getNBases(ISpline_W[[iW]]) - 1 + Intercept_W[iW]
				}
				neta0 <- sum(nWbasis)
				ieta0 <- 1:neta0 + tmp_nparam
				Ieta0 <- 1:neta0 
				First.eta0 <- tmp_nparam + 1
				first.eta0 <- 1
				iWend <- cumsum(nWbasis)
				iWbeg <- c(1, iWend[-nW]-1)
				
			}
			W<-as.matrix(newdes$W)
		} else {
			W <- NULL
			ISpline_W <- NULL
			Intercept_W <- NULL
			nW <- 0
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
		
		
		
# get coefficients in order baseline hazard | linear & non-linear effects | non Prop efefct | nonlin_nonprop effects 
		# coefficient <- param[des$param2coef]
		# param[des$param2coef] <- coefficients
		alltheparameters <- rep(0, length(object$coefficients))
		alltheparameters[des$param2coef] <- object$coefficients
		# computes linear predictors
		
		if ( type != "terms" ){
			if(!is_wce_model) {
				linpred <- .computeLinearPredictor_GA0B0AB(allparam=alltheparameters,
						Y=Y, X0=X0, X=X, Z=Z,
						nT0basis=nT0basis,
						Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
						ialpha0=ialpha0, nX0=nX0,
						ibeta0= ibeta0, nX=nX, 
						ialpha=ialpha, 
						ibeta= ibeta, 
						nTbasis=nTbasis,
						Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
						bhlink=bhlink,
						debug=10000)    
				if (se.fit) {
					stderr <- .computeStdErrorLinearPredictor_GA0B0AB(allparam=alltheparameters,
							var=object$var,
							Y=Y, X0=X0, X=X, Z=Z,
							nT0basis=nT0basis,
							Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
							ialpha0=ialpha0, nX0=nX0,
							ibeta0= ibeta0, nX=nX, 
							ialpha=ialpha, 
							ibeta= ibeta,
							listZSplineBasis = des$Spline_Z,
							nTbasis=nTbasis,
							Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
							bhlink=bhlink,
							debug=FALSE)
				}
			}
			else {
				# is_wce_model=TRUE
				if(!is_wce1add_model){
					linpred <- .computeLinearPredictor_GA0B0ABE0(allparam=alltheparameters,
							Y=Y, X0=X0, X=X, Z=Z, W=W,
							Id=Id, FirstId=FirstId, LastId=LastId, 
							
							nT0basis=nT0basis,
							Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
							ialpha0=ialpha0, nX0=nX0,
							ibeta0= ibeta0, nX=nX, 
							ialpha=ialpha, 
							ibeta= ibeta, 
							nTbasis=nTbasis,
							Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
							ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW, 
							ISpline_W =ISpline_W,
							Intercept_W=Intercept_W,
							bhlink=bhlink,
							debug=10000)    
					if (se.fit) {
						stderr <- .computeStdErrorLinearPredictor_GA0B0ABE0(allparam=alltheparameters,
								var=object$var,
								Y=Y, X0=X0, X=X, Z=Z, W=W,
								Id=Id, FirstId=FirstId, LastId=LastId, 
								
								nT0basis=nT0basis,
								Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
								ialpha0=ialpha0, nX0=nX0,
								ibeta0= ibeta0, nX=nX, 
								ialpha=ialpha, 
								ibeta= ibeta,
								listZSplineBasis = des$Spline_Z,
								nTbasis=nTbasis,
								Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
								ieta0=ieta0, iWbeg=iWbeg, iWend=iWend, nW=nW, 
								ISpline_W =ISpline_W,
								Intercept_W=Intercept_W,
								bhlink=bhlink,
								debug=FALSE)
						
					}
				}
				else {
					# is_wce1add_model = TRUE
					linpred <- .computeLinearPredictor_fromto_1wceadd(allparam=alltheparameters,
							Y=Y, X0=X0, X=X, Z=Z, W=W,
							Id=Id, FirstId=FirstId, LastId=LastId, 
							
							ialpha0=ialpha0, nX0=nX0,
							ibeta0= ibeta0, nX=nX, 
							ialpha=ialpha, 
							ibeta= ibeta, 
							nTbasis=nTbasis,
							Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
							ieta0=ieta0, 
							ISpline_W =ISpline_W,
							Intercept_W=Intercept_W,
							wcelink="identity",
							debug=10000)    
					if (se.fit) {
						# not yet implemnted
						se.fit <- FALSE
#						stderr <- .computeStdErrorLinearPredictor_GA0B0ABE0(allparam=alltheparameters,
#								var=object$var,
#								Y=Y, X0=X0, X=X, Z=Z, W=W,
#								Id=Id, FirstId=FirstId, LastId=LastId, 
#								
#								ialpha0=ialpha0, nX0=nX0,
#								ibeta0= ibeta0, nX=nX, 
#								ialpha=ialpha, 
#								ibeta= ibeta,
#								listZSplineBasis = des$Spline_Z,
#								nTbasis=nTbasis,
#								Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
#								ieta0=ieta0,  
#								ISpline_W =ISpline_W,
#								Intercept_W=Intercept_W,
#								bhlink=bhlink,
#								debug=FALSE)
#						
						
					}
				}
			}
		}
		else {
#      linearPredictors <- .computeTermsLinearPredictor_GA0B0AB(allparam=alltheparameters,
#                                                                            var=object$var,
#                                                               Y=Y, X0=X0, X=X, Z=Z,
#                                                               nT0basis=nT0basis,
#                                                               Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
#                                                               ialpha0=ialpha0, nX0=nX0,
#                                                               ibeta0= ibeta0, nX=nX, 
#                                                               ialpha=ialpha, 
#                                                               ibeta= ibeta, 
#                                                               nTbasis=nTbasis,
#                                                               Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
#                                                               terms = terms,
#                                                               debug=FALSE)
#      
#      if (!se.fit) {
#        stdErrorLinearPredictors <- .computeStdErrorTermsLinearPredictor_GA0B0AB(allparam=alltheparameters,
#                                                                                 Y=Y, X0=X0, X=X, Z=Z,
#                                                                                 nT0basis=nT0basis,
#                                                                                 Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
#                                                                                 ialpha0=ialpha0, nX0=nX0,
#                                                                                 ibeta0= ibeta0, nX=nX, 
#                                                                                 ialpha=ialpha, 
#                                                                                 ibeta= ibeta, 
#                                                                                 nTbasis=nTbasis,
#                                                                                 Spline_t = Spline_t, Intercept_t_NPH=Intercept_t_NPH,
#                                                                                 terms = terms,
#                                                                                 debug=FALSE)
#      }
		} 
		
		if (type=="risk"){
			if(!is_wce1add_model){
				if(bhlink == "log"){
					linearPredictors <- pmax(exp(linpred ), .Machine$double.eps)
					if (se.fit) {
						stdErrorLinearPredictors <- stderr * linearPredictors
					}
				} else {
					linearPredictors <- pmax(linpred[,1] * exp(linpred[,2] ), .Machine$double.eps)
					if (se.fit) {
						varerr <- attr(stderr, "varerr")
						stdErrorLinearPredictors <- linearPredictors * (varerr[,1]/linpred[,1]^2 + 2*varerr[,2]/linpred[,1] + varerr[,3])
					}
				}
			}
			else {
				linearPredictors <- pmax(linpred[,1] * exp(linpred[,2] ), .Machine$double.eps)
				if (se.fit) {
					varerr <- attr(stderr, "varerr")
					stdErrorLinearPredictors <- linearPredictors * (varerr[,1]/linpred[,1]^2 + 2*varerr[,2]/linpred[,1] + varerr[,3])
				}
				
			}
		} else {
			linearPredictors <- linpred
			if (se.fit) {
				stdErrorLinearPredictors <- stderr
			}
		}
		
		if (se.fit) {
			# structure similar to predic.lm() 
			pred <- list(fit=linearPredictors , se.fit=stdErrorLinearPredictors)
		} else {
			pred <- linearPredictors
		}
		
	}
	
	pred
	
	
}




