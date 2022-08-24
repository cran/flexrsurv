\name{flexrsurvclt}
\Rdversion{0.2}
\alias{Flexrsurvclt}
\alias{flexrsurvclt}
\alias{flexrsurvclt.ll}


\title{ Fit Relative Survival Model and Correct Life Tables}

\description{
   \code{flexrsurvclt} is used to fit relative survival regression model. 
transition package. 
}

\usage{
flexrsurvclt(formula=formula(data),
	formula.table=NULL, 
	data=parent.frame(),
	Id,
	baselinehazard=TRUE,
	firstWCEIadditive=FALSE,
	knots.Bh,
	degree.Bh=3,
	intercept.Bh=TRUE,			
	Spline=c("b-spline", "tp-spline", "tpi-spline"), 
	log.Bh=FALSE,
	bhlink=c("log", "identity"),
	Min_T=0,
	Max_T=NULL,
	model=c("additive","multiplicative"),
	rate, 
	logit_start, 
	logit_end,
	logit_start_byperiod = NULL,
	logit_end_byperiod = NULL,
	weights_byperiod = NULL, 
	Id_byperiod = NULL,
	contrasts.table = NULL,
	knots.table=c(-2.5,0,2.5),   
	degree.table=3,
	Spline.table=c("restricted B-splines"), 
	Spline.CLT=R2bBSplineBasis(knots=c(-2.5,0,2.5), degree=3),
	model_correction = c("cohort", "period"),
	weights=NULL,
	na.action=NULL,
	datacontrol=NULL,
	Idcontrol,
	ratecontrol, 
	logit_startcontrol, 
	logit_endcontrol,
	logit_start_byperiodcontrol = NULL, 
	logit_end_byperiodcontrol = NULL, 
	weights_byperiodcontrol = NULL,
	Id_byperiodcontrol = NULL,	
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
	control.glm=list(epsilon=1e-8, maxit=100, trace=FALSE, 
					epsilon.glm=.1, maxit.glm=25),
	vartype =  c("oim", "opg", "none"),
	varmethod = c("optim", "numDeriv.hessian", "numDeriv.jacobian"),
	numDeriv.method.args=list(eps=5e-7, d=0.001, 
			zero.tol=sqrt(.Machine$double.eps/7e-4), r=4, v=2),
	debug=FALSE
   )

flexrsurvclt.ll(formula=formula(data),
	formula.table=NULL, 
	data=parent.frame(),
	Id,
	baselinehazard=TRUE,
	firstWCEIadditive=FALSE,
	knots.Bh,
	degree.Bh=3,
	Spline=c("b-spline", "tp-spline", "tpi-spline"), 
	log.Bh=FALSE,
	bhlink=c("log", "identity"),
	intercept.Bh=TRUE,
	Min_T=0,
	Max_T=NULL,
	model=c("additive","multiplicative"),
	rate, 
	logit_start, 
	logit_end,
	logit_start_byperiod = NULL,
	logit_end_byperiod = NULL,
	weights_byperiod = NULL, 
	Id_byperiod = NULL,
	contrasts.table = NULL,
	knots.table=c(-2.5,0,2.5),   
	degree.table=3,
	Spline.table=c("restricted B-splines"), 
	Spline.CLT=R2bBSplineBasis(knots=c(-2.5,0,2.5), degree=3),
	model_correction = c("cohort", "period"),
	weights=NULL,
	na.action=NULL,
	datacontrol=NULL,
	Idcontrol,
	ratecontrol, 
	logit_startcontrol, 
	logit_endcontrol,
	logit_start_byperiodcontrol = NULL, 
	logit_end_byperiodcontrol = NULL, 
	weights_byperiodcontrol = NULL,
	Id_byperiodcontrol = NULL,	
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
	vartype =  c("oim", "opg", "none"),
	varmethod = c("optim", "numDeriv.hessian", "numDeriv.jacobian"),
	numDeriv.method.args=list(eps=5e-7, d=0.001, 
	 		zero.tol=sqrt(.Machine$double.eps/7e-4), r=4, v=2),
	debug=FALSE
   )

}

\arguments{
  \item{formula}{
      a formula object, with the response on the left of a ~ operator, and the terms on the
      right. The response must be a survival object as returned by the \code{\link[survival:Surv]{Surv}} function. 
}
  \item{formula.table}{
    a formula object, with empty left hand side, and the terms on the
right. This is the formula of the proportional part of the correction
model for the table table 
     
}
  \item{data}{
      a data.frame in which to interpret the variables named in the formulas.
}
  \item{Id}{
      vector whose unique values defines the Ids of the subjects.
}
\item{baselinehazard}{
  if FALSE, no baseline hazard in the model
}
\item{firstWCEIadditive}{
  if TRUE, the first WCEI term in the formula is considered as the baseline
}
  \item{knots.Bh}{
      the internal breakpoints that define the spline used to estimate the baseline hazard. 
      Typical values are the mean or median for one knot, quantiles for more knots.
}
  \item{degree.Bh}{
      degree of the piecewise polynomial of the baseline hazard. Default is 3 for cubic splines.
}
  \item{intercept.Bh}{
      TRUE if the first bases is included in the baseline hazard. Default is TRUE.
}
  \item{Spline}{ 
      a character string specifying the type of spline basis. "b-spline" for B-spline basis, 
	  "tp-spline" for truncated power basis and "tpi-spline" for monotone (increasing) truncated power basis.
}
  \item{log.Bh}{
      logical value: if TRUE, an additional basis equal to log(time) is added to the spline bases of time.
}
  \item{bhlink}{
  character string specifying the link function of the baseline hazard: 
  Default is \code{bhlink="log"} for including the baseline in the exponential; 
  if \code{bhlink="identity"}, the baseline hazard is out of the exponential.
}
\item{Min_T}{
      minimum of time period which is analysed. Default is \code{max(0.0, min(bands) )}.
}
   \item{Max_T}{
      maximum of time period which is analysed. Default is \code{max(c(bands, timevar))}
	  }
  \item{model}{
      character string specifying the type of model for both non-proportional and non linear effects.    
The model \code{method=="additive"} assumes effects as explained in Remontet et al.(2007), the model \code{method=="multiplicative"} assumes effects as explained in Mahboubi et al. (2011).
}
   \item{rate}{
      a vector of the background rate for a relevant comparative population to be used in the fitting process. Should be a numeric vector (for relative survival model). 
     \code{rate} is evaluated in the same way as variables in \code{formula}, that is first 
     in \code{data} and then in the environment of \code{formula}.
}
   \item{logit_start}{
      a vector of the logit of the cumulative hazard at the start of the
interval in the life table.
     \code{logit_start} is evaluated in the same way as variables in \code{formula}, that is first 
     in \code{data} and then in the environment of \code{formula}.
}
   \item{logit_end}{
      a vector of the logit of the cumulative hazard at the end of the
interval in the life table.
     \code{logit_end} is evaluated in the same way as variables in \code{formula}, that is first 
     in \code{data} and then in the environment of \code{formula}.
}
   \item{logit_start_byperiod, logit_end_byperiod, weights_byperiod, Id_byperiod}{
A REMPLIR   
}
  \item{knots.table}{
      the internal breakpoints on the logit scale that define the knots
      of the spline used to estimate the correction model of the life table. 
}
  \item{degree.table}{
    degree of the piecewise polynomial of the spline used to estimate
    the correction model of the life table.
    Default is 3 for cubic splines.
}
\item{contrasts.table}{
    an optional list. See the \code{contrasts.arg} of  \code{\link[stats:model.matrix]{model.matrix()}}.
}
  \item{Spline.table}{ 
      a character string specifying the type of spline basis of the the
     correction model of the life table. In this version, only "restricted B-splines" is available. 
     "restricted B-splines" are B-spline basis with linear extrapolation + 2nd derivative at boundaries == 0.
}
  \item{Spline.CLT}{ 
        a S4 object with method deriv() and evaluate().
        The spline basis of the correction of the life table can be specified either by the parameters (\code{knots.table}, \code{degree.table})
        or an S4 object that ca be used for this purpose.
        IMPORTANT : the coef of the first basis is constraints to one and \code{evaluate(deriv(spline_B), left_boundary_knots) == 1}  
}
  \item{model_correction}{
      character string specifying A COMPLETER.    
      \code{method_correction="cohort"} when the provided logit are those of the survival of individuals;
      \code{method_correction="period"} when the provided logit are those of the survival fuction of age distribution by period.
}   \item{weights}{
      an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. If not null, the total likelihood is the weighted sum of individual likelihood. 
} 
  \item{na.action}{
      a missing-data filter function, applied to the model.frame. If \code{NULL}, default is \code{options()$na.action}.
}
  \item{datacontrol}{
      a data.frame in which to interpret the variables named in the formula for the control group.
}
  \item{Idcontrol, ratecontrol, logit_startcontrol, logit_endcontrol, weightscontrol}{
  Id, rate, logit of the cumulative hazard at the start and the end of
  the intervalle in the life table, and weights for the control group
}
   \item{logit_start_byperiodcontrol, logit_end_byperiodcontrol, weights_byperiodcontrol, Id_byperiodcontrol}{
A REMPLIR   
}
  \item{int_meth}{
      character string specifying the the numerical integration method. Possible values are 
    "GL" for Gauss-Legendre method, "CAV_SIM" for Cavalieri-Simpson's rule, "SIM_3_8" for the Simpson's 3/8 rule, 
	  "BOOLE" for the Boole's rule, or "BANDS" for the midpoint rule with specified bands.
}
   \item{bands}{
      bands used to split data in the numerical integration when \code{int_meth="BANDS"}).
}
   \item{npoints}{
      number of points used in the numerical integration when \code{int_meth="GL"}).
}
   \item{stept}{
      scalar value of the time-step in numerical integration. It is required only when  \code{int_meth="CAV_SIM"} or \code{"SIM_3_8"} or \code{"BOOLE"}. 
    If no value is supplied,  \code{Max_T/500} is used.
}
  \item{init}{
      starting values of the parameters.  
}
  \item{initbyglm}{
      a logical value indicating indicating how are found or refined init values. If TRUE, the fitting method described in Remontet et al.(2007) is ued to find or refine starting values. This may speedup the fit. If FALSE, the maximisation of the likelihood starts at values given in \code{init}. If \code{init=NULL}, the starting values correspond to a constant net hazard equal to the ratio of the number of event over the total number of person-time.
}
   \item{initbands}{
      bands used to split data when \code{initbyglm=TRUE}.
}
  \item{optim.control}{
      a list of control parameters passed to the \code{\link[stats:optim]{optim()}} function.
}
  \item{optim_meth}{
       method to be used to optimize the likelihood. 
      See \code{\link[stats:optim]{optim}}.
}
  \item{Coptim.control}{
      a list of control parameters passed to the \code{\link[stats:constrOptim]{constrOptim()}} function 
      See \code{\link[stats:constrOptim]{constrOptim}}.
}
  \item{lower, upper}{
      Bounds on the variables for the "L-BFGS-B" method, or bounds in which to search for method "Brent". 
      See \code{\link[stats:optim]{optim}}.
}
  \item{control.glm}{
      a list of control parameters passed to the \code{\link[stats:glm]{glm()}} function when \code{method="glm"}. 
}
  \item{vartype}{
      character string specifying the type of variance matrix computed by \code{flexrsurv}: the inverse of the hessian matrix computed at the MLE estimate (ie. the inverse of the observed information matrix) if \code{vartype="oim"}, the inverse of the outer product of the gradients if \code{vartype="opg"}. The variance is not computed when \code{vartype="none"}.
}
  \item{varmethod}{
      character string specifying the method to compute the hessian matrix when \code{vartype="oim"}. If \code{varmethod="oim"}, the hessian matrixe is computed by \code{\link[stats:optim]{optim}}. If \code{varmethod="numDeriv.hessian"}, the hessian matrix is computed by \code{\link[numDeriv:hessian]{numDeriv:hessian}} with \code{method="Richardson"}. If \code{varmethod="numDeriv.jacobian"}, the hessian matrixe is computed by \code{\link[numDeriv:jacobian]{numDeriv:jacobian}} with \code{method="Richardson"}. 
}
  \item{numDeriv.method.args}{
      arguments passed to \code{\link[numDeriv:hessian]{numDeriv:hessian}} or \code{\link[numDeriv:jacobian]{numDeriv:jacobian}} when \code{varmethod="numDeriv.hessian"} or \code{varmethod="numDeriv.jacobian"}. Arguments not specified remain with their default values as specified in details. See \code{\link[numDeriv:grad]{numDeriv:grad}} for details about these parameters.
}
  \item{debug}{
    control the volum of intermediate output
}
}

\details{
   A full description of the additive and the multiplicative both non-linear and non-proportional models is given respectively in Remontet (2007) and Mahboubi (2011).
   
   \code{flexrsurv.ll} is the workhorse function: it is not normally called
  directly.
}

\value{
   \code{flexrsurv} returns an object of class \code{"flexrsurv"}.
   An object of class \code{"flexrsurv"} is a list containing at least the following components: 

	\item{coefficients}{
          a named vector of coefficients
	}
	\item{loglik}{
          the log-likelihood
	}
	\item{var}{
		estimated covariance matrix for the estimated coefficients
	}
	\item{informationMatrix}{
		estimated information matrix
	}
	\item{bhlink}{
          the linkk of baseline hazard:
	  if \code{"identity"} baseline = sum g0_i b_i(t);
	  if \code{"log"} log(baseline) = sum g0_i b_i(t); 
	}
	\item{init}{
          vector of the starting values supplied
	}
	\item{converged}{
          logical, Was the optimlizer algorithm judged to have converged? 
	}
	\item{linear.predictors}{
		the linear fit on link scale (not including the baseline hazard term if \code{bhlink = "identity"})
	}
	\item{fitted.values}{
		the estimated value of the hazard rate at each event time, obtained by transforming the linear predictors by the inverse of the link function
	}
	\item{cumulative.hazard}{
		the estimated value of the cumulative hazard in the time interval
	}
	\item{call}{
		the matched call
	}
	\item{formula}{
		the formula supplied
	}
	\item{terms}{
		the \code{\link{terms}} object used
	}
	\item{data}{
		the \code{data} argument
		}
	\item{rate}{
		the rate vector used
	}
	\item{time}{
		the time vector used
	}
	\item{workingformula}{
		the formula used by the fitter
	}
	\item{optim.control}{
		the value of the \code{optim.control} argument supplied
	}
	\item{control.glm}{
		the value of the \code{control.glm} argument supplied
	}   
    \item{method}{
	the name of the fitter function used
	}
}

\references{
    Mahboubi, A., M. Abrahamowicz, et al. (2011). "Flexible modeling of the effects of continuous prognostic factors in relative survival." Stat Med 30(12): 1351-1365. \doi{10.1002/sim.4208}

    Remontet, L., N. Bossard, et al. (2007). "An overall strategy based on regression models to estimate relative survival and model the effects of prognostic factors in cancer survival studies." Stat Med 26(10): 2214-2228. \doi{10.1002/sim.2656}

}


\seealso{
   \code{\link{print.flexrsurv}}, 
   \code{\link{summary.flexrsurv}}, 
   \code{\link{logLik.flexrsurv}}, 
   \code{\link{predict.flexrsurv}}, 
%%   \code{\link{plot.flexrsurv}}, 
   \code{\link{NPH}}, 
   \code{\link{NLL}}, and 
   \code{\link{NPHNLL}}.
}

                   

\examples{
\donttest{

if (requireNamespace("relsurv", quietly = TRUE) & requireNamespace("date", quietly = TRUE)) {

	library(date)
	# data from package relsurv
	data(rdata, package="relsurv")
	
	class(rdata$year)<-"integer"
	
	
	# rate table from package relsurv
	data(slopop, package="relsurv")
	
	# get the death rate at event (or end of followup) from slopop for rdata
	rdata$iage <- findInterval(rdata$age*365.24+rdata$time, attr(slopop, "cutpoints")[[1]])
	rdata$iyear <- findInterval(rdata$year+rdata$time, attr(slopop, "cutpoints")[[2]])
	therate <- rep(-1, dim(rdata)[1])
	for( i in 1:dim(rdata)[1]){
	  therate[i] <- slopop[rdata$iage[i], rdata$iyear[i], rdata$sex[i]]
	}
	
	rdata$slorate <- therate
	
	# get the logit_start and logit_end
	# logit start at age 18
	
	
	tmpsurv <- Surv(rep(0, length(rdata$time)), rdata$time, rdata$cens)
	
	
	HH <- getHazardFromTable(tmpsurv, startdate=rdata$year,
	       startage=rdata$age*365.25 , matchdata=rdata, ratetable=slopop,
	       age="age", year="year",
	       rmap=list(sex=sex),
	       agemin=18,
	       ratename = "poprate", cumrateendname ="cumrateend", cumrateentername ="cumrateenter"
	      ) 
	
	rdata$slorate <- HH$poprate
	rdata$logit_start <- log(exp(HH$cumrateenter)-1)
	rdata$logit_end <- log(exp(HH$cumrateend)-1)
	
	rdata$Id <- 1:dim(rdata)[1]
	
	
	
	
	# change sex coding
	rdata$sex01 <- rdata$sex -1
	
	# fit a relative survival model with a non linear effect of age 
	#   without correction of life table
	#   partial likelihood
	fit00 <- flexrsurvclt(Surv(time,cens)~sex01+NLL(age, Knots=60, Degree=3,
	                                           Boundary.knots = c(24, 95)), 
	                      rate=slorate, 
			      data=rdata,
			      knots.Bh=1850,  # one interior knot at 5 years
	                 degree.Bh=3,
	                 Max_T=5400,
	                 Spline = "b-spline",
	                 initbyglm=TRUE,
	                 initbands=seq(0, 5400, 100), 
	                 int_meth= "BANDS",
	                 bands=seq(0, 5400, 50)
	                 )
	summary(fit00)
	
	# fit a relative survival model with a non linear effect of age 
	#   without correction of life table
	#   full likelihood
	fit0 <- flexrsurvclt(Surv(time,cens)~sex01+NLL(age, Knots=60, Degree=3,
	                                           Boundary.knots = c(24, 95)), 
	                                           rate=slorate, 
					           logit_start=logit_start,
	                    logit_end=logit_end,
				data=rdata,
				Id=Id,
	                 knots.Bh=1850,  # one interior knot at 5 years
	                 degree.Bh=3,
	                 Max_T=5400,
	                 Spline = "b-spline",
	                 initbyglm=TRUE,
	                 initbands=seq(0, 5400, 100), 
	                 int_meth= "BANDS",
	                 bands=seq(0, 5400, 50)
	                 )
	summary(fit0)
	
	# fit a relative survival model with a non linear effect of age
	#   with correction of life table
	#   full likelihood
	fit1 <- flexrsurvclt(Surv(time,cens)~sex01+NLL(age, Knots=60, Degree=3,
	                                           Boundary.knots = c(24, 95)), 
	                 rate=slorate, 
					                     logit_start=logit_start,
	                    logit_end=logit_end,
				data=rdata,
				Id=Id,
	                 knots.Bh=1850,  # one interior knot at 5 years
	                 degree.Bh=3,
	                 Max_T=5400,
	                 Spline = "b-spline",
	             Spline.CLT=flexrsurv:::R2bBSplineBasis(knots=c(-2.5,0,2.5), degree=3),
	                 initbyglm=TRUE,
	                 initbands=seq(0, 5400, 100), 
	                 int_meth= "BANDS",
	                 bands=seq(0, 5400, 50)
	                 )
	summary(fit1)
	
	print(coef(fit1))
	
	# fit a relative survival model with a non linear effect of age
	#   with correction of life table, strabified by sex
	#   full likelihood
	fit2 <- flexrsurvclt(Surv(time,cens)~sex01+NLL(age, Knots=60, Degree=3,
	                                           Boundary.knots = c(24, 95)), 
		    formula.table= ~sex,
	                 rate=slorate, 
					                     logit_start=logit_start,
	                    logit_end=logit_end,
				data=rdata,
				Id=Id,
	                 knots.Bh=1850,  # one interior knot at 5 years
	                 degree.Bh=3,
	                 Max_T=5400,
	                 Spline = "b-spline",
	             Spline.CLT=flexrsurv:::R2bBSplineBasis(knots=c(-2.5,0,2.5), degree=3),
	                 initbyglm=TRUE,
	                 initbands=seq(0, 5400, 100), 
	                 int_meth= "BANDS",
	                 bands=seq(0, 5400, 50)
	                 )
	summary(fit2)
	
	AIC(fit0, fit1, fit2)
}
}
}

%
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
\keyword{ survival }
\keyword{ nonlinear }
\keyword{ models }


