\name{flexrsurv}
\alias{flexrsurv}
\alias{flexrsurv.ll}

\title{Fit Relative Survival Model}

\description{
   \code{flexrsurv} is used to fit relative survival regression model. 
Time dependent variables, non-proportionnal (time dependent) effects,
non-linear effects are implemented using Splines (B-spline and truncated power basis).
Simultaneously non linear and non proportional effects are implemented
using approaches developed by Remontet et al.(2007) and Mahboubi et al. (2011). 
}

\usage{
flexrsurv(formula=formula(data),
   data=parent.frame(), 
   knots.Bh,
   degree.Bh=3,
   Spline=c("b-spline", "tp-spline", "tpi-spline"), 
   log.Bh=FALSE,
   bhlink=c("log", "identity"),
   Min_T=0,
   Max_T=NULL,
   model=c("additive","multiplicative"),
   rate=NULL, 
   weights=NULL,
   na.action=NULL,
   int_meth=c("GL", "CAV_SIM", "SIM_3_8", "BOOLE", "BANDS"),
   npoints=20,   
   stept=NULL,              
   bands=NULL,
   init=NULL,
   initbyglm=TRUE,
   initbands=bands,
   optim.control=list(trace=100, REPORT=1, fnscale=-1, maxit=25), 
   optim_meth=c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "SANN", "Brent"),
   control.glm=list(epsilon=1e-8, maxit=100, trace=FALSE, epsilon.glm=1e-1, maxit.glm=25),
   vartype =  c("oim", "opg", "none"),
   debug=FALSE
   )


flexrsurv.ll(formula=formula(data), 
   data=parent.frame(), 
   knots.Bh=NULL,   
   degree.Bh=3,
   Spline=c("b-spline", "tp-spline", "tpi-spline"), 
   log.Bh=FALSE,
   bhlink=c("log", "identity"),
   Min_T=0,
   Max_T=NULL,
   model=c("additive","multiplicative"),
   rate=NULL, 
   weights=NULL,
   na.action=NULL, 
   int_meth=c("GL", "CAV_SIM", "SIM_3_8", "BOOLE", "GLM", "BANDS"),
   npoints=20,   
   stept=NULL,
   bands=NULL,
   init=NULL,
   optim.control=list(trace=100, REPORT=1, fnscale=-1, maxit=25), 
   optim_meth=c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "SANN", "Brent"),
   vartype =  c("oim", "opg", "none"),
   debug=FALSE
   )
}

\arguments{
  \item{formula}{
      a formula object, with the response on the left of a ~ operator, and the terms on the
      right. The response must be a survival object as returned by the \code{\link[survival:Surv]{Surv}} function. 
}
  \item{data}{
      a data.frame in which to interpret the variables named in the formula.
}
  \item{knots.Bh}{
      the internal breakpoints that define the spline used to estimate the baseline hazard. 
      Typical values are the mean or median for one knot, quantiles for more knots.
}
  \item{degree.Bh}{
      degree of the piecewise polynomial of the baseline hazard. Default is 3 for cubic splines.
}
  \item{Spline}{ 
      a character string specifying the type of spline basis. "b-spline" for B-spline basis, 
	  "tp-spline" for truncated power basis and "tpi-spline" for monotone (increasing) truncated power basis.
}
  \item{log.Bh}{
      logical value: if TRUE, an additional basis equal to log(time) is added to the spline bases of time.
}
  \item{bhlink}{
  logical value: if TRUE, log of baseline hazard is modelled, if FALSE, the baseline hazard is out of the log.
}
\item{Min_T}{
      minimum of time period which is analysed. Default is \code{max(0.0, min(bands) )}.
}
   \item{Max_T}{
      maximum of time period which is analysed. Default is \code{max(c(bands, timevar))}
	  }
  \item{model}{
      character string specifying the type of model for both non-proportionnal and non linear effects. 
      The model \code{method=="additive"} assumes effects as explained in Remontet et al.(2007), 
      the model \code{method=="multiplicative"} assumes effects as explained in Mahboubi et al. (2011).
}
   \item{rate}{
      an optional vector of the background rate for a relevant comparative population to be used in the fitting process. 
      Should be a numeric vector (for relative survival model). 
     \code{rate} is evaluated in the same way as variables in \code{formula}, that is first 
     in \code{data} and then in the environment of \code{formula}.
}
   \item{weights}{
      an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. 
      If not null, the total likelihood is the weighted sum of individual likelihood. 
} 
  \item{na.action}{
      a missing-data filter function, applied to the model.frame, after any subset argument has been used. 
      Default is options()$na.action.
}
  \item{int_meth}{
      character string specifying the the numerical integration method. Possible values are 
	  "GL" for Gauss-Legendre quadrature, 
      "CAV_SIM" for Cavalieri-Simpson's rule, "SIM_3_8" for the Simpson's 3/8 rule, 
	  "BOOLE" for the Boole's rule, or "BANDS" for the midpoint rule with specified bands.
}
   \item{npoints}{
      number of points used in the Gauss-Legendre quadrature (when \code{int_meth="GL"}).
}
   \item{stept}{
      scalar value of the time-step in numerical integration. It is required only when  \code{int_meth="CAV_SIM"} or \code{"SIM_3_8"} or \code{"BOOLE"}. 
    If no value is supplied,  \code{Max_T/500} is used.
}
   \item{bands}{
      bands used to split data in the numerical integration when \code{int_meth="BANDS"}.
}
  \item{init}{
      starting values of the parameters.  
}
  \item{initbyglm}{
      a logical value indicating indicating how are found or refined init values. If TRUE, the fitting method described in Remontet 
	  et al.(2007) is ued to find or refine starting values. This may speedup the fit. If FALSE, the maximisation of the likelihood 
	  starts at values given in \code{init}. If \code{init=NULL}, the starting values correspond to a constant net hazard equal to 
	  the ratio of the number of event over the total number of person-time.
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
 \item{control.glm}{
      a list of control parameters passed to the \code{\link[stats:glm]{glm()}} function when \code{method="glm"}. 
}
  \item{vartype}{
      character string specifying the type of variance matrix computed by \code{flexrsurv}: the inverse of the hessian matrix computed 
	  at the MLE estimate (ie. the inverse of the observed information matrix) if \code{vartype="oim"}, the inverse of the outer 
	  product of the gradients if \code{vartype="opg"}. The variance is not computed when \code{vartype="none"}.
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

if (requireNamespace("relsurv", quietly = TRUE)) {

	# data from package relsurv
	data(rdata, package="relsurv")
	
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
	
	# change sex coding
	rdata$sex01 <- rdata$sex -1
	
	# fit a relative survival model with a non linear effect of age
	fit <- flexrsurv(Surv(time,cens)~sex01+NLL(age, Knots=60, Degree=3,
	                                           Boundary.knots = c(24, 95)), 
	                 rate=slorate, data=rdata,
	                 knots.Bh=1850,  # one interior knot at 5 years
	                 degree.Bh=3,
	                 Max_T=5400,
	                 Spline = "b-spline",
	                 initbyglm=TRUE,
	                 initbands=seq(0, 5400, 100), 
	                 int_meth= "BANDS",
	                 bands=seq(0, 5400, 50)
	                 )
	summary(fit)
	
	# fit a relative survival model with a non linear & non proportional effect of age
	fit2 <- flexrsurv(Surv(time,cens)~sex01+NPHNLL(age, time, Knots=60,
	                                               Degree=3,
	                                               Knots.t = 1850, Degree.t = 3), 
	                 rate=slorate, data=rdata,
	                 knots.Bh=1850,  # one interior knot at 5 years
	                 degree.Bh=3,
	                 Spline = "b-spline",
	                 initbyglm=TRUE, 
	                 int_meth= "BOOLE",
	                 step=50
	                 )
	summary(fit2, correlation=TRUE)
	
}
}
}

%
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
\keyword{ survival }
\keyword{ nonlinear }
\keyword{ models }


