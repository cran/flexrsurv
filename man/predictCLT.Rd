\name{predictCLT}
\alias{predictCLT}
\alias{predictCLT.default}
\alias{predictCLT.flexrsurvclt}

\title{Predictions for relational life table model}

\description{
  Predict the relational model for a life table correction model fitted by \code{flexrsuvclt} or for specified knots, degree and coefficients }

\usage{
  predictCLT (...)

  \method{predictCLT}{flexrsurvclt}(object, newdata= NULL,
   type = c("clt", "correction"),
   se.fit=FALSE, na.action=na.pass, newcoef = NULL, ...)
   
   \method{predictCLT}{default}(knots, degree, newdata, newcoef, ...)
}

\arguments{
\item{object}{
the results of a flexrsurvclt fit. 
}
\item{newdata}{
  Optional new vector of logarithm of the cumulative distribution odds (LCDO) at which to do predictions.
  If absent predictions are for values of the LCDO used in the original fit (\code{logit_end} parameter in the call to \code{flexrsuvclt})). 
}
\item{newcoef}{
  Optional new coefficients for which to do predictions.
  If absent predictions are for the coefficients of the fitted model in \code{object}. 
}
\item{type}{
the type of predicted value. 
Choices are \code{"clt"} or \code{"correction"} to compute the corrected logarithm of the cumulative distribution odds. 
}
\item{se.fit}{
  if TRUE, pointwise standard errors are produced for the predictions
  (not yet implemented).
}
\item{knots, degree}{
  knots and degree of the relational model.
}
\item{na.action}{
  function determining what should be done with missing values in
  \code{newdata}. The default is to predict \code{NA}.
}
\item{\dots}{For future methods}

}

\value{
a vector or a list containing the predicted relational model (element \code{"fit"}) and their
standard errors (element \code{"se.fit"}) if the se.fit option is TRUE. 
}

\details{
  \code{predictCLT} with \code{knots} and \code{degree} arguments computes corrected values of  . 
}



  \seealso{
    \code{\link{predict.flexrsurvclt}}, \code{\link{flexrsurv}}, \code{\link{flexrsurvclt}}  }

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
	
	
	
	
	corrected_logit_end <- predictCLT(fit1)
	
	
	
	try_logit_end <- predictCLT(knots=c(-2.5,0,2.5), degree=3, newcoef = c(0.5, 2), 
		newdata = rdata$logit_end  )
	
	plot(rdata$logit_end, corrected_logit_end)
	points(rdata$logit_end, try_logit_end, col = 2)
	

}
}
}