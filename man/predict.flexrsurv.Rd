\name{predict.flexrsurv}
\alias{predict.flexrsurv}
\alias{predict.flexrsurvclt}

\title{Predictions for a relative survival model}

\description{
  Predict linear predictors, hazard and cumulative
hazard for a model fitted by \code{flexrsuv} }

\usage{
  \method{predict}{flexrsurv}(object, newdata= NULL,
   type = c("lp", "link", "risk", "hazard", "hazardrate",
             "rate", "loghazard", "log", "lograte",
             "cumulative.rate", "cumulative.hazard", "cumulative", "cum",
             "survival", "surv", "netsurv"),
   se.fit=FALSE, na.action=na.pass, ...)
}

\arguments{
\item{object}{
the results of a flexrsurv fit. 
}
\item{newdata}{
  Optional new data at which to do predictions.
  If absent predictions are for the data frame used in the original fit. 
}
\item{type}{
the type of predicted value. 
Choices are the linear predictor (\code{"lp"}, \code{"log"},
\code{"loghazard"}, \code{"lograte"}), the hazard (\code{"rate"}, \code{"hazard"}, \code{"hazardrate"}, \code{"risk"}) or the cumulative
hazard (\code{"cum"}, \code{"cumulative.hazard"}, \code{"cumulative"}). 
}
\item{se.fit}{
  if TRUE, pointwise standard errors are produced for the predictions
  (not available for cumulative hazard).
}
\item{na.action}{
  function determining what should be done with missing values in
  \code{newdata}. The default is to predict \code{NA}.
}
\item{\dots}{For future methods}

}

\value{
a vector or a list containing the predictions (element \code{"fit"}) and their
standard errors (element \code{"se.fit"}) if the se.fit option is TRUE. 
}

\details{
  For cumulative hazard, the cumulative hazard is computed from 0 until
the given end time. The cumulative hazard is computed using the same numerical
integration method as the one used to fit the model. 
}

\note{
To work correctly, arguments \code{Boundary.knots} and
\code{Boundary.knots.t} must be included in the call to \code{NPH()}, \code{NLL()} and
\code{NPHNLL()} in the formula of \code{flexrsurv} 
  }

  \seealso{
    \code{\link{predict}}, \code{\link{flexrsurv}}, \code{\link{flexrsurvclt}}  }

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
	# centering age
	rdata$agec <- rdata$age- 60
	
	# fit a relative survival model with a non linear effect of age
	fit <- flexrsurv(Surv(time,cens)~sex01+NLL(age, Knots=60, Degree=3,
	                                           Boundary.knots = c(24, 95)), 
	                 rate=slorate, data=rdata,
	                 knots.Bh=1850,  # one interior knot at 5 years
	                 degree.Bh=3,
	                 Spline = "b-spline",
	                 initbyglm=TRUE, 
	                 int_meth= "BOOLE",
	                 step=50
	                 )
	summary(fit, correlation=TRUE)
	
	
	
	newrdata <- rdata
	newrdata$age <- rep(60, length(rdata$age))
	newrdata$sex <- factor(newrdata$sex, labels=c("m", "f"))
	
	
	
	linpred <- predict(fit, newdata=newrdata, type="lp", se.fit=TRUE )
	predhazard <- predict(fit, newdata=newrdata, type="hazard" , se.fit=TRUE )
	predcumhazard <- predict(fit, newdata=newrdata, type="cum", se.fit=TRUE)
	
	
	require(ggplot2)
	tmp <- cbind(newrdata, linpred)
	glp <- ggplot(tmp, aes(time, colour=sex))
	glp + geom_ribbon(aes(ymin = fit-2*se.fit, ymax = fit + 2*se.fit, fill=sex)) +
	   geom_line(aes(y=fit)) +
	   scale_fill_manual(values = alpha(c("blue", "red"), .3))
	
	
	tmp <- cbind(newrdata, predhazard)
	glp <- ggplot(tmp, aes(time, colour=sex))
	glp + geom_ribbon(aes(ymin = fit-2*se.fit, ymax = fit + 2*se.fit, fill=sex)) +
	   geom_line(aes(y=fit))  +
	   scale_fill_manual(values = alpha(c("blue", "red"), .3))
	
	
	tmp <- cbind(newrdata, predcumhazard)
	glp <- ggplot(tmp, aes(time, colour=sex))
	glp + geom_ribbon(aes(ymin = fit-2*se.fit, ymax = fit + 2*se.fit, fill=sex)) +
	   geom_line(aes(y=fit)) +
	   scale_fill_manual(values = alpha(c("blue", "red"), .3))
}
}
}