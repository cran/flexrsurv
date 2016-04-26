\name{summary.flexrsurv}
\alias{summary.flexrsurv}
\alias{print.summary.flexrsurv}

\title{ Summarizing Flexible Relative Survival Model Fits}

\description{
   \code{summary} methods for class \code{flexrsurv}.
   Produces and prints summaries of the results of a fitted Relative Survival Model 
}

\usage{
\method{summary}{flexrsurv}(object, correlation = FALSE, symbolic.cor = FALSE, ...)

\method{print}{summary.flexrsurv}(x, digits = max(3L, getOption("digits") - 3L), 
      symbolic.cor = x$symbolic.cor,
      signif.stars = getOption("show.signif.stars"), ...)

}

\arguments{
  \item{object}{
	an object of class "flexrsurv", usually, a result of a call to \code{\link{flexrsurv}}.
      }
    \item{x}{an object of class \code{"summary.flexrsurv"}, usually, a result of a call to \code{summary.flexrsurv}.}    
  \item{correlation}{logical; if \code{TRUE}, the correlation matrix of
    the estimated parameters is returned and printed.}
  \item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in
    a symbolic form (see \code{\link{symnum}}) rather than as numbers.}
  \item{digits}{
	the number of significant digits to use when printing.
  }
  \item{signif.stars}{
	logical. If TRUE,'significance stars' are printed for each coefficient.
      }
  \item{...}{
    further arguments passed to or from other methods.
      }
}


\details{
   \code{print.summary.glm} tries to be smart about formatting the coefficients, standard errors, etc. 
   and additionally gives `significance stars' if \code{signif.stars} is \code{TRUE}. 

   Correlations are printed to two decimal places (or symbolically): to see the actual correlations 
   print \code{summary(object)$correlation} directly.

   The dispersion of a GLM is not used in the fitting process, but it is needed to find standard 
   errors. If dispersion is not supplied or NULL, the dispersion is taken as 1 for the binomial and
   Poisson families, and otherwise estimated by the residual Chisquared statistic (calculated from 
   cases with non-zero weights) divided by the residual degrees of freedom.

   

}
\value{
   The function summary.flexrsurv computes and returns a list of summary statistics of the fitted flexible relative survival model given in \code{object}.
   The returned value is an object of class "\code{summary.flexrsurv}", which a list with components:
\item{call}{the "\code{call}" component from \code{object}.}
\item{terms}{the "\code{terms}" component from \code{object}.}
\item{coefficients}{the matrix of coefficients, standard errors, z-values and p-values.}
\item{cov}{the estimated covariance matrix of the estimated coefficients.}
\item{correlation}{(only if \code{correlation} is true.) the estimated correlations of the estimated coefficients.}
\item{symbolic.cor}{(only if \code{correlation} is true.) the value of the argument symbolic.cor.}
\item{loglik}{the "\code{loglik}" component from object.}
\item{df.residual}{the "\code{df.residual}" component from object.}
}




\examples{
\dontrun{
# data from package relsurv
data(rdata, package="relsurv")
# rate table from package relsurv
data(slopop, package="relsurv")

# get the death rate from slopop for rdata
rdata$iage <- findInterval(rdata$age*365.24, attr(slopop, "cutpoints")[[1]])
rdata$iyear <- findInterval(rdata$year, attr(slopop, "cutpoints")[[2]])
therate <- rep(-1, dim(rdata)[1])
for( i in 1:dim(rdata)[1]){
  therate[i] <- slopop[rdata$iage[i], rdata$iyear[i], rdata$sex[i]]
}

rdata$slorate <- therate


# change sex coding
rdata$sex01 <- rdata$sex -1

# fit a relative survival model with a non linear effetc of age

fit <- flexrsurv(Surv(time,cens)~sex01+NLL(age, Knots=60, Degree=3), 
                 rate=slorate, data=rdata,
                 knots.Bh=1850,  # one interior knot at 5 years
                 degree.Bh=3,
                 Spline = "b-spline",
                 initbyglm=TRUE, 
                 initbands=seq(from=0, to=5400, by=200)
                 int_meth= "CAV_SIM",
                 step=50
                 )

summary(fit)
}
}


