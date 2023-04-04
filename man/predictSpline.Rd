\name{predictSpline}
\alias{predictSpline}
\alias{predictSpline.default}
\alias{predictSpline.character}
\alias{predictSpline,character,numeric-method}

\title{Generic method for prediction of spline function}

\description{
  Predict a spline function by specifying its type, knots, degree and coefficients }

\usage{
  predictSpline (object, x, \dots)
   
   \method{predictSpline}{default}(object=c("b-spline", "tp-spline"), 
		x, knots, degree, keep.duplicates = FALSE, coef=1, ... )
}

\arguments{
\item{object}{
the type of spline to be predicted ("b-spline", the default, or "tp-spline") 
}
\item{x}{
  Vector of values at wich to predict the spline function. 
}
\item{knots, degree}{
  knots and degree of the relational model.
}
\item{keep.duplicates}{
  Should duplicate interior knots be kept or removed. Defaults is \code{FALSE}, which removes
  duplicate knots with a warning if duplicate interior knots are found.
}
\item{coef}{
  vector of coefficient of the spline function. 
}

\item{\dots}{not used}

}

\value{
	A vector the evaluated spline function with same length as x. 
}

\details{
  \code{predictSpline}  . 
}



  \seealso{
    \code{\link{predict.flexrsurvclt}}, \code{\link{flexrsurv}}, \code{\link{flexrsurvclt}}  }

\examples{

\donttest{
	
	predspline <- predictSpline("b-spline", 
		x= seq(from=-3, to = 3, by=.1), 
		coef = .5 * 1:5, 
		knots=c(-3,0,3), degree=3)
	plot(seq(from=-3, to = 3, by=.1), predspline)

}
}
