\name{NLLbeta}
\alias{NLLbeta}
\alias{NPHalpha}

\title{Non Log-Linear effect and non proportional effect}

\description{
Internal functions not intended for users. 
}

\usage{

NLLbeta(y, x, 
    Spline = c("b-spline", "tp-spline", "tpi-spline"), 
    Knots = NULL, 
    Degree = 3, 
    Intercept = FALSE, 
    Boundary.knots = range(x), 
    Keep.duplicates = TRUE, 
    outer.ok = TRUE, 
    ...)

NPHalpha(x, 
    timevar, 
    Spline = c("b-spline", "tp-spline", "tpi-spline"),
    Knots.t = NULL, 
    Degree.t = 3, 
    Intercept.t = TRUE, 
    Boundary.knots.t = c(0, max(timevar)), 
    Keep.duplicates.t = TRUE, 
    outer.ok = TRUE, 
    ...)
}

\arguments{
  \item{x}{
the predictor variable.
}
  \item{timevar}{
the time variable.
}
\item{y}{
the name of variable for which tests NLL effect.
}
  \item{Spline}{
type of spline basis. "b-spline" for B-spline basis, "tp-spline" for truncated power basis 
and "tpi-spline" for monotone (increasing) truncated power basis.
}
  \item{Knots}{
the internal breakpoints that define the spline used to estimate the NLL effect. By default there 
are none.
}
  \item{Degree}{
degree of splines which are considered.
}
  \item{Intercept}{
a logical value indicating whether intercept/first basis of spline should be considered.
}
  \item{Boundary.knots}{
range of variable which is analysed.
}
  \item{Keep.duplicates}{
  Should duplicate interior knots be kept or removed. Defaults is \code{FALSE}, which removes
  duplicate knots with a warning if duplicate interior knots are found.
}
  \item{Knots.t}{
the internal breakpoints that define the spline used to estimate the NPH effect. By default 
there are none.
}
  \item{Degree.t}{
degree of splines which are considered.
}
  \item{Intercept.t}{
a logical value indicating whether intercept/first basis of spline should be considered.
}
  \item{Boundary.knots.t}{
range of time period which is analysed. By default it is \code{c(0, max(timevar))}.
}
  \item{Keep.duplicates.t}{
Should duplicate interior knots be kept or removed. Defaults is \code{FALSE}, which removes duplicate knots with a warning if duplicate interior knots are found.
}
  \item{outer.ok}{
 logical indicating how are managed \code{timevar} or \code{x} values outside the knots. If \code{FALSE},
  return \code{NA}, if \code{TRUE}, return \code{0} for the corresponding \code{timevar} or \code{x} values.
}
  \item{\dots}{
    not used
  }
}

\details{
Internal functions.
}
\value{
  \code{NLLbeta(x, y, ...)} returns \code{y * NLL(x, ...)}. 
  
  \code{NPH(x, timevar, ...)} is equal to \code{x * NPHalpha(x, timevar, ...)}. 
}

\seealso{
   \code{\link[=NPH]{NPH}}, 
   \code{\link[=NLL]{NLL}}, and 
   \code{\link[=NPHNLL]{NPHNLL}}.
}

\keyword{internal}
