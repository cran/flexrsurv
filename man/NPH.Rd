\name{NPH}
\alias{NPH}


\title{Non Proportional Hazard effect}

\description{
Generate the design matrix of spline basis for non proportional effect. 
}

\usage{
NPH(x, 
    timevar, 
    Spline = c("b-spline", "tp-spline", "tpi-spline"),
    Knots.t = NULL, 
    Degree.t = 3, 
    Intercept.t = TRUE, 
    Boundary.knots.t = range(timevar), 
    Keep.duplicates.t = TRUE, 
    outer.ok = TRUE, 
    ...)

}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the predictor variable.
}
  \item{timevar}{
the time variable.
}
  \item{Spline}{
a character string specifying the type of spline basis. "b-spline" for B-spline basis, "tp-spline" for truncated power basis and 
"tpi-spline" for monotone (increasing) truncated power basis.
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
range of time period which is analysed. By default it is the range of time variable.
}
  \item{Keep.duplicates.t}{
Should duplicate interior knots be kept or removed. Defaults is \code{FALSE}, which removes duplicate knots with a warning if duplicate interior knots are found.
}
  \item{outer.ok}{
logical indicating how are managed \code{timevar} values outside the knots. If \code{FALSE}, return \code{NA}, if \code{TRUE}, return \code{0} for the corresponding \code{timevar} values.
}
  \item{\dots}{
    not used
  }}

\details{
   \code{NPH} is based on package \code{\link[orthogonalsplinebasis:SplineBasis]{orthogonalsplinebasis}}
}

\seealso{
   \code{\link[=NLL]{NLL}}, and 
   \code{\link[=NPHNLL]{NPHNLL}}.
}

