\name{WCEI}
\alias{WCEI}


\title{Weighted cumulative exposure index}

\description{
Generate the spline basis matrix for Weighted cumulative exposure index. 
}

\usage{
WCEI(x, 
    timevar, 
    fromT=0,
    Spline.WCEI=NULL, 
    Spline = c("m-spline", "b-spline", "tp-spline", "tpi-spline"), 
    Knots.t = NULL, 
    Degree.t = 3, 
    Intercept.t = TRUE,
    Boundary.knots.t = range(c(timevar, fromT)),
    Keep.duplicates.t = TRUE, 
    outer.ok = TRUE, 
    ...)


}

\arguments{
  \item{x}{
the exposure variable.
}
  \item{timevar}{
the time variable.
}
  \item{fromT}{
Time at which starts exposure
}
  \item{Spline.WCEI}{
        a S4 object with method deriv(), evaluate() and predict().
}
  \item{Spline}{
a character string specifying the type of spline basis. "b-spline" for B-spline basis, 
"tp-spline" for truncated power basis and "tpi-spline" for monotone (increasing) truncated power basis.
}
  \item{Knots.t}{
the internal breakpoints that define the spline used to estimate the WCEI effect. By default there 
are none.
}
  \item{Degree.t}{
degree of splines which are considered.
}
  \item{Intercept.t}{
a logical value indicating whether intercept/first basis of spline should be considered.
}
  \item{Boundary.knots.t}{
range of variable which is analysed.
}
  \item{Keep.duplicates.t}{
  Should duplicate interior knots be kept or removed. Defaults is \code{FALSE}, which removes
  duplicate knots with a warning if duplicate interior knots are found.
}
  \item{outer.ok}{
  logical indicating how are managed \code{x} values outside the knots. If \code{FALSE},
  return \code{NA}, if \code{TRUE}, return \code{0} for the corresponding \code{x} values.
}
  \item{\dots}{
    not used
  }
}

\details{
   \code{WCEI} is based on package \code{\link[orthogonalsplinebasis:SplineBasis]{orthogonalsplinebasis}}
}
\seealso{
   \code{\link[=NLL]{NLL}} \code{\link[=NPH]{NPH}} and 
   \code{\link[=NPHNLL]{NPHNLL}}.
}

