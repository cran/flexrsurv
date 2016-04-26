\name{NLL}
\alias{NLL}


\title{Non Log-Linear effect}

\description{
Generate the spline basis matrix for non log-linear effect. 
}

\usage{
NLL(x, 
    Spline = c("b-spline", "tp-spline", "tpi-spline"), 
    Knots = NULL, 
    Degree = 3, 
    Intercept = FALSE,
    Boundary.knots = range(x), 
    Keep.duplicates = TRUE, 
    outer.ok = TRUE, 
    ...)


}

\arguments{
  \item{x}{
the predictor variable.
}
  \item{Spline}{
a character string specifying the type of spline basis. "b-spline" for B-spline basis, 
"tp-spline" for truncated power basis and "tpi-spline" for monotone (increasing) truncated power basis.
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
  \item{outer.ok}{
  logical indicating how are managed \code{x} values outside the knots. If \code{FALSE},
  return \code{NA}, if \code{TRUE}, return \code{0} for the corresponding \code{x} values.
}
  \item{\dots}{
    not used
  }
}

\details{
   \code{NLL} is based on package \code{\link[orthogonalsplinebasis:SplineBasis]{orthogonalsplinebasis}}
}
\seealso{
   \code{\link[=NPH]{NPH}} and 
   \code{\link[=NPHNLL]{NPHNLL}}.
}

