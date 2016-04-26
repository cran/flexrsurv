\name{NPHNLL}
\alias{NPHNLL}

\title{Non Proportional Hazard and Non Log-Linear effect}

\description{
Generate the design matrix of spline basis for both non log-linear and non proportional effect. 
}

\usage{
NPHNLL(x, 
       timevar, 
       model = c("additive", "multiplicative"), 
       Spline = c("b-spline", "tp-spline", "tpi-spline"), 
       Knots = NULL, 
       Degree = 3, 
       Intercept = FALSE, 
       Boundary.knots = range(x), 
       Knots.t = NULL, 
       Degree.t = 3,
       Intercept.t = (model == "multiplicative"), 
       Boundary.knots.t = range(timevar), 
       outer.ok = TRUE, 
       Keep.duplicates = TRUE, 
       xdimnames = ":XxXxXXxXxX ", 
       tdimnames = ":TtTtTTtTtT ")
}

\arguments{
  \item{x}{
the predictor variable.
}
  \item{timevar}{
the time variable.
}
  \item{model}{
character string specifying the type of model for both non-proportionnal and non linear effects.    
The model \code{method=="additive"} assumes effects as explained in Remontet et al.(2007), the model \code{method=="multiplicative"} assumes effects as explained in Mahboubi et al. (2011).
}
  \item{Spline}{
a character string specifying the type of spline basis. "b-spline" for B-spline basis, 
"tp-spline" for truncated power basis and "tpi-spline" for monotone (increasing) truncated power basis.
}
  \item{Knots}{
the internal breakpoints that define the spline used to estimate the NLL part of effect. By default 
there are none.
}
  \item{Degree}{
degree of splines of variable which are considered.
}
  \item{Intercept}{
a logical value indicating whether intercept/first basis of spline should be considered.
}
  \item{Boundary.knots}{
range of variable which is analysed.
}
  \item{Knots.t}{
the internal breakpoints that define the spline used to estimate the NPH part of effect. By default 
there are none.
}
  \item{Degree.t}{
degree of splines of time variable which are considered.
}
  \item{Intercept.t}{
a logical value indicating whether intercept/first basis of spline should be considered.
}
  \item{Boundary.knots.t}{
range of time period which is analysed.
}
  \item{Keep.duplicates}{
  Should duplicate interior knots be kept or removed. Defaults is \code{FALSE}, which removes
  duplicate knots with a warning if duplicate interior knots are found.
}
  \item{outer.ok}{
  logical indicating how are managed \code{timevar} or \code{x} values outside the knots. If \code{FALSE},
  return \code{NA}, if \code{TRUE}, return \code{0} for the corresponding \code{timevar} or \code{x} values.
}
  \item{xdimnames}{
  string to build dimnames of \code{x} bases
}
  \item{tdimnames}{
  string to build dimnames of \code{timevar} bases
}
}
\details{
   \code{NPHNLL} is based on package \code{\link[orthogonalsplinebasis:SplineBasis]{orthogonalsplinebasis}}
}

\references{
    Mahboubi, A., M. Abrahamowicz, et al. (2011). "Flexible modeling of the effects of continuous prognostic factors in relative survival." Stat Med 30(12): 1351-1365. \doi{DOI:10.1002/sim.4208}

    Remontet, L., N. Bossard, et al. (2007). "An overall strategy based on regression models to estimate relative survival and model the effects of prognostic factors in cancer survival studies." Stat Med 26(10): 2214-2228. \doi{10.1002/sim.2656}

}


\seealso{
   \code{\link[=NPH]{NPH}} and 
   \code{\link[=NLL]{NLL}}.
}
