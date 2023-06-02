\name{flexrsurv-package}
\Rdversion{2.0.11}
\alias{flexrsurv-package}
\alias{flexrsurvpackage}
\alias{Flexrsurv-package}
\alias{Flexrsurvpackage}
\docType{package}
\title{
Package for flexible relative survival analyses
}
\description{
\code{flexrsurv} is a package for parametric relative survival analyses. The package implements non-linear, 
non-proportional effects and both non proportional and non linear effects, using splines (B-spline and truncated power basis), 
Weighted Cumulative Index of Exposure effect, with correction model for the life table. 
Both non proportional and non linear effects are described in Remontet et al. (2007) \doi{10.1002/sim.2656} and 
Mahboubi et al. (2011) \doi{10.1002/sim.4208}. 
}
\details{
\tabular{ll}{
Package: \tab flexrsurv\cr
Type: \tab Package\cr
Version: \tab 2.0.11\cr
Date: \tab 2022-08-23\cr
License: \tab GPL(>)2.0)\cr
}

The main function is \code{flexrsurv()}

}
\author{
Michel Grzebyk and Isabelle Clerc-Urmès, with contributions from the CENSUR working survival group. \cr

Maintainer: <michel.grzebyk@inrs.fr>
}
\references{
    Mahboubi, A., M. Abrahamowicz, et al. (2011). "Flexible modeling of
the effects of continuous prognostic factors in relative survival." Stat
Med 30(12): 1351-1365. \doi{10.1002/sim.4208} \cr

    Remontet, L., N. Bossard, et al. (2007). "An overall strategy based
on regression models to estimate relative survival and model the effects
of prognostic factors in cancer survival studies." Stat Med 26(10):
2214-2228. \doi{10.1002/sim.2656} \cr
}
\keyword{package}
\keyword{model}
\keyword{nonlinear}
\keyword{survival}

\seealso{
   \code{\link[=flexrsurv]{flexrsurv}}
   }
