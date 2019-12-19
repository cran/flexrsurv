\name{logLik.flexrsurv}
\alias{logLik.flexrsurv}
\alias{nobs.flexrsurv}

\title{Log-Likelihood and the number of observations for a \code{flexrsuv} fit.}
\usage{
\method{logLik}{flexrsurv}(object, \dots)

\method{nobs}{flexrsurv}(object, \dots)
}
\arguments{
 \item{object}{any object of class \code{flexrsuv} results of a \code{\link{flexrsurv}} fit.}
 \item{\dots}{not used}
}
\description{
Function to extract Log-Likelihood and the number of observations from a \code{flexrsuv} or \code{flexrsuvclt} fit.
}


%\details{
%
% }
\value{
  \code{logLik} returns a standard \code{logLik} object (see \code{\link[stats:logLik]{logLik}})
  
  \code{nobs} returns a single number, normally an integer. 

}
\seealso{
  \code{\link[stats:logLik]{logLik}}, \code{\link[stats:nobs]{nobs}}.
}

%\author{
%  Michel Grzebyk
%}
