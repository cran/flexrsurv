\name{print.flexrsurv}
\alias{print.flexrsurv}

\title{ Print a Short Summary of a Relative Survival Model}

\description{
  Print number of observations, number of events, the formula, the
estimated coeficients and the log likelihood. 
}
\usage{
  \method{print}{flexrsurv}(x, 
      digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{
     the result of a call to the \code{flexrsuv} function.
}
  \item{digits}{
     the minimum number of significant digits to be printed in values, see \code{\link{print.default}}.
   }
     \item{...}{
    other options 
  }

}


\seealso{
   The default method \code{\link{print.default}}, and help for the
function \code{\link{flexrsurv}}.

}

