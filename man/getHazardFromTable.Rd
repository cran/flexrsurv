\name{getHazardFromTable}
\alias{getHazardFromTable}
\title{ computes expected hazards with respect to a reference life table
}
\description{
returns the cumulative hazard and the hazard rate of subjects in a reference life table
}
\usage{
getHazardFromTable(Y, startdate, startage, matchdata = NULL, 
        ratetable = survival::survexp.us, 
		age = 1, year = 2, rmap, agemin = 16, scale = 365.25, 
		ratename = "rateend", 
		cumrateendname = "cumrateend", 
		cumrateentername = "cumrateenter", 
		origin = "01/01/1970", format = "\%d/\%m/\%Y",
		left.open = FALSE, verbose=FALSE)
}
\arguments{
  \item{Y}{
An object with interval data. It can be an object of class \code{\link[survival:Surv]{Surv}} with arguments \code{time}, \code{time2} and \code{event} 
or a two-column matrix with starting time in the first column and ending time in the second column.
}
  \item{startdate}{
a numeric vector such that \code{as.Date(startadate)} is interpreted as the date of the start (when Y[,]==0).
}
  \item{startage}{
a numeric vector of age in days the start (when Y[,]==0).
}
  \item{matchdata}{
an optional data.frame in which to interpret the additional variables to be mapped to the \code{ratetable} variables.
}
  \item{ratetable}{
an object of class \code{\link[survival:ratetable]{ratetable}}, ie a table of event rates.
}
  \item{age, year}{
character values of the names of the age and period variables in the rate table.
}
  \item{rmap}{
	an optional list that maps data set names to the ratetable names. See \code{\link[survival:survexp]{survexp}} an example bellow.
}
  \item{agemin}{
numeric value of the age at which the cumulative hazard starts.
}
  \item{scale}{
numeric value to scale \code{agemin}. 
}
  \item{ratename, cumrateendname, cumrateentername}{
names of the returned variables
}
  \item{origin, format}{
arguments passed \code{\link[base:as.Date]{as.Date}}
}
  \item{left.open}{
logical, passed to \code{\link[base:findInterval]{findInterval}}
}
  \item{verbose}{
logical, if \code{true} the progression of the computation is output.
}

}
\details{
The cumulative rates are computed using \code{\link[survival:survexp]{survexp}}.
}
\value{
A \code{data.frame} with 3 columns with the rate at the ending time, the cumulative rate from agemin up to the starting time and upt to the ending time.
}

\seealso{
\code{\link{getPseudoHazardFromTable}} for the cumulative hazard in each period of a reference life table. 
\code{\link{getBrassHazardFromTable}} for the cumulative hazard in a corrected reference life table. 
\code{\link[survival:survexp]{survexp}}
}

\examples{
\donttest{

if (requireNamespace("relsurv", quietly = TRUE) & requireNamespace("date", quietly = TRUE)) {

	library(date)
	# data from package relsurv
	data(rdata, package="relsurv")
	
	# rate table from package relsurv
	data(slopop, package="relsurv")
	
	tmpsurv <- Surv(rep(0, length(rdata$time)), rdata$time, rdata$cens)
	
	HH <- getHazardFromTable(tmpsurv, startdate=rdata$year,
	               startage=rdata$age*365.24 , matchdata=rdata, ratetable=slopop,
	               age="age", year="year",
	               rmap=list(sex=sex),
	               agemin=18, scale=365.24,
	               ratename = "poprate", 
	               cumrateendname ="cumrateend", 
	               cumrateentername ="cumrateenter"
	               ) 
	summary(HH)
}
}						
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }