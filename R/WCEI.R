# weighted cumulative exposure index
# returned value : matrix of dim length(timevar) X nb of bases of Spline object (with/without intercept)
#                 of the contribution at timevar of exposure of subject with increment exposure X at Tfrom
# length(x) = length(timevar)
# length(FromT) = length(timevar) or 1
# spline basis can be specified in Spline.WCEI as an object of class .SplineBasis
#                               or by specifying the type, knots, degree, boundary
WCEI <-function(x,timevar,
		fromT=0,
		Spline.WCEI=NULL, 
		Spline=c("m-spline","b-spline","tp-spline", "tpi-spline"),
		Knots.t=NULL, 
		Degree.t=3, 
		Intercept.t=TRUE, 
		Boundary.knots.t = range(c(timevar, fromT)), 
		Keep.duplicates.t = TRUE,
		outer.ok=TRUE,
		...){
# x is the increment of the exposure
	
	
	if (is.null(Spline.WCEI)){
		Spline <- match.arg(Spline)
		if (Spline=="m-spline") {
			Spline.WCEI  <- MSplineBasis(knots=c(Boundary.knots.t[1], Knots.t, Boundary.knots.t[2]),
					degree=Degree.t,
					keep.duplicates=Keep.duplicates.t,
					log=FALSE)
		}
		else if (Spline=="b-spline") {
			Spline.WCEI  <- BSplineBasis(knots=c(Boundary.knots.t[1], Knots.t, Boundary.knots.t[2]),
					degree=Degree.t,
					keep.duplicates=Keep.duplicates.t,
					log=FALSE)
		}
		else if (Spline=="tp-spline") {
			Spline.WCEI  <- TPSplineBasis(knots=Knots.t,
					degree=Degree.t,  
					min=Boundary.knots.t[1],
					max=Boundary.knots.t[2],
					log=FALSE,
					type="standard")
		}
		else if (Spline=="tpi-spline") {
			Spline.WCEI  <- TPSplineBasis(knots=Knots.t,
					degree=Degree.t,  
					min=Boundary.knots.t[1],
					max=Boundary.knots.t[2],
					log=FALSE,
					type="standard")
		}
	}
	IS <- integrate(Spline.WCEI)
	evaluate(object=IS, x=timevar - fromT, intercept=Intercept.t, outer.ok=TRUE) * x
	
}


