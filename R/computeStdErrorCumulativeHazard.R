.computeStdErrorCumulativeHazard_GA0B0AB<-function(GA0B0AB, var,
                                                   Y, X0, X, Z, 
                                                   step, Nstep, 
                                                   intTD=intTD_NC, intweightsfunc=intweights_CAV_SIM,
                                                   intTD_base=intTD_base_NC,
                                                   nT0basis,
                                                   Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                                                   ialpha0, nX0,
                                                   ibeta0, nX,
                                                   ialpha, ibeta,
                                                   nTbasis,
                                                   Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                                                   listZSplineBasis,
                                                   Intercept_t_NPH=rep(TRUE, nX), 
                                                   debug=FALSE,  ...){
  # compute std error of a cumulative hazard (log rate) if the model
  # rate = exp(f(t)%*%gamma + X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
  # using the delta method
  # var(Cumhaz) = t(grad(Cumahaz)) var(param) grad(CumaHaz)
  #
  # code from gr_ll_flexrsurv_fromto_GA0B0AB or gr_ll_flexrsurv_GA0B0AB 
  #
  #################################################################################################################
  #################################################################################################################
  #  the coef of the first t-basis is constraint to 1 for nat-spline, and n-sum(other beta) if BS using expand() method
  #################################################################################################################
  #################################################################################################################
  #################################################################################################################
  # GA0B0AB ; vector of all coefs
                                        # gamma0 = GA0B0AB[1:nY0basis]
                                        # alpha0= GA0B0AB[ialpha0]
                                        # beta0= matrix(GA0B0AB[ibeta0], ncol=nX, nrow=nTbasis)
                                        # alpha= diag(GA0B0AB[ialpha])
                                        # beta= expand(matrix(GA0B0AB[ibeta], ncol=nZ, nrow=nTbasis-1))
                                        # beta does not contains coef for the first t-basis
  #################################################################################################################
  # Y : object of class Surv
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : objesct of class DeSignMatrixLPHNLL of time dépendent variables (spline basis expended)
  # step : lag of subinterval for numerical integration fr each observation
  # Nstep : number of lag for each observation
  # intTD : function to perform numerical integration 
  # intweightfunc : function to compute weightsfor numerical integration
  # Knots_t0=NULL,Intercept_t0=FALSE, degree_t0=3, Boundary.knots_t0 time spline parameters for baseline hazard
  # Knots_t=NULL,Intercept_t=FALSE, degree_t0=, Boundary.knots_t  time spline parameters for time-dependant effects (same basis for each TD variable)
  # nT0basis : number of spline basis for NPHLIN effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  # nTbasis : number of time spline basis
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z@signature
  # returned value : the log liikelihood of the model
  
if (debug) cat("# computing gradient of the cumulative hazard: .computeStdErrorCumulativeHazard_GA0B0AB\n")

if(dim(Y)[2] == 2){
gradcumhaz <-gr_cumhaz_flexrsurv_GA0B0AB(GA0B0AB=GA0B0AB, var=var,
                                             Y=Y, X0=X0, X=X, Z=Z, 
                                             step=step, Nstep=Nstep, 
                                             intTD=intTD, intweightsfunc=intweightsfunc,
                                             intTD_base=intTD_base,
                                             nT0basis=nT0basis,
                                             Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                                             ialpha0=ialpha0, nX0=nX0,
                                             ibeta0=ibeta0, nX=nX,
                                             ialpha=ialpha, ibeta=ibeta,                             
                                             nTbasis=nTbasis,
                                             Spline_t = Spline_t,
                                             Intercept_t_NPH=Intercept_t_NPH,
                                             debug=debug)
} else {
gradcumhaz <-gr_cumhaz_flexrsurv_fromto_GA0B0AB(GA0B0AB=GA0B0AB, var=var,
                                             Y=Y, X0=X0, X=X, Z=Z, 
                                             step=step, Nstep=Nstep, 
                                             intTD=intTD, intweightsfunc=intweightsfunc,
                                             intTD_base=intTD_base,
                                             nT0basis=nT0basis,
                                             Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                                             ialpha0=ialpha0, nX0=nX0,
                                             ibeta0=ibeta0, nX=nX,
                                             ialpha=ialpha, ibeta=ibeta,                             
                                             nTbasis=nTbasis,
                                             Spline_t = Spline_t,
                                             Intercept_t_NPH=Intercept_t_NPH,
                                             debug=debug)

}

    return(sqrt(apply(gradcumhaz * tcrossprod(gradcumhaz, var), 1, sum)))

  
}


.computeStdErrorCumulativeHazard_GA0B0AB_bh<-function(GA0B0AB, var,
                                                   Y, X0, X, Z, 
                                                   step, Nstep, 
                                                   intTD=intTD_NC, intweightsfunc=intweights_CAV_SIM,
                                                   intTD_base=intTD_base_NC,
                                                   nT0basis,
                                                   Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                                                   ialpha0, nX0,
                                                   ibeta0, nX,
                                                   ialpha, ibeta,
                                                   nTbasis,
                                                   Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                                                   listZSplineBasis,
                                                   Intercept_t_NPH=rep(TRUE, nX), 
                                                   debug=FALSE,  ...){
  # compute std error of a cumulative hazard (log rate) if the model
  # rate = f(t)%*%gamma exp( X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
  # using the delta method
  # var(Cumhaz) = t(grad(Cumahaz)) var(param) grad(CumaHaz)
  #
  # code from gr_ll_flexrsurv_fromto_GA0B0AB_bh or gr_ll_flexrsurv_GA0B0AB_bh 
  #
  #################################################################################################################
  #################################################################################################################
  #  the coef of the first t-basis is constraint to 1 for nat-spline, and n-sum(other beta) if BS using expand() method
  #################################################################################################################
  #################################################################################################################
  #################################################################################################################
  # GA0B0AB ; vector of all coefs
                                        # gamma0 = GA0B0AB[1:nY0basis]
                                        # alpha0= GA0B0AB[ialpha0]
                                        # beta0= matrix(GA0B0AB[ibeta0], ncol=nX, nrow=nTbasis)
                                        # alpha= diag(GA0B0AB[ialpha])
                                        # beta= expand(matrix(GA0B0AB[ibeta], ncol=nZ, nrow=nTbasis-1))
                                        # beta does not contains coef for the first t-basis
  #################################################################################################################
  # Y : object of class Surv
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : objesct of class DeSignMatrixLPHNLL of time dépendent variables (spline basis expended)
  # step : lag of subinterval for numerical integration fr each observation
  # Nstep : number of lag for each observation
  # intTD : function to perform numerical integration 
  # intweightfunc : function to compute weightsfor numerical integration
  # Knots_t0=NULL,Intercept_t0=FALSE, degree_t0=3, Boundary.knots_t0 time spline parameters for baseline hazard
  # Knots_t=NULL,Intercept_t=FALSE, degree_t0=, Boundary.knots_t  time spline parameters for time-dependant effects (same basis for each TD variable)
  # nT0basis : number of spline basis for NPHLIN effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  # nTbasis : number of time spline basis
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z@signature
  # returned value : the log liikelihood of the model
  
if (debug) cat("# computing gradient of the cumulative hazard: .computeStdErrorCumulativeHazard_GA0B0AB_bh\n")


if(dim(Y)[2] == 2){
gradcumhaz <-gr_cumhaz_flexrsurv_GA0B0AB_bh(GA0B0AB=GA0B0AB, var=var,
                                             Y=Y, X0=X0, X=X, Z=Z, 
                                             step=step, Nstep=Nstep, 
                                             intTD=intTD, intweightsfunc=intweightsfunc,
                                             intTD_base=intTD_base,
                                             nT0basis=nT0basis,
                                             Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                                             ialpha0=ialpha0, nX0=nX0,
                                             ibeta0=ibeta0, nX=nX,
                                             ialpha=ialpha, ibeta=ibeta,                             
                                             nTbasis=nTbasis,
                                             Spline_t = Spline_t,
                                             Intercept_t_NPH=Intercept_t_NPH,
                                             debug=debug)
} else {
  gradcumhaz <-gr_cumhaz_flexrsurv_fromto_GA0B0AB_bh(GA0B0AB=GA0B0AB, var=var,
                                             Y=Y, X0=X0, X=X, Z=Z, 
                                             step=step, Nstep=Nstep, 
                                             intTD=intTD, intweightsfunc=intweightsfunc,
                                             intTD_base=intTD_base,
                                             nT0basis=nT0basis,
                                             Spline_t0=Spline_t0, Intercept_t0=Intercept_t0,
                                             ialpha0=ialpha0, nX0=nX0,
                                             ibeta0=ibeta0, nX=nX,
                                             ialpha=ialpha, ibeta=ibeta,                             
                                             nTbasis=nTbasis,
                                             Spline_t = Spline_t,
                                             Intercept_t_NPH=Intercept_t_NPH,
                                             debug=debug)
}
    return(sqrt(apply(gradcumhaz * tcrossprod(gradcumhaz, var), 1, sum)))

  
}







