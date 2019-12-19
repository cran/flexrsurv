gr_link_flexrsurv_GA0B0AB<-function(GA0B0AB, Y, X0, X, Z, 
                                    nT0basis,
                                    Spline_t0=BSplineBasis(knots=NULL, degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                                    ialpha0, nX0,
                                    ibeta0, nX,
                                    ialpha, ibeta,                             
                                    nTbasis,
                                    Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                                    Intercept_t_NPH=rep(TRUE, nX),
                                    debug=FALSE,  ...){
  # compute gradient of the linear parts of the rate (on the link scale) of the relatice survival model
  # rate = invlink(f(t)%*%gamma) exp(  X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
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
  # Y : object of class Surv (with ncol=2 or more)
  #                the time at which the predictors are computed is Y[,1] if ncol=2, Y[,2] if ncol>2
  #
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : objesct of class DeSignMatrixLPHNLL of time dépendent variables (spline basis expended)
  # Knots_t0=NULL,Intercept_t0=FALSE, degree_t0=3, Boundary.knots_t0 time spline parameters for baseline hazard
  # Knots_t=NULL,Intercept_t=FALSE, degree_t0=, Boundary.knots_t  time spline parameters for time-dependant effects (same basis for each TD variable)
  # nT0basis : number of spline basis for NPHLIN effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  # nTbasis : number of time spline basis
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z@signature
  # returned value : the log liikelihood of the model
  
  if (debug) cat("# computing gradient: gr_linkrate_flexrsurv_GA0B0AB\n")

  colEndTime <- ifelse(ncol(Y)==2, 1, 2)

  if(is.null(Z)){
    nZ <- 0
  } else {
    nZ <- Z@nZ
    }

  if(Intercept_t0){
    tmpgamma0 <- GA0B0AB[1:nT0basis]
  }
  else {
    tmpgamma0 <- c(0, GA0B0AB[1:nT0basis])
  }

  # contribution of time d?pendant effect
  # parenthesis are important for efficiency
  if(nZ) {
    # add a row for the first basis
    tBeta <- t(ExpandAllCoefBasis(GA0B0AB[ibeta], ncol=nZ,  value=1))
    # Zalpha est la matrice des alpha(Z)
    # parenthesis important for speed ?
    Zalpha <- Z@DM %*%( diag(GA0B0AB[ialpha]) %*% Z@signature )
  }  
  YT0 <- fevaluate(Spline_t0, Y[,colEndTime], intercept=Intercept_t0)
  YT <- fevaluate(Spline_t, Y[,colEndTime], intercept=TRUE)
  indx_without_intercept <- 2:getNBases(Spline_t)
  
#####################################################################"
# now computes the mean score


# d<dgamma0
    dLdgamma0 <- YT0 

    if (nX0) {
      dLdalpha0 <- X0
    }
    else {
      dLdalpha0 <- NULL
    }

    if (nX){
#  traiter les Intercept_t_NPH
      dLdbeta0 <- NULL
      for(i in 1:nX){
        if ( Intercept_t_NPH[i] ){
          dLdbeta0 <- cbind(dLdbeta0,  X[,i] *  YT)
        }
        else {
          dLdbeta0 <- cbind(dLdbeta0, X[,i]  * YT[,indx_without_intercept])
        }
      }
    }
    else {
      dLdbeta0 <- NULL
    }
    
    if (nZ) { 
      YTbeta <- YT %*% t(tBeta)
      indZ <- getIndex(Z)

      dLdalpha <- NULL
      dLdbeta <- NULL
      
      for(iZ in 1:nZ){
        dLdalpha <- cbind( dLdalpha, Z@DM[,indZ[iZ,1]:indZ[iZ,2]] * YTbeta[,iZ] )
        dLdbeta <- cbind(dLdbeta, YT[,-1, drop=FALSE] * Zalpha[,iZ])
      }
    }
    else {
      dLdalpha <- NULL
      dLdbeta <- NULL
    }

  rep <- cbind(dLdgamma0,          
           dLdalpha0,          
           dLdbeta0,          
           dLdalpha,          
           dLdbeta )

  rep
  
}
