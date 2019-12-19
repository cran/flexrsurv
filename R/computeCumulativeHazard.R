# when dim(Y)==2, Y[,1] is timeToEvent, Y[,2] is event
.computeCumulativeHazard_GA0B0AB<-function(GA0B0AB, Y, X0, X, Z, 
                      step, Nstep,
                      intTD=intTD_NC, intweightsfunc=intweights_CAV_SIM,
                      nT0basis,
                      Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      ialpha0, nX0,
                      ibeta0, nX,
                      ialpha, ibeta,                             
                      nTbasis,
                      Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                      Intercept_t_NPH=rep(TRUE, nX), 
                      debug=FALSE,  ...){
  # compute the cumulative hazard frm 0 to Y[,1]
  #################################o################################################################################
  #################################################################################################################
  #  the coef of the first t-basis is constraint to 1 for nat-spline, and n-sum(other beta) if bs using expand() method
  #################################################################################################################
  #################################################################################################################
  #################################################################################################################
  # GA0B0AB ; vector of all coefs
                                        # gamma0 = GA0B0AB[1:nY0basis]
                                        # alpha0= GA0B0AB[ialpha0]
                                        # beta0= matrix(GA0B0AB[ibeta0], ncol=nX, nrow=nTbasis)
                                        # alpha= diag(GA0B0AB[ialpha])
                                        # beta= expand(matrix(GA0B0AB[ibeta], ncol=Z@nZ, nrow=nTbasis-1))
                                        # beta does not contains coef for the first t-basis
  #################################################################################################################
  # Y : object of class Surv
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
  # expected_rate : expected rate at event time T
  # step : object of class "NCLagParam" or "GLMLagParam"
  # intTD : function to perform numerical integration 
  # intweightfunc : function to compute weightsfor numerical integration
  # nT0basis : number of spline basis 
  #  Spline_t0, spline object for baseline hazard, with evaluate() méthod
  #  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
  # nTbasis : number of time spline basis for NPH or NLL effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  #  Spline_t, spline object for time dependant effects,  with evaluate() méthod
  # Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
  # returned value : the log liikelihood of the model
  
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
    # add a row of one for the first T-basis 
    Beta <- t(ExpandAllCoefBasis(GA0B0AB[ibeta], ncol=nZ,  value=1))
    # parenthesis important for speed ?
    Zalphabeta <- Z@DM %*%( diag(GA0B0AB[ialpha]) %*% Z@signature  %*% Beta )
    if(nX) {
    # add a row of 0 for the first T-basis when !Intercept_T_NPH
      Zalphabeta <- Zalphabeta + X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                                         ncol=nX,
                                                         splinebasis=Spline_t,
                                                         expand=!Intercept_t_NPH,
                                                         value=0))
    }
  } else {
    if(nX) {
      Zalphabeta <- X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                            ncol=nX,
                                            splinebasis=Spline_t,
                                            expand=!Intercept_t_NPH,
                                            value=0))
    }
  }
  
  if(nX + nZ) {
    NPHterm <- intTD(rateTD_gamma0alphabeta, intTo=Y[,1], intToStatus=Y[,2],
                     step=step, Nstep=Nstep,
                     intweightsfunc=intweightsfunc, 
                     gamma0=GA0B0AB[1:nT0basis], Zalphabeta=Zalphabeta, 
                     Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0,
                     Spline_t = Spline_t, Intercept_t=TRUE)
  } else {
    NPHterm <- intTD(rateTD_gamma0, intTo=Y[,1], intToStatus=Y[,2], 
                     step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
                     gamma0=GA0B0AB[1:nT0basis],
                     Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0)
  }


  # contribution of non time dependant variables
  if( nX0){
    ret <- exp(X0 %*% GA0B0AB[ialpha0]) * NPHterm 
  } else {
    ret <- NPHterm 
  }
ret
}



.computeCumulativeHazard_GA0B0AB_bh<-function(GA0B0AB, Y, X0, X, Z, 
                      step, Nstep,
                      intTD=intTD_NC, intweightsfunc=intweights_CAV_SIM,
                      nT0basis,
                      Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      ialpha0, nX0,
                      ibeta0, nX,
                      ialpha, ibeta,                             
                      nTbasis,
                      Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                      Intercept_t_NPH=rep(TRUE, nX), 
                      debug=FALSE,  ...){
  # compute the cumulative hazard frm 0 to Y[,1]
  # rate = (f(t)%*%gamma) *  exp( X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
  #################################o################################################################################
  #################################################################################################################
  #################################################################################################################
  #################################################################################################################
  #################################################################################################################
  # GA0B0AB ; vector of all coefs
                                        # gamma0 = GA0B0AB[1:nY0basis]
                                        # alpha0= GA0B0AB[ialpha0]
                                        # beta0= matrix(GA0B0AB[ibeta0], ncol=nX, nrow=nTbasis)
                                        # alpha= diag(GA0B0AB[ialpha])
                                        # beta= expand(matrix(GA0B0AB[ibeta], ncol=Z@nZ, nrow=nTbasis-1))
                                        # beta does not contains coef for the first t-basis
  #################################################################################################################
  # Y : object of class Surv
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
  # expected_rate : expected rate at event time T
  # step : object of class "NCLagParam" or "GLMLagParam"
  # intTD : function to perform numerical integration 
  # intweightfunc : function to compute weightsfor numerical integration
  # nT0basis : number of spline basis 
  #  Spline_t0, spline object for baseline hazard, with evaluate() méthod
  #  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
  # nTbasis : number of time spline basis for NPH or NLL effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  #  Spline_t, spline object for time dependant effects,  with evaluate() méthod
  # Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
  # returned value : the log liikelihood of the model
  
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
    # add a row of one for the first T-basis 
    Beta <- t(ExpandAllCoefBasis(GA0B0AB[ibeta], ncol=nZ,  value=1))
    # parenthesis important for speed ?
    Zalphabeta <- Z@DM %*%( diag(GA0B0AB[ialpha]) %*% Z@signature  %*% Beta )
    if(nX) {
    # add a row of 0 for the first T-basis when !Intercept_T_NPH
      Zalphabeta <- Zalphabeta + X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                                         ncol=nX,
                                                         splinebasis=Spline_t,
                                                         expand=!Intercept_t_NPH,
                                                         value=0))
    }
  } else {
    if(nX) {
      Zalphabeta <- X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                            ncol=nX,
                                            splinebasis=Spline_t,
                                            expand=!Intercept_t_NPH,
                                            value=0))
    }
  }
  
  if(nX + nZ) {
    NPHterm <- intTD(rateTD_bh_alphabeta, intTo=Y[,1], intToStatus=Y[,2],
                     step=step, Nstep=Nstep,
                     intweightsfunc=intweightsfunc, 
                     gamma0=GA0B0AB[1:nT0basis], Zalphabeta=Zalphabeta, 
                     Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0,
                     Spline_t = Spline_t, Intercept_t=TRUE)
  } else {
#    NPHterm <- intTD(rateTD_bh_gamma0, intTo=Y[,1], intToStatus=Y[,2], 
#                     step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
#                     gamma0=GA0B0AB[1:nT0basis],
#                     Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0)
   NPHterm <- predict(integrate(Spline_t0*tmpgamma0), Y[,1], intercep=Intercept_t0) 
  }


  # contribution of non time dependant variables
  if( nX0){
    ret <- exp(X0 %*% GA0B0AB[ialpha0]) * NPHterm 
  } else {
    ret <- NPHterm 
  }
ret
}



# when dim(Y)==3, Y[,1] is beginTime, Y[,2] is endTime, Y[,3] is event
.computeCumulativeHazard_fromto_GA0B0AB<-function(GA0B0AB, Y, X0, X, Z, 
                      step, Nstep,
                      intTD=intTDft_NC, intweightsfunc=intweights_CAV_SIM,
                      nT0basis,
                      Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      ialpha0, nX0,
                      ibeta0, nX,
                      ialpha, ibeta,                             
                      nTbasis,
                      Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                      Intercept_t_NPH=rep(TRUE, nX), 
                      debug=FALSE,  ...){
  # compute the cumulative hazard frm Y[,1] to Y[,2]
  # rate = exp( f(t)%*%gamma + X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
  #################################################################################################################
  #################################################################################################################
  #  the coef of the first t-basis is constraint to 1 for nat-spline, and n-sum(other beta) if bs using expand() method
  #################################################################################################################
  #################################################################################################################
  #################################################################################################################
  # GA0B0AB ; vector of all coefs
                                        # gamma0 = GA0B0AB[1:nY0basis]
                                        # alpha0= GA0B0AB[ialpha0]
                                        # beta0= matrix(GA0B0AB[ibeta0], ncol=nX, nrow=nTbasis)
                                        # alpha= diag(GA0B0AB[ialpha])
                                        # beta= expand(matrix(GA0B0AB[ibeta], ncol=Z@nZ, nrow=nTbasis-1))
                                        # beta does not contains coef for the first t-basis
  #################################################################################################################
  # Y : object of class Surv with beginning and end of interval
  #
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
  # expected_rate : expected rate at event time T
  # step : object of class "NCLagParam" or "GLMLagParam"
  # intTD : function to perform numerical integration 
  # intweightfunc : function to compute weightsfor numerical integration
  # nT0basis : number of spline basis 
  #  Spline_t0, spline object for baseline hazard, with evaluate() méthod
  #  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
  # nTbasis : number of time spline basis for NPH or NLL effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  #  Spline_t, spline object for time dependant effects,  with evaluate() méthod
  # Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
  # returned value : the log liikelihood of the model
  

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
    # add a row of one for the first T-basis 
    Beta <- t(ExpandAllCoefBasis(GA0B0AB[ibeta], ncol=nZ,  value=1))
    # parenthesis important for speed ?
    Zalphabeta <- Z@DM %*%( diag(GA0B0AB[ialpha]) %*% Z@signature  %*% Beta )
    if(nX) {
    # add a row of 0 for the first T-basis when !Intercept_T_NPH
      Zalphabeta <- Zalphabeta + X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                                         ncol=nX,
                                                         splinebasis=Spline_t,
                                                         expand=!Intercept_t_NPH,
                                                         value=0))
    }
  } else {
    if(nX) {
      Zalphabeta <- X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                            ncol=nX,
                                            splinebasis=Spline_t,
                                            expand=!Intercept_t_NPH,
                                            value=0))
    }
  }
  
  if(nX + nZ) {
    NPHterm <- intTD(rateTD_gamma0alphabeta, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],
                     step=step, Nstep=Nstep,
                     intweightsfunc=intweightsfunc, 
                     gamma0=GA0B0AB[1:nT0basis], Zalphabeta=Zalphabeta, 
                     Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0,
                     Spline_t = Spline_t, Intercept_t=TRUE)
  } else {
    NPHterm <- intTD(rateTD_gamma0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],
                     step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
                     gamma0=GA0B0AB[1:nT0basis],
                     Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0)
  }

  # contribution of non time dependant variables
  if( nX0){
    ret <- exp(X0 %*% GA0B0AB[ialpha0]) * NPHterm 
  } else {
    ret <- NPHterm 
  }
  ret
}

.computeCumulativeHazard_fromto_GA0B0AB_bh<-function(GA0B0AB, Y, X0, X, Z, 
                      step, Nstep,
                      intTD=intTDft_NC, intweightsfunc=intweights_CAV_SIM,
                      nT0basis,
                      Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      ialpha0, nX0,
                      ibeta0, nX,
                      ialpha, ibeta,                             
                      nTbasis,
                      Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                      Intercept_t_NPH=rep(TRUE, nX), 
                      debug=FALSE,  ...){
  # compute the cumulative hazard frm Y[,1] to Y[,2]
  # rate = (f(t)%*%gamma) *  exp( X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))
  #################################o################################################################################
  #################################################################################################################
  #################################################################################################################
  #################################################################################################################
  #################################################################################################################
  # GA0B0AB ; vector of all coefs
                                        # gamma0 = GA0B0AB[1:nY0basis]
                                        # alpha0= GA0B0AB[ialpha0]
                                        # beta0= matrix(GA0B0AB[ibeta0], ncol=nX, nrow=nTbasis)
                                        # alpha= diag(GA0B0AB[ialpha])
                                        # beta= expand(matrix(GA0B0AB[ibeta], ncol=Z@nZ, nrow=nTbasis-1))
                                        # beta does not contains coef for the first t-basis
  #################################################################################################################
  # Y : object of class Surv
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
  # expected_rate : expected rate at event time T
  # step : object of class "NCLagParam" or "GLMLagParam"
  # intTD : function to perform numerical integration 
  # intweightfunc : function to compute weightsfor numerical integration
  # nT0basis : number of spline basis 
  #  Spline_t0, spline object for baseline hazard, with evaluate() méthod
  #  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
  # nTbasis : number of time spline basis for NPH or NLL effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  #  Spline_t, spline object for time dependant effects,  with evaluate() méthod
  # Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
  # returned value : the log liikelihood of the model
  
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
    # add a row of one for the first T-basis 
    Beta <- t(ExpandAllCoefBasis(GA0B0AB[ibeta], ncol=nZ,  value=1))
    # parenthesis important for speed ?
    Zalphabeta <- Z@DM %*%( diag(GA0B0AB[ialpha]) %*% Z@signature  %*% Beta )
    if(nX) {
    # add a row of 0 for the first T-basis when !Intercept_T_NPH
      Zalphabeta <- Zalphabeta + X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                                         ncol=nX,
                                                         splinebasis=Spline_t,
                                                         expand=!Intercept_t_NPH,
                                                         value=0))
    }
  } else {
    if(nX) {
      Zalphabeta <- X %*% t(ExpandCoefBasis(GA0B0AB[ibeta0],
                                            ncol=nX,
                                            splinebasis=Spline_t,
                                            expand=!Intercept_t_NPH,
                                            value=0))
    }
  }
  
  if(nX + nZ) {
    NPHterm <- intTD(rateTD_bh_alphabeta, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],
                     step=step, Nstep=Nstep,
                     intweightsfunc=intweightsfunc, 
                     gamma0=GA0B0AB[1:nT0basis], Zalphabeta=Zalphabeta, 
                     Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0,
                     Spline_t = Spline_t, Intercept_t=TRUE)
  } else {
#    NPHterm <- intTD(rateTD_gamma0_bh, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],
#                     step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
#                     gamma0=GA0B0AB[1:nT0basis],
#                     Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0)
   NPHterm <- predict(integrate(Spline_t0*tmpgamma0), Y[,2], intercep=Intercept_t0) -
              predict(integrate(Spline_t0*tmpgamma0), Y[,1], intercep=Intercept_t0)
  }


  # contribution of non time dependant variables
  if( nX0){
    ret <- exp(X0 %*% GA0B0AB[ialpha0]) * NPHterm 
  } else {
    ret <- NPHterm 
  }
ret
}

# with WCE effect
# when dim(Y)==3, Y[,1] is beginTime, Y[,2] is endTime, Y[,3] is event
.computeCumulativeHazard_fromto_GA0B0ABE0<-function(GA0B0ABE0, Y, X0, X, Z,  W, 
                                                    Id, FirstId, LastId,
                                                    step, Nstep,
                                                    intTD=intTDft_NC, intweightsfunc=intweights_CAV_SIM,
                                                    nT0basis,
                                                    Spline_t0=BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                                                    ialpha0, nX0,
                                                    ibeta0, nX,
                                                    ialpha, ibeta,                             
                                                    nTbasis,
                                                    ieta0, iWbeg, iWend, nW, 
                                                    Spline_t =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                                                    Intercept_t_NPH=rep(TRUE, nX), 
                                                    ISpline_W =BSplineBasis(knots=NULL,  degree=3,   keep.duplicates=TRUE),
                                                    Intercept_W=TRUE,
                                                    debug=FALSE,  ...){
  # compute the cumulative hazard frm Y[,1] to Y[,2]
  # rate = exp( f(t)%*%gamma + X0%*%alpha0 + X%*%beta0(t) + sum( alphai(zi)betai(t) ))+ sum ( wce(Wi , eta0i)(t))
  #################################################################################################################
  #################################################################################################################
  #  the coef of the first t-basis is constraint to 1 for nat-spline, and n-sum(other beta) if bs using expand() method
  #################################################################################################################
  #################################################################################################################
  #################################################################################################################
  # GA0B0ABE0 ; vector of all coefs
                                        # gamma0 = GA0B0ABE0[1:nY0basis]
                                        # alpha0= GA0B0ABE0[ialpha0]
                                        # beta0= matrix(GA0B0ABE0[ibeta0], ncol=nX, nrow=nTbasis)
                                        # alpha= diag(GA0B0ABE0[ialpha])
                                        # beta= expand(matrix(GA0B0ABE0[ibeta], ncol=Z@nZ, nrow=nTbasis-1))
                                        # beta does not contains coef for the first t-basis
                                        # eta0 = GA0B0ABE0E0[ieta0]
  #################################################################################################################
  # Y : object of class Surv but the matrix has 4 columns :
  # Y[,1] beginning(1) , fromT
  # Y[,2] end(2), toT,
  # Y[,3] status(3) intToStatus
  # Y[,4] end of followup(4) 
  #     end of followup is assumed constant by Id
  # X0 : non-time dependante variable (may contain spline bases expended for non-loglinear terms)
  # X : log lineair but time dependante variable 
  # Z : object of class "DesignMatrixNPHNLL" time dependent variables (spline basis expended)
   # W : Exposure variables used in Weighted Cumulative Exposure Models
  # Id : varibale indicating individuals Id, lines with the same Id are considered to be from the same individual
  # FirstId : all lines in FirstId[iT]:iT in the data comes from the same individual 
 # expected_rate : expected rate at event time T
  # step : object of class "NCLagParam" or "GLMLagParam"
  # intTD : function to perform numerical integration 
  # intweightfunc : function to compute weightsfor numerical integration
  # nT0basis : number of spline basis 
  #  Spline_t0, spline object for baseline hazard, with evaluate() méthod
  #  Intercept_t0=FALSE, option for evaluate, = TRUE all the basis, =FALSE all but first basis 
  # nTbasis : number of time spline basis for NPH or NLL effects
  # nX0   : nb of PH variables dim(X0)=c(nobs, nX0)
  # nX    : nb of NPHLIN variables dim(X)=c(nobs, nX)
  #  Spline_t, spline object for time dependant effects,  with evaluate() méthod
  # Intercept_t_NPH vector of intercept option for NPH spline (=FALSE when X is NLL too, ie in case of remontet additif NLLNPH)
  # nW    : nb of WCE variables dim(W)=c(nobs, nW)
  # iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
  #  ISpline_W, list of nW spline object for WCE effects,  with evaluate() méthod
  #  ... not used args
  # the function do not check the concorcance between length of parameter vectors and the number of knots and the Z.signature
  # returned value : the log liikelihood of the model
  

  if(is.null(Z)){
    nZ <- 0
  } else {
    nZ <- Z@nZ
  }

  if(Intercept_t0){
    tmpgamma0 <- GA0B0ABE0[1:nT0basis]
  }
  else {
    tmpgamma0 <- c(0, GA0B0ABE0[1:nT0basis])
  }

  # contribution of time d?pendant effect
  # parenthesis are important for efficiency
  if(nZ) {
    # add a row of one for the first T-basis 
    Beta <- t(ExpandAllCoefBasis(GA0B0ABE0[ibeta], ncol=nZ,  value=1))
    # parenthesis important for speed ?
    Zalphabeta <- Z@DM %*%( diag(GA0B0ABE0[ialpha]) %*% Z@signature  %*% Beta )
    if(nX) {
    # add a row of 0 for the first T-basis when !Intercept_T_NPH
      Zalphabeta <- Zalphabeta + X %*% t(ExpandCoefBasis(GA0B0ABE0[ibeta0],
                                                         ncol=nX,
                                                         splinebasis=Spline_t,
                                                         expand=!Intercept_t_NPH,
                                                         value=0))
    }
  } else {
    if(nX) {
      Zalphabeta <- X %*% t(ExpandCoefBasis(GA0B0ABE0[ibeta0],
                                            ncol=nX,
                                            splinebasis=Spline_t,
                                            expand=!Intercept_t_NPH,
                                            value=0))
    }
  }
  
  if(nW){
    IS_W<- ISpline_W
    eta0 <- GA0B0ABE0[ieta0]
    for(iW in 1:nW){
      if(Intercept_W[[iW]]){
        IS_W[[iW]] <- ISpline_W[[iW]] * eta0[iWbeg[iW]:iWend[iW]]
      }
      else {
        IS_W[[iW]]<- ISpline_W[[iW]] * c(0, eta0[iWbeg[iW]:iWend[iW]])
      }
    }
    if(nX + nZ) {
      NPHterm <- intTD(rateTD_gamma0alphabetaeta0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],
                       step=step, Nstep=Nstep,
                       intweightsfunc=intweightsfunc, 
                       fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
                       gamma0=GA0B0ABE0[1:nT0basis], Zalphabeta=Zalphabeta,
                       nW = nW, W = W, eta0=GA0B0ABE0[ieta0], iWbeg=iWbeg, iWend=iWend,
                       Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0,
                       Spline_t = Spline_t, Intercept_t=TRUE,
                       ISpline_W = IS_W)
    } else {
      NPHterm <- intTD(rateTD_gamma0eta0, intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],
                       step=step, Nstep=Nstep,
                       intweightsfunc=intweightsfunc, 
                       fromT=Y[,1], toT=Y[,2], FirstId=FirstId, LastId=LastId,
                       nW = nW, W = W, eta0=GA0B0ABE0[ieta0], iWbeg=iWbeg, iWend=iWend,
                       Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0,
                       ISpline_W = IS_W)
    }
  }
  else {
    # no VCE effect, same NPH term than ll_flexrsurv_fromto_GA0B0ABE0
    if(nX + nZ) {
      NPHterm <- intTD(rateTD_gamma0alphabeta,  intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],
                       step=step, Nstep=Nstep,
                       intweightsfunc=intweightsfunc, 
                       gamma0=GA0B0ABE0[1:nT0basis], Zalphabeta=Zalphabeta, 
                       Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0,
                       Spline_t = Spline_t, Intercept_t=TRUE)
    } else {
      NPHterm <- intTD(rateTD_gamma0,  intFrom=Y[,1], intTo=Y[,2], intToStatus=Y[,3],
                       step=step, Nstep=Nstep, intweightsfunc=intweightsfunc, 
                       gamma0=GA0B0ABE0[1:nT0basis],
                       Spline_t0=Spline_t0*tmpgamma0, Intercept_t0=Intercept_t0)
    }
  }
  # contribution of non time dependant variables
  if( nX0){
    ret <- exp(X0 %*% GA0B0ABE0[ialpha0]) * NPHterm 
  } else {
    ret <- NPHterm 
  }
  ret
}











