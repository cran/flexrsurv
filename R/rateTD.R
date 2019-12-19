# the spline basis of the baseline hazard are scaled by the parameter gamma0,
# the spline basis of the NPH NPHNLL effect are the normalised spline basis

rateTD_beta0alphabeta<- function(T, iT, gamma0, Zbeta0, Zalphabeta, 
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE, ...){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
  # spline bases for baseline hazard
  # all T basis for the NPH/td effects are the same (Spline_t)
  # Zbeta0     = X %*% beta0 
  # Zalphabeta = f(Z,alpha) %*% beta 


  # predicted log-baseline hazard
    YT0Gamma0 <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
  # spline bases for each TD effect
#  YT  <- bs(T, knots=Knots_t, intercept=Intercept_t, degree=degree_t, Boundary.knots =  Boundary.knots_t)
    if(!is.null(Zbeta0)) {
      if(!is.null(Zalphabeta)) {
        exp(YT0Gamma0 + 
            fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] +
            fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
      }
      else {
        exp(YT0Gamma0 + 
            fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zbeta0[iT,] )
      }
    }
    else if(!is.null(Zalphabeta)) {
        exp(YT0Gamma0 + 
            fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE) %*% Zalphabeta[iT,])
      }
}

rateTD_gamma0alphabeta<- function(T, iT, gamma0, Zalphabeta, 
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE, ...){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
   # all T basis for the NPH/td effects are the same (Spline_t)
  # Zalphabeta = X %*% beta0 + f(Z,alpha) %*% beta 
  # Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
  #           : thus no nead to multiply each spline coordinate by its coef

  # predicted log-baseline hazard
    YT0Gamma0 <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
  # spline bases for each TD effect
    YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)

  # returned value
  exp(YT0Gamma0 + YT %*% Zalphabeta[iT,])
  
}

rateTD_alphabeta<- function(T, iT, Zalphabeta, 
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE, ...){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
  # all T basis for the NPH/td effects are the same (Spline_t)
  # Zalphabeta = X %*% beta0 + f(Z,alpha) %*% beta 
 
  # spline bases for each TD effect
    YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)

  # returned value 
  exp(YT %*% Zalphabeta[iT,])
  
}

rateTD_gamma0<- function(T, iT, gamma0, 
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE, ...){
  # compute the contribution of the baseline hazard rate to the rate 
  # of relative survival model for patient iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration 
  # Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
  #           : thus no nead to multiply each spline coordinate by its coef

  # spline bases for baseline hazard
    YT0Gamma0 <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
  # returned value

  exp(YT0Gamma0)
  
}


rateTD<- function(T, iT, ... ){
  # compute the contribution of the baseline hazard rate when there are no time dependent effect
  # and when the baseline hazard is null
 
  # returned value

  rep(1, length(T))
  
}



# computes the contribution of time dependent termes in the rate (baseline, NPH, NPHNLL and WCE effects 
rateTD_gamma0alphabetaeta0<- function(T, iT,
                                      fromT, toT, FirstId, LastId,
                                      gamma0, Zalphabeta,
                                      nW, W, eta0, iWbeg, iWend,
                      Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                      Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE,
                                      ISpline_W, Intercept_W=rep(TRUE,nW), ...){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for line iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration
  # fromT : begining of the time intervals
  # toT   : end of the time intervals
  # all T basis for the NPH/td effects are the same (Spline_t)
  # Zalphabeta = X %*% beta0 + f(Z,alpha) %*% beta 
  # FirstId : all lines in FirstId[iT]:iT of fromT, toT, and Zalphabeta comes from the same individual 
  # nW number of cols in W (number of WCE effects
  # W matrix of exposure INCREMENT variables W[FirstId[iT]:iT, k] is the vectore of exposure increment x_il - x_(i-1)l for patient l, expo variable k    
  # eta0 : vector all the coef of WCE
  # iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
  # ISpline_W : list of the nW integrated splines parameters for the WCE effects scaled by eta0
  #           : thus no nead to multiply each spline coordinate by its coef
  # Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
  #           : thus no nead to multiply each spline coordinate by its coef

#  print("rateTD_gamma0alphabetaeta0")
  # spline bases for baseline hazard
    WCE <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
  # spline bases for each TD effect
    YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)


    for(iW in 1:nW){
      wce <- predictwce(object=ISpline_W[[iW]], t=T, Increment=W[,iW], fromT=fromT, tId=rep(iT, length(T)),
                           FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
      WCE = WCE + wce

#      for(iId in FirstId[iT]:iT){
#        WCE <- WCE + W[iId, iW] * predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)  
#      }
    }
  # returned value 
  exp(WCE + YT %*% Zalphabeta[iT,])
  
}



# computes the contribution of gamma0 and eta0 (baseline & WCE)
rateTD_gamma0eta0<- function(T, iT,
                             fromT, FirstId, LastId,
                             gamma0, 
                             nW, W, 
                             Spline_t0=SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t0=TRUE,
                             ISpline_W, Intercept_W=rep(TRUE,nW), ...){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT 
  # at a vector of T (useful to compute numerical integration
  # fromT : begining of the time intervals
  # toT   : end of the time intervals
  # FirstId : all lines in FirstId[iT]:iT of fromT, toT, and Zalphabeta comes from the same individual 
 # nW number of cols in W (number of WCE effects
  # W matrix of exposure INCREMENT variables W[FirstId[iT]:iT, k] is the vectore of exposure increment x_il - x_(i-1)l for patient l, expo variable k    
  # eta0 : vector all the coef of WCE
  # iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
  # ISpline_W : list of the nW integrated splines parameters for the WCE effects SCALED by eta0
  #           : thus no nead to multiply each spline coordinate by its coef
  # Spline_t0 : splines parameters for the baseline hazard multiplied by gamma0
  #           : thus no nead to multiply each spline coordinate by its coef

#  print("rateTD_gamma0eta0")
  
      
    WCE <- predictSpline(Spline_t0, T, intercept=Intercept_t0, outer.ok=TRUE)
    for(iW in 1:nW){
      wce <- predictwce(object=ISpline_W[[iW]], t=T, Increment=W[,iW], fromT=fromT, tId=rep(iT, length(T)),
                        FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
      WCE = WCE + wce
       
#       for(iId in FirstId[iT]:iT){
#         WCE <- WCE + W[iId, iW] * predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)  
#       }
    }
  exp(WCE)
  
}



# computes the contribution eta0 (WCE)
rateTD_eta0<- function(T, iT,
                       fromT, FirstId, LastId,
                       nW, W, 
                       ISpline_W, Intercept_W=rep(TRUE,nW), ...){
  # compute the contribution of the time dependent variables to the rate 
  # of relative survival model for patient iT 
  # at a vector of T (useful to compute numerical integration
  # fromT : begining of the time intervals
  # toT   : end of the time intervals
  # FirstId : all lines in FirstId[iT]:iT of fromT, toT, and Zalphabeta comes from the same individual 
 # nW number of cols in W (number of WCE effects
  # W matrix of exposure INCREMENT variables W[FirstId[iT]:iT, k] is the vectore of exposure increment x_il - x_(i-1)l for patient l, expo variable k    
  # eta0 : vector all the coef of WCE
  # iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
  # ISpline_W : list of the nW integrated splines parameters for the WCE effects SCALED by eta0
  #           : thus no nead to multiply each spline coordinate by its coef

#  print("rateTD_eta0")
  
      
    WCE <- predictwce(object=ISpline_W[[1]], t=T, Increment=W[,1], fromT=fromT, tId=rep(iT, length(T)),
                           FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
    iW <- 2
    while(iW <=nW){
      wce <- predictwce(object=ISpline_W[[iW]], t=T, Increment=W[,iW], fromT=fromT, tId=rep(iT, length(T)),
                        FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
      WCE = WCE + wce
#      for(iId in FirstId[iT]:iT){
#        WCE <- WCE + W[iId, iW] * predictSpline(ISpline_W[[iW]], T-fromT[iId], intercept=Intercept_W[[iW]], outer.ok=TRUE)  
#      }
       iW <- iW+1
    }
  # returned value
  exp(WCE)
  
}






# computes the contribution of time dependent termes in the rate of additive proportionamle WCE model (WCE, NPH, NPHNLL)
rateTD_alphabeta_1addwce<- function(T, iT,
                                      fromT, toT, FirstId, LastId,
                                      Zalphabeta,
                                      W, 
                                        Spline_t =SplineBasis(knots=NULL,  order=4,   keep.duplicates=TRUE), Intercept_t=TRUE,
                                      ISpline_W, Intercept_W=TRUE, ...){
  # compute the contribution of the time dependent variables to the rate WCE(W,T)*exp(NPH + NPHNLL)
  # of relative survival model for line iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration
  # fromT : begining of the time intervals
  # toT   : end of the time intervals
  # all T basis for the NPH/td effects are the same (Spline_t)
  # Zalphabeta = X %*% beta0 + f(Z,alpha) %*% beta 
  # FirstId : all lines in FirstId[iT]:iT of fromT, toT, and Zalphabeta comes from the same individual 
  # W vecteur of exposure INCREMENT variables W[FirstId[iT]:iT] is the vectore of exposure increment x_il - x_(i-1)l for patient l
  # eta0 : vector all the coef of WCE
  # iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
  # ISpline_W : a single integrated splines parameters for the WCE effects scaled by eta0
  #           : thus no nead to multiply each spline coordinate by its coef

#  print("rateTD_alphabeta_1addwceeta0")
  # spline bases for baseline hazard
    WCE <- predictwce(object=ISpline_W, t=T, Increment=W, fromT=fromT, tId=rep(iT, length(T)),
                           FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
  # spline bases for each TD effect
#    YT <- fevaluate(Spline_t, T, intercept=Intercept_t, outer.ok=TRUE)
    YT <- evaluatelc(Spline_t, T, Zalphabeta[iT,], intercept=Intercept_t, outer.ok=TRUE)


  # returned value 
#  WCE * exp(YT %*% Zalphabeta[iT,])
  WCE * exp(YT)
  
}

# computes the contribution of additive proportionamle WCE model (no NPH, NPHNLL)
rateTD_1addwce<- function(T, iT,
                              fromT, toT, FirstId, LastId,
                              W, 
                              ISpline_W, Intercept_W=TRUE, ...){
  # compute the contribution of the time dependent variables to the rate WCE(W,T)*exp(NPH + NPHNLL)
  # of relative survival model for line iT with Zalphabeta[iT, ]
  # at a vector of T (useful to compute numerical integration
  # fromT : begining of the time intervals
  # toT   : end of the time intervals
  # all T basis for the NPH/td effects are the same (Spline_t)
  # Zalphabeta = X %*% beta0 + f(Z,alpha) %*% beta 
  # FirstId : all lines in FirstId[iT]:iT of fromT, toT, and Zalphabeta comes from the same individual 
  # W vecteur of exposure INCREMENT variables W[FirstId[iT]:iT] is the vectore of exposure increment x_il - x_(i-1)l for patient l
  # eta0 : vector all the coef of WCE
  # iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
  # ISpline_W : a single integrated splines parameters for the WCE effects scaled by eta0
  #           : thus no nead to multiply each spline coordinate by its coef

#  print("rateTD_alphabeta_1addwceeta0")
  # spline bases for baseline hazard
    predictwce(object=ISpline_W, t=T, Increment=W, fromT=fromT, tId=rep(iT, length(T)),
                           FirstId=FirstId, LastId=LastId, intercept=Intercept_W, outer.ok=TRUE)
  
}

