intTDft_base_NC <- function(func=function(x) return(x), intFrom, intTo,
                            Spline,
                            step, Nstep, intweightsfunc = intweights_CAV_SIM,
                            intToStatus=NULL,
                            debug=TRUE,
                            ...){
  # compute numerical integral of func*base_i(t)  in [intFrom , intTo] following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # intFrom : lower bound (vector)
  # intTo   : upper bound (vector)
  # Spline : Spline parameters
  # step : vector of the steps (one row per T)
  # Nstep : vector of the number of steps ((intTo - intFrom) = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # intToStatus : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
  for(i in 1:length(intTo)){
    # vector of evaluated t
    theT <- intFrom[i] + (0:Nstep[i])*step[i]
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, ...)
    TBase <- fevaluate(Spline, theT, intercept=TRUE)
    # weights 1 * nt matrix
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i,] <- crossprod(w*FF, TBase)
  }
  res * step 
}

intTDft_base2_NC <- function(func=function(x) return(x), intFrom, intTo, FromT,
                            Spline,
                            step, Nstep, intweightsfunc = intweights_CAV_SIM,
                            intToStatus=NULL,
                            debug=TRUE,
                            ...){
  #similar to intTDft_base_NC but
  # compute numerical integral of func*base_i(t - FromT)  in [intFrom , intTo] following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # intFrom : lower bound (vector)
  # intTo   : upper bound (vector)
  # Spline : Spline parameters
  # step : vector of the steps (one row per T)
  # Nstep : vector of the number of steps ((intTo - intFrom) = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # intToStatus : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
  for(i in 1:length(intTo)){
    # vector of evaluated t
    theT <- intFrom[i] + (0:Nstep[i])*step[i]
    # evaluate spline basis at t - intFrom
    TBase <- fevaluate(Spline, theT - FromT[i], intercept=TRUE)
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, ...)
    # weights 1 * nt matrix
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i,] <- crossprod(w*FF, TBase)
  }
#  cat("outinintTD_NC\n")
  res * step 
}


intTDft_base_GL <- function(func=function(x) x, intFrom, intTo,
                            Spline,
                            step, Nstep,
                            intweightsfunc = NULL,
                            intToStatus=NULL,
                            ...){
  # compute numerical integral of func*base_i(intFrom)  in [intFrom , intTo] Gauss Legendre quadrature
  # func : (vector of) function to integrate, func(t, ...)
  # intFrom : lower bound (vector)
  # intTo   : upper bound (vector)
  # Spline : Spline parameters
  # step : points of the quadrature
  # Nstep : weights of the quadrature
  # intweightfunc : unused
  # intToStatus : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-matrix(0, nrow = length(intFrom), ncol = Spline@nbases + Spline@log)
  Tmid <- (intTo + intFrom)/2
  dT   <- (intTo - intFrom)/2
  for(i in 1:length(intFrom)){
    # vector of evaluated t
    theT <- dT[i] * step + Tmid[i]
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, ...)
    TBase <- fevaluate(Spline, theT, intercept=TRUE)
# numerical integration
    res[i,] <- crossprod(Nstep*FF, TBase)
  }
  res * dT
  
}



intTDft_WCEbase_NC <- function(func=function(x) return(x), intFrom, intTo, 
                               Spline, intercept,
                               theW, fromT, FirstId, LastId,
                               step, Nstep, intweightsfunc = intweights_CAV_SIM,
                               debug=TRUE,
                               ...){
  #similar to intTDft_base_NC but
  # compute numerical integral_intFrom^intTo {sum_o=firstid^j  theW[o] base_i(t - fromT[o])} func(t) dt    following Newton_Cote method
  # ie      numerical integral_intFrom^intTo {sum_o=firstid^j  gradwce(t, theW, fromTo) func(t) dt    following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # intFrom : lower bound (vector)
  # intTo   : upper bound (vector)
  # fromT : begining of the time intervalle of the time-to-event exposure
  # FirstId : all lines in FirstId[iT]:iT in the data comes from the same individual 
  # Spline : integrated Spline parameters of the wce to integrate
  # intercept : =FALSE if intercept is removed
  # theW : vectot of increment of exposure 
  # step : vector of the steps (one row per T)
  # Nstep : vector of the number of steps ((intTo - intFrom) = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # intToStatus : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
  for(i in 1:length(intTo)){
    # vector of evaluated t
    theT <- intFrom[i] + (0:Nstep[i])*step[i]
    TBase <- 0
    # gradient of WCE at theT
    TBase <- gradientwce(object=Spline, t=theT, Increment=theW, fromT=fromT, tId=rep(i, Nstep[i]+1),
                           FirstId=FirstId, LastId=LastId, intercept=intercept, outer.ok=TRUE)
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, FirstId=FirstId, ...)
    # weights 1 * nt matrix
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i,] <- crossprod(w*FF, TBase)

  }
#  cat("outinintTD_NC\n")
  res * step 
}

intTDft_WCEbase_GL <- function(func=function(x) return(x), intFrom, intTo,
                               Spline, intercept,
                               theW, fromT, toT, FirstId, LastId,
                               step, Nstep, intweightsfunc = intweights_CAV_SIM,
                               debug=TRUE,
                               ...){
  #similar to intTDft_base_LG but
  # compute numerical integral_intFrom^intTo {sum_o=firstid^j  theW[o] base_i(t - fromT[o])} func(t) dt    following Gauss Legendre quadrature
  # ie      numerical integral_intFrom^intTo {sum_o=firstid^j  gradwce(t, theW, fromTo) func(t) dt    following  Gauss Legendre quadrature
  # func : (vector of) function to integrate, func(t, ...)
  # intFrom : lower bound (vector)
  # intTo   : upper bound (vector)
  # FirstId : all lines in FirstId[iT]:iT in the data comes from the same individual 
  # Spline : integrated Spline parameters of the wce to integrate
  # intercept : =FALSE if intercept is removed
  # theW : vectot of increment of exposure 
  # step : points of the quadrature
  # Nstep : weights of the quadrature
  # intweightfunc : unused
  # intToStatus : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
  Tmid <- (intTo + intFrom)/2
  dT   <- (intTo - intFrom)/2
  npoints <- length(step)
  for(i in 1:length(intTo)){
    # vector of evaluated t
    theT <- dT[i] * step + Tmid[i]
    TBase <- 0
    # gradient of WCE at theT
    TBase <- gradientwce(object=Spline, t=theT, Increment=theW, fromT=fromT, tId=rep(i, npoints),
                           FirstId=FirstId, LastId=LastId, intercept=intercept, outer.ok=TRUE)
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, FirstId=FirstId, ...)

# numerical integration
    res[i,] <- crossprod(Nstep*FF, TBase)
  }
  dT * res  
}

slowintTDft_WCEbase_NC <- function(func=function(x) return(x), intFrom, intTo, 
                                   Spline, intercept,
                                   theW, fromT, toT, FirstId,
                                   step, Nstep, intweightsfunc = intweights_CAV_SIM,
                                   intToStatus=NULL,
                                   debug=TRUE,
                            ...){
  #similar to intTDft_base_NC but
  # compute sum_o=firstid^j  theW[o] numerical integral_intFrom^intTo func(t) *base_i(t - fromT[o])   following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # intFrom : lower bound (vector)
  # intTo   : upper bound (vector)
  # FirstId : all lines in FirstId[iT]:iT in the data comes from the same individual 
  # Spline : Spline parameters of the base to integrate
  # intercept : =FALSE if intercept is removed
  # theW : vectot of increment of exposure 
  # step : vector of the steps (one row per T)
  # Nstep : vector of the number of steps ((intTo - intFrom) = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # intToStatus : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
  for(i in 1:length(intTo)){
    # vector of evaluated t
    theT <- intFrom[i] + (0:Nstep[i])*step[i]
    TBase <- 0
    for(iId in FirstId[i]:i){
      # evaluate spline basis at t - fromT[o]
      TBase <- TBase +  theW[iId] * fevaluate(Spline, theT-fromT[i], intercept=intercept)
    }
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, FirstId=FirstId, ...)
    # weights 1 * nt matrix
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i,] <- crossprod(w*FF, TBase)

  }
#  cat("outinintTD_NC\n")
  res * step 
}

intTDft_base_NC_debug<- function(func=function(x) return(x), intFrom, intTo,
                                 Spline,
                                 step, Nstep, intweightsfunc = intweights_CAV_SIM,
                                 intToStatus=NULL,
                                 debug=TRUE,
                                 ...){
  # compute numerical integral of func*base_i(t)  in [intFrom , intTo] following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # intFrom : lower bound (vector)
  # intTo   : upper bound (vector)
  # Spline : Spline parameters
  # step : vector of the steps (one row per T)
  # Nstep : vector of the number of steps ((intTo - intFrom) = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # intToStatus : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  cat("inintTD_NC_debug\n")
  cat("lengthT lengthNstep step \n")
  cat(length(T), length(Nstep), length(step))
  cat("\n")
  print(cbind(T,Nstep, step)[1:20,])
  cat("\n")
  res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
  for(i in 1:length(intTo)){
    # vector of evaluated t
    theT <- intFrom[i] + (0:Nstep[i])*step[i]
    TBase <- fevaluate(Spline, theT, intercept=TRUE)
    # vector of the evaluated functions
    FF <- func(theT, i, ...)
    # weights
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i,] <- crossprod(w*FF, TBase)
  }
  cat("outinintTD_NC\n")
  res * step 
}



fastintTDft_base_GLM <- function(func=function(x) return(x), intFrom, intTo,
                                 Spline,
                               step, Nstep, intweightsfunc=NULL,
                               intToStatus,
                               debug=FALSE,
                               ...){
  # compute numerical integral of func*b_i(t) in [intFrom , intTo] for equivalence with the poisson GLM trick
  # func : function to integrate, func(t, ...)
  # intFrom : lower bound (vector)
  # intTo   : upper bound (vector)
  # Spline : Spline parameters
  # step : object of class GLMStepParam
  # Nstep : index of the first and last complete band ( intFrom[i] < step@cuts[Nstep[i,1]] <= step@cuts[Nstep[i,2]+1] < intTo)
  #                                                    ( intFrom[i] < step@points[Nstep[i,1]] <= step@points[Nstep[i,2]] < intTo)
  # intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
  # intToStatus : status at intTo
  # ... : parameters of func()
  if(debug>200) {
    cat("fastinintTD_base_glm\n")
  }
  res<-matrix(0, nrow = length(intFrom), ncol = Spline@nbases+Spline@log)
                     # matrix of bases evaluated at the points and T
  allTBase <- fevaluate(Spline,step@points , intercept=TRUE)
  Tpoints <- ifelse(intToStatus, intTo,  (step@cuts[1+Nstep[,2]]+intTo)/2) 
  TBaseatintTo <- fevaluate(Spline, Tpoints , intercept=TRUE)
  TBaseatintFrom <- fevaluate(Spline, (step@cuts[Nstep[,1]]+intFrom)/2 , intercept=TRUE)
  for(i in 1:length(intTo)){
    # vector of evaluated t
    if(Nstep[i,2]>= Nstep[i,1]){
      # at least one complete step
        theT <- c((step@cuts[Nstep[i,1]]+intFrom[i])/2,
                  step@points[Nstep[i,1]:Nstep[i,2]] ,
                  Tpoints[i] )
                                        # vector of the evaluated functions
        FF <- func(theT, i, ...)
                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - intFrom[i],
            step@steps[Nstep[i,1]:Nstep[i,2]],
            intTo[i]-step@cuts[1+Nstep[i,2]])
                                        # numerical integration of the complete bands
        res[i,] <- crossprod(w*FF,
                             rbind( TBaseatintFrom[i,],
                                   allTBase[Nstep[i,1]:Nstep[i,2],, drop=FALSE],
                                   TBaseatintTo[i,]))
      }
    else if(Nstep[i,2] - Nstep[i,1] == -1L){
# intFrom and intTo are in 2 successive bands
      # Nstep[i,2] + 1 =  Nstep[i,1]
        theT <- c((step@cuts[Nstep[i,1]]+intFrom[i])/2,
                  Tpoints[i] )
                                        # vector of the evaluated functions
        FF <- func(theT, i, ...)
                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - intFrom[i],
            intTo[i]-step@cuts[Nstep[i,1]])
                                        # numerical integration of the complete bands
        res[i,] <- crossprod(w*FF,
                             rbind( TBaseatintFrom[i,], TBaseatintTo[i,]))
      }
    else { #if((Nstep[i,2] - Nstep[i,1]) == -2L){
# intFrom and intTo are in the same band
      if( intToStatus[i]!=0 ){
        res[i,] <- (intTo[i]- intFrom[i]) * func(intTo[i], i, ...) * TBaseatintTo[i,]
      }
      else {
        res[i,] <-  ((intTo[i] - intFrom[i]) * func((intTo[i] + intFrom[i])/2, i,  ...)) %*% fevaluate(Spline, (intTo[i] + intFrom[i])/2 , intercept=TRUE) #[,,drop=TRUE] 
      }
    }
  }
  res  
}



fastintTDft_base2_GLM <- function(func=function(x) return(x), intFrom, intTo, fromT, toT,
                                 Spline,
                               step, Nstep, intweightsfunc=NULL,
                               intToStatus,
                               debug=FALSE,
                               ...){
  #similar to fastintTDft_base2_GLM but
  # compute numerical integral of func*b_i(t-FromT) in [intFrom , intTo] for equivalence with the poisson GLM trick
  # func : function to integrate, func(t, ...)
  # intFrom : lower bound (vector)
  # intTo   : upper bound (vector)
  # Spline : Spline parameters
  # step : object of class GLMStepParam
  # Nstep : index of the first and last complete band ( intFrom[i] < step@cuts[Nstep[i,1]] <= step@cuts[Nstep[i,2]+1] < intTo)
  #                                                    ( intFrom[i] < step@points[Nstep[i,1]] <= step@points[Nstep[i,2]] < intTo)
  # intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
  # intToStatus : status at intTo
  # ... : parameters of func()
  if(debug>200) {
    cat("fastinintTD_base_glm\n")
  }
  res<-matrix(0, nrow = length(intFrom), ncol = Spline@nbases+Spline@log)
                     # matrix of bases evaluated at the points and T
  Tpoints <- ifelse(intToStatus, intTo,  (step@cuts[1+Nstep[,2]]+intTo)/2) 
  TBaseatintTo <- fevaluate(Spline, Tpoints - fromT, intercept=TRUE)
  TBaseatintFrom <- fevaluate(Spline, (step@cuts[Nstep[,1]] - fromT)/2 , intercept=TRUE)
  for(i in 1:length(intTo)){
    # vector of evaluated t
    if(Nstep[i,2]>= Nstep[i,1]){
      # at least one complete step
        theT <- c((step@cuts[Nstep[i,1]]+intFrom[i])/2,
                  step@points[Nstep[i,1]:Nstep[i,2]] ,
                  Tpoints[i] )
        #evaluated bases 
    allTBase <- fevaluate(Spline, theT - fromT[i], intercept=TRUE)
                                        # vector of the evaluated functions
        FF <- func(theT, i, ...)

                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - intFrom[i],
            step@steps[Nstep[i,1]:Nstep[i,2]],
            intTo[i]-step@cuts[1+Nstep[i,2]])
                                        # numerical integration of the complete bands
        res[i,] <- crossprod(w*FF, allTBase)
      }
    else if(Nstep[i,2] - Nstep[i,1] == -1L){
# intFrom and intTo are in 2 successive bands
      # Nstep[i,2] + 1 =  Nstep[i,1]
        theT <- c((step@cuts[Nstep[i,1]]+intFrom[i])/2,
                  Tpoints[i] )
        #evaluated bases 
        allTBase <- fevaluate(Spline, theT - fromT[i], intercept=TRUE)
                                        # vector of the evaluated functions
        FF <- func(theT, i, ...)
                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - intFrom[i],
            intTo[i]-step@cuts[Nstep[i,1]])
                                        # numerical integration of the complete bands
        res[i,] <- crossprod(w*FF, allTBase)
      }
    else { #if((Nstep[i,2] - Nstep[i,1]) == -2L){
# intFrom and intTo are in the same band
      if( intToStatus[i]!=0 ){
        res[i,] <- (intTo[i]- intFrom[i]) * func(intTo[i], i, ...) %*% fevaluate(Spline, intTo[i] - fromT[i], intercept=TRUE) #[,,drop=TRUE]
      }
      else {
        res[i,] <- ((intTo[i] - intFrom[i]) * func((intTo[i] + intFrom[i])/2, i,  ...)) %*% fevaluate(Spline, (intTo[i] - fromT[i])/2 , intercept=TRUE) #[,,drop=TRUE] 
      }
    }
  }
  res  
}



fastintTDft_WCEbase_GLM <- function(func=function(x) return(x), intFrom, intTo, fromT, toT, FirstId,
                                    Spline, intercept, theW, 
                                    step, Nstep, intweightsfunc=NULL,
                                    intToStatus,
                                    debug=FALSE,
                                    ...){
  #similar to fastintTDft_base2_GLM but
  # compute sum_o=firstid^j  theW[o] numerical integral_intFrom^intTo func(t) *base_i(t - fromT[o]) for equivalence with the poisson GLM trick
  # func : function to integrate, func(t, ...)
  # intFrom : lower bound (vector)
  # intTo   : upper bound (vector)
  # FirstId : all lines in FirstId[iT]:iT in the data comes from the same individual 
  # Spline : Spline parameters of the base to integrate
  # intercept : =FALSE if intercept is removed
  # theW : vectot of increment of exposure 
  # step : object of class GLMStepParam
  # Nstep : index of the first and last complete band ( intFrom[i] < step@cuts[Nstep[i,1]] <= step@cuts[Nstep[i,2]+1] < intTo)
  #                                                    ( intFrom[i] < step@points[Nstep[i,1]] <= step@points[Nstep[i,2]] < intTo)
  # intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
  # intToStatus : status at intTo
  # ... : parameters of func()
  if(debug>200) {
    cat("fastinintTD_base_glm\n")
  }
  res<-matrix(0, nrow = length(intFrom), ncol = Spline@nbases+Spline@log)
                     # matrix of bases evaluated at the points and T
  Tpoints <- ifelse(intToStatus, intTo,  (step@cuts[1+Nstep[,2]]+intTo)/2) 
  TBaseatintTo <- fevaluate(Spline, Tpoints - fromT, intercept=intercept)
  TBaseatintFrom <- fevaluate(Spline, (step@cuts[Nstep[,1]] - fromT)/2 , intercept=intercept)
  for(i in 1:length(intTo)){
    # vector of evaluated t
    if(Nstep[i,2]>= Nstep[i,1]){
      # at least one complete step
        theT <- c((step@cuts[Nstep[i,1]]+intFrom[i])/2,
                  step@points[Nstep[i,1]:Nstep[i,2]] ,
                  Tpoints[i] )
        #evaluated bases 
    allTBase <- 0
    for(iId in FirstId[i]:i){
      # evaluate spline basis at t - fromT[o]
      allTBase <- allTBase +  theW[iId] * fevaluate(Spline, theT-fromT[iId], intercept=intercept)
    }
                                        # vector of the evaluated functions
        FF <- func(theT, i, FirstId=FirstId, ...)

        print("theT - fromT")
#        print(theT-fromT[i])
        print("rate FF")
#        print(log(FF))
        print("allbase")
        print(cbind(theT-fromT[i], log(FF), allTBase))
                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - intFrom[i],
            step@steps[Nstep[i,1]:Nstep[i,2]],
            intTo[i]-step@cuts[1+Nstep[i,2]])
                                        # numerical integration of the complete bands
        res[i,] <- crossprod(w*FF, allTBase)
      }
    else if(Nstep[i,2] - Nstep[i,1] == -1L){
# intFrom and intTo are in 2 successive bands328.4213
      # Nstep[i,2] + 1 =  Nstep[i,1]
        theT <- c((step@cuts[Nstep[i,1]]+intFrom[i])/2,
                  Tpoints[i] )
        #evaluated bases 
    allTBase <- 0
    for(iId in FirstId[i]:i){
      # evaluate spline basis at t - fromT[o]
      allTBase <- allTBase +  theW[iId] * fevaluate(Spline, theT-fromT[iId], intercept=intercept)
    }
                                        # vector of the evaluated functions
        FF <- func(theT, i, ...)
                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - intFrom[i],
            intTo[i]-step@cuts[Nstep[i,1]])
                                        # numerical integration of the complete bands
        res[i,] <- crossprod(w*FF, allTBase)
      }
    else { #if((Nstep[i,2] - Nstep[i,1]) == -2L){
# intFrom and intTo are in the same band
      if( intToStatus[i]!=0 ){
        theT <- intTo[i]
      }
      else {
        theT <- (intTo[i] + intFrom[i])/2
      }
        #evaluated bases 
      allTBase <- 0
      for(iId in FirstId[i]:i){
      # evaluate spline basis at t - fromT[o]
        allTBase <- allTBase +  theW[iId] * fevaluate(Spline, theT-fromT[iId], intercept=intercept)
      }
      res[i,] <- (intTo[i]- intFrom[i]) * func(intTo[i], i, ...) %*% allTBase #[,,drop=TRUE]
    }
  }
  res  
}







