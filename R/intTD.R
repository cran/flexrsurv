# numerical integration in [0, intTo], of functions with arguments : func(t, i, ...)
# where i is th index 

intTD_NC <- function(func=function(x) return(x), intTo,
                     step, Nstep,
                     intweightsfunc = intweights_CAV_SIM,
                     intToStatus=NULL,
                     debug=FALSE,
                     ...){
  # compute numerical integral of func in [0 , intTo] following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # intTo    : upper bound (vector)
  # step : vector of the steps (one row per intTo)
  # Nstep : vector of the number of steps (T = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # intToStatus : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-vector("numeric", length(intTo))
  for(i in 1:length(intTo)){
    # vector of evaluated t
    theT <- (0:Nstep[i])*step[i]
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, ...)
    # weights 1 * nt matrix
    w<-intweightsfunc(Nstep[i])
# numerical integration
    res[i] <- crossprod(w , FF)
  }
  res * step 
}

intTD_NC_debug<- function(func=function(x) return(x), intTo, step, Nstep, intweightsfunc = intweights_CAV_SIM,
                     intToStatus=NULL,
                      debug=FALSE, ...){
  # compute numerical integral of func in [0 , intTo] following Newton_Cote method
  # func : function to integrate, func(t, ...)
  # intTo    : upper bound (vector)
  # step : vector of the steps (one row per intTo)
  # Nstep : vector of the number of steps (T = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # intToStatus : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-vector("numeric", length(intTo))
  for(i in 1:length(intTo)){
    # vector of evaluated t
    theT <- (0:Nstep[i])*step[i]
    # vector of the evaluated functions
    FF <- func(theT, i, ...)
    # weights
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i] <- crossprod(w , FF)
  }
  res * step 
}

intTD_SIM3_8 <- function(func=function(x) x, intTo, step, Nstep, ...){
  # compute numerical integral of func in [0 , intTo] following cavalieri Simpson method
  # func : function to integrate
  # intTo    : upper bound (vector)
  # step : vector of the steps (one lig per intTo)
  # Nstep : vector of the number of steps (T = Nstep * step), Nstep = 3 * k
  # weights are (1 3 3 2 3 3 2 3 3  ... 3 3 2 3 3 2 3 3 1 )* step * 3 / 8
  res<-rep(0, length(intTo))
           for(i in 1:length(intTo)){
             res[i] <- 3 * sum(func( (1:(Nstep[i]-1))*step, i, ...)) -
                           sum(func((3*(1:(Nstep[i]/3-1)))*step, i, ...))
                               
           }
           (res + func(0, i, ...) + func(intTo, i, ...) )* step * 3 / 8
}


intTD_BOOLE <- function(func=function(x) x, intTo, step, Nstep, ...){
  # compute numerical integral of func in [0 , intTo] following cavalieri Simpson method
  # func : function to integrate
  # intTo    : upper bound (vector)
  # step : vector of the steps (one lig per intTo)
  # Nstep : vector of the number of steps (T = Nstep * step), Nstep = 4 * k
  # weights are (7 32 12 32 14 32 12 32 14 ... 14 32 12 32 14 32 12 32 7 ) * 4 step / 90 
  res<-rep(0, length(intTo))
           for(i in 1:length(intTo)){
             res[i] <- 32 * sum(func((1:(Nstep[i]/2-1))*step, ...)) + 
                       12 * sum(func((4*(1:(Nstep[i]/4))-2)*step, ...)) +
                       14 * sum(func((4*(1:(Nstep[i]/4-1)))*step, ...))
           }
           (res + 7 * (func(0, ...) + func(intTo, ...)) )* step / 90
}



intTD_GL <- function(func=function(x) x, intTo,
                     step, Nstep,
                     intweightsfunc = NULL,
                     intToStatus=NULL,
                     ...){
  # compute numerical integral of func in [0 , intTo] following Gauss Legendre quadrature
  # func : (vector of) function to integrate, func(t, ...)
  # intTo   : upper bound (vector)
  # Spline : Spline parameters
  # step : points of the quadrature
  # Nstep : weights of the quadrature
  # intweightfunc : unused
  # intToStatus : unused but present for compatibility with inTD_GLM
  res<-vector("numeric", length(intTo))
  # from  0 to intTo then (b-a)/2 = (b+a)/2 = intTo/2
  intTo2 <- intTo/2
  step1 <- step+1
  for(i in 1:length(intTo)){
    res[i] <- sum(Nstep * func(intTo2[i] * step1 , i, ...)) 
  }
  intTo2 * res  
}


intTD_integrate <- function(func=function(x) x, intTo,
                     step, Nstep,
                     intweightsfunc = NULL,
                     intToStatus=NULL,
                     ...){
  # compute numerical integral of func in [0 , intTo] using stats::integrate()
  # func : (vector of) function to integrate, func(t, ...)
  # intTo   : upper bound (vector)
  # Spline : Spline parameters
  # step : points of the quadrature
  # Nstep : weights of the quadrature
  # intweightfunc : unused
  # intToStatus : unused but present for compatibility with inTD_GLM
  # from  0 to intTo then (b-a)/2 = (b+a)/2 = intTo/2
  res<-rep(0, length(intTo))
  for(i in 1:length(intTo)){
    res[i] <-   stats::integrate(func, 0, intTo[i], subdivisions = 100L,
			rel.tol = .Machine$double.eps^0.25, abs.tol = .Machine$double.eps^0.25,
			stop.on.error = TRUE, ...)$value
  }
  res
}





intTD_GLM <- function(func=function(x) return(x), intTo, step, Nstep, intweightsfunc=NULL,
                      intToStatus,
                      debug=FALSE, #Zalphabeta,
                     ...){
  # compute numerical integral of func in [0 , intTo] for equivalence with the poisson GLM trick
  # func : function to integrate, func(t, ...)
  # intTo    : upper bound (vector)
  # step : object of class GLMStepParam
  # Nstep : number of complete bands (< intTo)
  # intTo is in the Nstep'th band
  # intweightfunc function for computing weights : 
  #     weights are (b_i - b_(i-1))
  #                  with b0=0 and bn = intTo[j]
  # intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
  # intToStatus : statuts at intTo
  # ... : parameters of func()
  res<-vector("numeric", length(intTo))
  for(i in 1:length(intTo)){
    # vector of evaluated t
    if(Nstep[i]>0){
      if( intToStatus[i] !=0 ){
        theT <- c(step@points[1:Nstep[i]] , intTo[i])
      }
      else {
        theT <- c(step@points[1:Nstep[i]] , (step@cuts[1+Nstep[i]]+intTo[i])/2)
      }
                                        # vector of the evaluated functions
      FF <- func(theT, i, ...)
                                        # weights
      w<- c(step@steps[1:Nstep[i]], intTo[i]-step@cuts[1+Nstep[i]])

     # numerical integration
      res[i] <- crossprod(w , FF)

    }
    else {
      if( intToStatus[i] != 0 ){
        res[i] <- intTo[i]* func(intTo[i], i,  ...)
      }
      else {
        res[i] <- intTo[i]* func(intTo[i]/2, i, ...)
      }
    }
  }

 res  
}


