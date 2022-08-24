intTDft_NC <- function(func=function(x) return(x), intFrom, intTo,
		step, Nstep, intweightsfunc = intweights_CAV_SIM,
		intToStatus=NULL,
		debug=FALSE,
		...){
	# compute numerical integral of func in [intFrom , intTo] following Newton_Cote method
	# func : (vector of) function to integrate, func(t, ...)
	# intFrom : lower bound (vector)
	# intTo   : upper bound (vector)
	# step : vector of the steps (one row per intTo)
	# Nstep : vector of the number of steps ((intTo - intFrom) = Nstep * step), Nstep is even
	# intweightfunc function for computing weights : 
	#     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
	#     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
	#     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
	# intToStatus : unused but present for compatibility with inTD_GLM
	# ... : parameters of func()
	res<-vector("numeric", length(intTo))
	for(i in 1:length(intTo)){
		# vector of evaluated t
		theT <- intFrom[i] + (0:Nstep[i])*step[i]
		# matrix of the evaluated functions (nt row, nfunc col)
		FF <- func(theT, i, ...)
		# weights 1 * nt matrix
		w<-intweightsfunc(Nstep[i],step[i])
		
# numerical integration
		res[i] <- crossprod(w , FF)
	}
	res 
}

intTDft_NC0 <- function(func=function(x) return(x), intFrom, intTo,
		step, Nstep, degree = 4L, intweightsfunc = intweights_CAV_SIM,
		intToStatus=NULL,
		debug=FALSE,
		...){
	# compute numerical integral of func in [intFrom , intTo] following Newton_Cote method
	# func : (vector of) function to integrate, func(t, ...)
	# intFrom : lower bound (vector)
	# intTo   : upper bound (vector)
	# step : vector of the steps (one row per intTo)
	# Nstep : vector of the number of steps ((intTo - intFrom) = Nstep * step), Nstep is even
	# intweightfunc function for computing weights : 
	#     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is 2 (even
	#     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
	#     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
	# intToStatus : unused but present for compatibility with inTD_GLM
	# ... : parameters of func()
	func <- match.fun(func)
	ff <- function(x, i) func(x, i, ...)
#    res <- .External(C_call_intTDft_NC, ff, rho = environment(), 
#                     as.double(intFrom), as.double(intTo),
#                     as.double(step), as.integer(Nstep),
#                     as.integer(intweightsfunc), 
#                     as.integer(debug))
	res <- .Call(C_intTDft_NC, ff, 
			as.double(intFrom), as.double(intTo),
			as.double(step), as.integer(Nstep), as.integer(max(Nstep)),
			as.integer(degree), environment())
	res 
}

intTDft_NC_debug<- function(func=function(x) return(x),  intFrom, intTo,
		step, Nstep, intweightsfunc = intweights_CAV_SIM,
		intToStatus=NULL,
		debug=FALSE, ...){
	# compute numerical integral of func in [intFrom , intTo] following Newton_Cote method
	# func : function to integrate, func(t, ...)
	# intFrom : lower bound (vector)
	# intTo   : upper bound (vector)
	# step : vector of the steps (one row per T)
	# Nstep : vector of the number of steps ((intTo - intFrom) = Nstep * step), Nstep is even
	# intweightfunc function for computing weights : 
	#     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
	#     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
	#     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
	# intToStatus : unused but present for compatibility with inTD_GLM
	# ... : parameters of func()
	res<-vector("numeric", length(intTo))
	for(i in 1:length(intTo)){
		# vector of evaluated t
		theT <- intFrom[i] + (0:Nstep[i])*step[i]
		# vector of the evaluated functions
		FF <- func(theT, i, ...)
		# weights
		w<-intweightsfunc(Nstep[i], step[i])
		
# numerical integration
		res[i] <- crossprod(w , FF)
	}
	res 
}

intTDft_SIM3_8 <- function(func=function(x) x, intFrom, intTo,
		step, Nstep, ...){
	# compute numerical integral of func in [intFrom , intTo] following cavalieri Simpson method
	# func : function to integrate
	# intFrom : lower bound (vector)
	# intTo   : upper bound (vector)
	# step : vector of the steps (one lig per T)
	# Nstep : vector of the number of steps ((intTo - intFrom) = Nstep * step), Nstep = 3 * k
	# weights are (1 3 3 2 3 3 2 3 3  ... 3 3 2 3 3 2 3 3 1 )* step * 3 / 8
	res<-rep(0, length(intTo))
	for(i in 1:length(intTo)){
		res[i] <- 3 * sum(func( intFrom + (1:(Nstep[i]-1))*step, i, ...)) -
				sum(func(intFrom + (3*(1:(Nstep[i]/3-1)))*step, i, ...))
		
	}
	(res + func(intFrom, i, ...) + func(intTo, i, ...) )* step * 3 / 8
}


intTDft_BOOLE <- function(func=function(x) x, intFrom, intTo,
		step, Nstep, ...){
	# compute numerical integral of func in [intFrom , intTo] following cavalieri Simpson method
	# func : function to integrate
	# intFrom : lower bound (vector)
	# intTo   : upper bound (vector)
	# step : vector of the steps (one ligne per T)
	# Nstep : vector of the number of steps ((intTo - intFrom) = Nstep * step), Nstep = 4 * k
	# weights are (7 32 12 32 14 32 12 32 14 ... 14 32 12 32 14 32 12 32 7 ) * 4 step / 90 
	res<-vector("numeric", length(intTo))
	for(i in 1:length(intTo)){
		res[i] <- 32 * sum(func(intFrom + (1:(Nstep[i]/2-1))*step, ...)) + 
				12 * sum(func(intFrom + (4*(1:(Nstep[i]/4))-2)*step, ...)) +
				14 * sum(func(intFrom + (4*(1:(Nstep[i]/4-1)))*step, ...))
	}
	(res + 7 * (func(intFrom, ...) + func(intTo, ...)) )* step / 90
}

intTDft_GL <- function(func=function(x) x, intFrom, intTo,
		step, Nstep,
		intweightsfunc = NULL,
		intToStatus=NULL,
		...){
	# compute numerical integral of func in [intFrom , intTo] following Gauss Legendre quadrature
	# func : (vector of) function to integrate, func(t, ...)
	# intFrom : lower bound (vector)
	# intTo   : upper bound (vector)
	# Spline : Spline parameters
	# step : points of the quadrature
	# Nstep : weights of the quadrature
	# intweightfunc : unused
	# intToStatus : unused but present for compatibility with inTD_GLM
	# evaluation points are (b-a)/2 * step + (b+a)/2 =  dT * step + Tmid
	res<-vector("numeric", length(intTo))
	Tmid <- (intTo + intFrom)/2
	dT   <- (intTo - intFrom)/2
	for(i in 1:length(intTo)){
		res[i] <- sum(Nstep * func(dT[i] * step + Tmid[i], i, ...)) 
	}
	dT * res
}




intTDft_GLM <- function(func=function(x) return(x), intFrom, intTo,
		step,
		Nstep, 
		intweightsfunc=NULL,
		intToStatus,
		debug=FALSE, #Zalphabeta,
		...){
	# compute numerical integral of func in [0 , T] for equivalence with the poisson GLM trick
	# func : function to integrate, func(t, ...)
	# intFrom : lower bound (vector)
	# intTo   : upper bound (vector)
	# step : object of class GLMStepParam
	# Nstep : index of the first and last complete band ( intFrom[i] < step@cuts[Nstep[i,1]] <= step@cuts[Nstep[i,2]+1] < intTo)
	#                                                    ( intFrom[i] < step@points[Nstep[i,1]] <= step@points[Nstep[i,2]] < intTo)
	# Nstep  : index of the first complete band 
	# intTo is in the Nstep'th band
	# intFrom is in the first band
	# intweightfunc function for computing weights : 
	#     weights are (b_i - b_(i-1))
	#                  with b0=0 and bn = T[j]
	# intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
	# intToStatus : statuts at intTo
	# ... : parameters of func()
	res<-vector("numeric", length(intTo))
	for(i in 1:length(intTo)){
		# vector of evaluated t
		if(Nstep[i,2]>= Nstep[i,1]){
			# at least one complete step
			if( intToStatus[i] ){
				theT <- c((step@cuts[Nstep[i,1]]+intFrom[i])/2,
						step@points[Nstep[i,1]:Nstep[i,2]] ,
						intTo[i])
			}
			else {
				theT <- c((step@cuts[Nstep[i,1]]+intFrom[i])/2,
						step@points[Nstep[i,1]:Nstep[i,2]] ,
						(step@cuts[1+Nstep[i,2]]+intTo[i])/2)
			}
			# vector of the evaluated functions
			FF <- func(theT, i, ...)
			# weights
			w<- c(step@cuts[Nstep[i,1]] - intFrom[i],
					step@steps[Nstep[i,1]:Nstep[i,2]],
					intTo[i]-step@cuts[1+Nstep[i,2]])
			
			# numerical integration
			res[i] <- crossprod(w , FF)
		}
		else if((Nstep[i,2] - Nstep[i,1]) == -1L){
# intFrom and intTo are in 2 successive bands
			# Nstep[i,2] + 1 =  Nstep[i,1]
			if( intToStatus[i] ){
				theT <- c((step@cuts[Nstep[i,1]]+intFrom[i])/2,
						intTo[i])
			}
			else {
				theT <- c((step@cuts[Nstep[i,1]]+intFrom[i])/2,
						(step@cuts[Nstep[i,1]]+intTo[i])/2)
			}
			# vector of the evaluated functions
			FF <- func(theT, i, ...)
			# weights
			w<- c(step@cuts[Nstep[i,1]] - intFrom[i],
					intTo[i]-step@cuts[Nstep[i,1]])
			
			# numerical integration
			res[i] <- crossprod(w , FF)
		}
		else {   #if((Nstep[i,2] - Nstep[i,1]) == -2L){
# intFrom and intTo are in the same band
			if( intToStatus[i] ){
				res[i] <- (intTo[i] - intFrom[i]) * func(intTo[i], i,  ...)
			}
			else {
				res[i] <- (intTo[i] - intFrom[i]) * func((intTo[i] + intFrom[i])/2, i,  ...)
			}
		}
	}
	
	
	res  
}




