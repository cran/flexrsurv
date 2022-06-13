intTD_base_NC <- function(func=function(x) return(x), intTo, Spline,
		step, Nstep, intweightsfunc = intweights_CAV_SIM,
		intToStatus=NULL,
		debug=TRUE,
		...){
	# compute numerical integral of func*base_i(t)  in [0 , intTo] following Newton_Cote method
	# func : (vector of) function to integrate, func(t, ...)
	# intTo    : upper bound (vector)
	# Spline : Spline parameters
	# step : vector of the steps (one row per intTo)
	# Nstep : vector of the number of steps (T = Nstep * step), Nstep is even
	# intweightfunc function for computing weights : 
	#     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
	#     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
	#     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
	# intToStatus : unused but present for compatibility with inTD_GLM
	# ... : parameters of func()
	res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
	for(i in 1:length(intTo)){
		# vector of evaluated t
		theT <- (0:Nstep[i])*step[i]
		# matrix of the evaluated functions (nt row, nfunc col)
		FF <- func(theT, i, ...)
		TBase <- fevaluate(Spline, theT, intercept=TRUE)
		# weights 1 * nt matrix
		w<-intweightsfunc(Nstep[i], step[i])
		
# numerical integration
		res[i,] <- crossprod(w*FF, TBase)
	}
	res  
}

intTD_base_NC_debug<- function(func=function(x) return(x), intTo, Spline,
		step, Nstep, intweightsfunc = intweights_CAV_SIM,
		intToStatus=NULL,
		debug=TRUE,
		...){
	# compute numerical integral of func*base_i(t) in [0 , intTo] following Newton_Cote method
	# func : function to integrate, func(t, ...)
	# intTo    : upper bound (vector)
	# Spline : Spline parameters
	# step : vector of the steps (one row per intTo)
	# Nstep : vector of the number of steps (T = Nstep * step), Nstep is even
	# intweightfunc function for computing weights : 
	#     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
	#     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
	#     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
	# intToStatus : unused but present for compatibility with inTD_GLM
	# ... : parameters of func()
	res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
	for(i in 1:length(intTo)){
		# vector of evaluated t
		theT <- (0:Nstep[i])*step[i]
		TBase <- fevaluate(Spline, theT, intercept=TRUE)
		# vector of the evaluated functions
		FF <- func(theT, i, ...)
		# weights
		w<-intweightsfunc(Nstep[i], step[i])
		
# numerical integration
		res[i,] <- crossprod(w*FF, TBase)
	}
	res  
}



intTD_base_GL <- function(func=function(x) x, intTo, Spline,
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
	res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
	# from  0 to intTo then (b-a)/2 = (b+a)/2 = intTo/2
	intTo2 <- intTo/2
	for(i in 1:length(intTo)){
		# vector of evaluated t
		theT <- intTo2[i] * (step + 1)
		# matrix of the evaluated functions (nt row, nfunc col)
		FF <- func(theT, i, ...)
		TBase <- fevaluate(Spline, theT, intercept=TRUE)
# numerical integration
		res[i,] <- crossprod(Nstep*FF, TBase)
	}
	intTo2 * res
	
}



intTD_base_GLM <- function(func=function(x) return(x), intTo, Spline,
		step, Nstep, intweightsfunc=NULL,
		intToStatus,
		debug=FALSE, ...){
	# compute numerical integral of func*b_i(t) in [0 , intTo] for equivalence with the poisson GLM trick
	# func : function to integrate, func(t, ...)
	# intTo    : upper bound (vector)
	# Spline : Spline parameters
	# step : object of class GLMStepParam
	# Nstep : number of complete bands (< intTo)
	# intTo is in the (Nstep+1)'th band
	# intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
	# intToStatus : statuts at intTo
	#... : parameters of func()
	res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
	for(i in 1:length(intTo)){
		# vector of evaluated t
		if(Nstep[i]>0){
			if( intToStatus[i] ){
				theT <- c(step@points[1:Nstep[i]] , intTo[i])
			}
			else {
				theT <- c(step@points[1:Nstep[i]] , (step@cuts[1+Nstep[i]]+intTo[i])/2)
			}
			# vector of the evaluated functions
			FF <- func(theT, i, ...)
			# matrix of bases evaluated at theT
			TBase <- fevaluate(Spline, theT, intercept=TRUE)
			# weights
			w<- c(step@steps[1:Nstep[i]], intTo[i]-step@cuts[1+Nstep[i]])
			# numerical integration
			res[i,] <- crossprod(w*FF, TBase)
		}
		else {
# Nstep[i]==0
			if( intToStatus[i] ){
				res[i,] <- intTo[i]* func(intTo[i], i, ...)%*% fevaluate(Spline, intTo[i], intercept=TRUE)
			}
			else {
				res[i,] <- intTo[i]* func(intTo[i]/2, i, ...)%*% fevaluate(Spline, intTo[i]/2, intercept=TRUE)
			}
		}
	}
	res  
}

fastintTD_base_GLM <- function(func=function(x) return(x), intTo, Spline,
		step, Nstep, intweightsfunc=NULL,
		intToStatus,
		debug=FALSE,
		...){
	# compute numerical integral of func*b_i(t) in [0 , intTo] for equivalence with the poisson GLM trick
	# func : function to integrate, func(t, ...)
	# intTo    : upper bound (vector)
	# Spline : Spline parameters
	# step : object of class GLMStepParam
	# Nstep : number of complete bands (< intTo)
	# intTo is in the (Nstep+1)'th band
	# intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
	# intToStatus : statuts at intTo
	# ... : parameters of func()
	res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
	# matrix of bases evaluated at the points and intTo
	allTBase <- fevaluate(Spline,step@points , intercept=TRUE)
	Tpoints <- ifelse(intToStatus, intTo,  (step@cuts[1+Nstep]+intTo)/2) 
	TBaseatT <- fevaluate(Spline,Tpoints , intercept=TRUE)
	for(i in 1:length(intTo)){
		# vector of evaluated t
		if(Nstep[i]>0){
			theT <- c(step@points[1:Nstep[i]], Tpoints[i] )
			# vector of the evaluated functions
			FF <- func(theT, i, ...)
			# weights
			w<- c(step@steps[1:Nstep[i]], intTo[i]-step@cuts[1+Nstep[i]] )
			# numerical integration of the complete bands
			res[i,] <- crossprod(w*FF, rbind(allTBase[1:Nstep[i],, drop=FALSE],  TBaseatT[i,]))
			# last bands
		}
		else {
# Nstep[i]==0
			res[i,] <- (intTo[i]-step@cuts[1+Nstep[i]]) * func(Tpoints[i], i, ...) * TBaseatT[i,]
		}
	}
	res  
}

# moins rapide
fastintTD_base_GLM0 <- function(func=function(x) return(x), intTo, Spline,
		step, Nstep, intweightsfunc=NULL,
		intToStatus,
		debug=FALSE,
		...){
	# compute numerical integral of func*b_i(t) in [0 , intTo] for equivalence with the poisson GLM trick
	# func : function to integrate, func(t, ...)
	# intTo    : upper bound (vector)
	# Spline : Spline parameters
	# step : object of class GLMStepParam
	# Nstep : number of complete bands (< intTo)
	# intTo is in the (Nstep+1)'th band
	# intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
	# intToStatus : statuts at intTo
	# ... : parameters of func()
	res<-matrix(0, nrow = length(intTo), ncol = Spline@nbases + Spline@log)
	# matrix of bases evaluated at the points and intTo
	allTBase <- fevaluate(Spline,step@points , intercept=TRUE)
	Tpoints <- ifelse(intToStatus, intTo,  (step@cuts[1+Nstep]+intTo)/2) 
	TBaseatT <- fevaluate(Spline,Tpoints , intercept=TRUE)
	for(i in 1:length(intTo)){
		# vector of evaluated t
		if(Nstep[i]>0){
			theT <- c(step@points[1:Nstep[i]] )
			# vector of the evaluated functions
			FF <- func(theT, i, ...)
			# weights
			w<- c(step@steps[1:Nstep[i]] )
			# numerical integration of the complete bands
			res[i,] <- crossprod(w*FF, allTBase[1:Nstep[i],, drop=FALSE])
			# last bands
			res[i,] <- res[i,] +
					(intTo[i]-step@cuts[1+Nstep[i]])*func(Tpoints[i], i, ...) *  TBaseatT[i,]
		}
		else {
# Nstep[i]==0
			res[i,] <- (intTo[i]-step@cuts[1+Nstep[i]]) * func(Tpoints[i], i, ...) * TBaseatT[i,]
		}
	}
	res  
}















