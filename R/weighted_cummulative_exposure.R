################
#    WCE at theT

weighted_cummulative_exposure <- function(Increment,
		fromT, toT, FirstId, LastId,
		theT, tId,
		eta0=NULL, iWbeg, iWend,
		ISpline_W, Intercept_W){
	
	# computes the weighted cummulative exposures at theT for subject tId
	#  it is the cumsum of the contributions of all the lines by individuals
	# the return value has the sqame length as fromT
	# the actual WCE at final T of each indivudual is the last value of the individual 
	#   WCE_id(theT) = sum_ID ( Icrement ISpline_W(theT - fromT) )
	#
	# tId : rows with same Id are from the same individual
	# tId : index such that each thet is in ] fromT[tId], toT[tId]]
	#  it is assumed that for each Id, the profile is such that the last interval is right-open [fromT[last] + infty[
	#             (ie that t <= toT[last])  
	# Increment :  matrix of exposure INCREMENT variables 
	#  fromT, toT, beginning, end  at the end of ]fromT, toT] for the row
	# eta0 : vector all the coef of WCE, if NULL, ISpline_W are assumed to be scaled by the coef of the weghting function
	# iWbeg, iWend : coef of the ith WCE variable is eta0[iWbeg[i]:iWend[i]]
	# ISpline_W : list of the integrated splines parameters for the WCE effects
	IS_W <- ISpline_W
	if(!is.null(eta0)){
		for(iW in 1:length(ISpline_W)){
			if(Intercept_W[[iW]]){
				IS_W[[iW]] <- ISpline_W[[iW]] * eta0[iWbeg[iW]:iWend[iW]]
			}
			else {
				IS_W[[iW]]<- ISpline_W[[iW]] * c(0, eta0[iWbeg[iW]:iWend[iW]])
			}
		}
	}
	
#    wce <- matrix(0, nrow=dim(Increment)[1], ncol = dim(Increment)[2])
	WCE <- matrix(0, nrow=dim(Increment)[1], ncol = length(ISpline_W))
	for(iW in 1:length(ISpline_W)){
		WCE[,iW] <- predictwce(object=IS_W[[iW]], t=theT, Increment= Increment[,iW], fromT=fromT, tId=tId,
				FirstId=FirstId, LastId=LastId, intercept=Intercept_W[iW], outer.ok=TRUE)
	}
	
	return(WCE)
}


