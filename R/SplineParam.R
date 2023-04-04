# Methods for classes *SplineBasis
# classes *SplineBasis are defined in AllClass.R 


###################################################################################
#######          createurs
###################################################################################

# knots are the full set of knots used to define the basis functions.
BSplineBasis0<-function(knots, degree=3, keep.duplicates=FALSE, log=FALSE, clog=ifelse(log,1.0, 0.0) ) {
	# for M-splines,
	#
	# the bases are not defined outside of [kmin, kMax]
	#
	# knots are 1 boundary min knots, interior knots, 1 boundary max knots
	# boundary knots are not tested to be duplicated
	# code from orthogonalsplinebasis
	order<-degree+1
	n<-length(knots)
	if (any(table(knots[order:(n - order + 1)]) > 1) && !keep.duplicates) {
		warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
		knots <- unique(knots[order:(n - order + 1)])
		knots <- knots[c(rep(1, order-1), seq(length(knots)), rep(length(knots), order-1))]
	}
	# recompute n number of knots
	n <- length(knots)
	q <- n-order
	
	SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=order, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	new("BSplineBasis", knots=knots, min=knots[1], max=knots[length(knots)],
			degree=as.integer(degree), nbases=as.integer(q), Matrices=M, SplineBasis=SB, log=log, clog=clog)
}

BSplineBasis<-function(knots, degree=3, keep.duplicates=FALSE, log=FALSE, clog=ifelse(log,1.0, 0.0)) {
	# for M-splines,
	#
	# the bases are not defined outside of [kmin, kMax]
	#
	# knots are 1 boundary min knots, interior knots, 1 boundary max knots
	# boundary knots are not tested to be duplicated
	# code from orthogonalsplinebasis
	order<-degree+1
	n<-length(knots)
	if( n ==2 ){
# no interior knots
		knots<-c(rep(knots[1], order), rep(knots[n], order))
	}
	else {
		if(any(table(knots[2:(n-1)])>1)&&!keep.duplicates){
			warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
			# modif MGk 06/06/2011 to keep the (order-1) first and last knots 
			oknots <- knots
			iknots<-unique(oknots[2:(n-1)])
			knots<-c(rep(oknots[1], order), iknots, rep(oknots[n], order))
		}
		else {
			# just duplicates first and last knots
			oknots <- knots
			iknots<-oknots[2:(n-1)]
			knots<-c(rep(oknots[1], order), iknots, rep(oknots[n], order))
		}            
	}
	# recompute n number of knots
	n <- length(knots)
	q <- n-order
	
	SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=order, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	new("BSplineBasis", knots=knots, min=knots[1], max=knots[length(knots)],
			degree=as.integer(degree), nbases=as.integer(q), Matrices=M, SplineBasis=SB, log=log, clog=clog)
}



MSplineBasis<-function(knots, degree=3, keep.duplicates=FALSE, log=FALSE, clog=ifelse(log,1.0, 0.0) ) {
	# for M-splines basis int_boundinf^^boudsup b_i(t)=1
	# just scaled BSplineBasis
	# knots are interior knots; uplicated knots are allowed for discontinuity conditions
	# boundary knots are not tested to be duplicated
	# code using thogonalsplinebasis
	tmpobj <- BSplineBasis(knots=knots, degree=degree, keep.duplicates=keep.duplicates, log=FALSE) 
	xscale<-evaluate(integrate(tmpobj), max(knots), intercept=TRUE)
	tmpobj<- tmpobj * (1/as.numeric(xscale))
	new("MSplineBasis", knots=tmpobj@knots, min=tmpobj@min, max=tmpobj@max,
			degree=tmpobj@degree, nbases=tmpobj@nbases, Matrices=tmpobj@Matrices, SplineBasis=tmpobj@SplineBasis, log=log, clog=clog)
}

maketpdegrees <- function(knots, order){
	order - unlist(lapply(table(knots), function(x) x:1))
}



EBSplineBasis<-function(knots, degree=3L, keep.duplicates=FALSE) {
	# for 0-extrapolated B-splines,
	# b_i(x) is assuemed to be 0 before min(knots) and after max(knots)
	# no regularity conditions at boundary knots. (excepetd that
	# if the coef of the first and the last basis are 0 then the spline are continuous
	# at the boundary knots with b_i(kmin)=b_u(kmax) = 0
	#
	# when integrating (and derivating), the 0-extrapolation is integrated 
	#
	# knots are 1 boundary min knots, interior knots, 1 boundary max knots
	# boundary knots are not tested to be duplicated
	# code from orthogonalsplinebasis
	BS<-BSplineBasis(knots=knots, degree=degree, keep.duplicates=keep.duplicates)
	
	dims<- dim(BS@Matrices)
	dims2<-dims
	dims2[3]<-dims[3]+1
	
# ori	M2<-array(        dim=c(order,q,n-2*order+1))
# in init of extended Spline bases, the matrix correspnding to the interval [ max, +inf[ is 0
# 
	M2<-array(data=0, dim=dims2)
	for(i in seq(1, dims[3])) {  #Identifying interior intervals
		M2[,,i]<-BS@Matrices[,,i]
	}
	new("EBSplineBasis", knots=BS@knots, min=BS@min, max=BS@max,
			degree=BS@degree, nbases=BS@nbases, Matrices=M2)
}


EMSplineBasis<-function(knots, degree=3L, keep.duplicates=FALSE) {
	# for 0-extrapolated M-splines,
	# b_i(x) is assuemed to be 0 before min(knots) and after max(knots)
	# no regularity conditions at boundary knots. (excepetd that
	# if the coef of the first and the last basis are 0 then the spline are continuous
	# at the boundary knots with b_i(kmin)=b_u(kmax) = 0
	#
	# when integrating (and derivating), the 0-extrapolation is integrated 
	#
	# knots are 1 boundary min knots, interior knots, 1 boundary max knots
	# boundary knots are not tested to be duplicated
	# code from orthogonalsplinebasis
	MS<-MSplineBasis(knots=knots, degree=degree, keep.duplicates=keep.duplicates)
	
	dims<- dim(MS@Matrices)
	dims2<-dims
	dims2[3]<-dims[3]+1
	
# ori	M2<-array(        dim=c(order,q,n-2*order+1))
# in init of extended Spline bases, the matrix correspnding to the interval [ max, +inf[ is 0
# 
	M2<-array(data=0, dim=dims2)
	for(i in seq(1, dims[3])) {  #Identifying interior intervals
		M2[,,i]<-MS@Matrices[,,i]
	}
	new("EMSplineBasis", knots=MS@knots, min=MS@min, max=MS@max,
			degree=MS@degree, nbases=MS@nbases, Matrices=M2)
}


LEBSplineBasis<-function(knots, degree=3, keep.duplicates=FALSE, log=FALSE) {
	# for linearly extended B-splines,
	# knots are 1 boundary min knots, interior knots, 1 boundary max knots
	# boundary knots are not tested to be duplicated
	# code from orthogonalsplinebasis
	BS<-BSplineBasis(knots=knots, degree=degree, keep.duplicates=keep.duplicates)
	dims<- dim(BS@Matrices)
	# coef of the linear extrapolation
	linexinf <- matrix(0, ncol=dims[2], nrow=dims[1])
	linexinf[1,] <- evaluate(BS, knots[1], intercept=TRUE)
	linexinf[2,] <- evaluate(deriv(BS), knots[1], intercept=TRUE)
	
	linexsup <-  matrix(0, ncol=dims[2], nrow=dims[1])
	linexsup[1,] <- evaluate(BS, BS@max, intercept=TRUE)
	linexsup[2,] <- evaluate(deriv(BS), BS@max, intercept=TRUE)
	
	new("LEBSplineBasis", knots=BS@knots, min=BS@min, max=BS@max,
			degree=BS@degree, nbases=BS@nbases, linexinf=linexinf, linexsup=linexsup,
			orderextrapol = 1L, 
			Matrices=BS@Matrices, SplineBasis=BS@SplineBasis, log=FALSE)
}


LEBSplineBasis0<-function(knots, degree=3, keep.duplicates=FALSE, log=FALSE) {
	# for linearly extended M-splines,
	# knots are 1 boundary min knots, interior knots, 1 boundary max knots
	# boundary knots are not tested to be duplicated
	# code from orthogonalsplinebasis
	order<-degree+1
	n<-length(knots)
	if( n ==2 ){
# no interior knots
		knots<-c(rep(knots[1], order), rep(knots[n], order))
	}
	else {
		if(any(table(knots[2:(n-1)])>1)&&!keep.duplicates){
			warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
			# modif MGk 06/06/2011 to keep the (order-1) first and last knots 
			oknots <- knots
			iknots<-unique(oknots[2:(n-1)])
			knots<-c(rep(oknots[1], order), iknots, rep(oknots[n], order))
		}
		else {
			# just duplicates first and last knots
			oknots <- knots
			iknots<-oknots[2:(n-1)]
			knots<-c(rep(oknots[1], order), iknots, rep(oknots[n], order))
		}            
	}
	# recompute n number of knots
	n <- length(knots)
	q <- n-order
	
	SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=order, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	
	# coef of the linear extrapolation
	linexinf <- matrix(0, ncol=q, nrow=order)
	linexinf[1,] <- evaluate(SB, knots[1])
	linexinf[2,] <- evaluate(deriv(SB), knots[1])
	
	linexsup <-  matrix(0, ncol=q, nrow=order)
	linexsup[1,] <- evaluate(SB, knots[n])
	linexsup[2,] <- evaluate(deriv(SB), knots[n])
	
	new("LEBSplineBasis", knots=knots, min=knots[1], max=knots[length(knots)], linexinf=linexinf, linexsup=linexsup,
			orderextrapol = 1L, 
			degree=as.integer(degree), nbases=as.integer(q), Matrices=M, SplineBasis=SB, log=log)
}


# Restricted Bsplinebasis : Linear-Tailed Bspline with C2 derivatbility at boundary-knots
#                            linear extrapolation + derivative(of order 2 at boundaries) == 0 and derivative o order <=2  continuity condition at boundary knots
# if order = 4, restricted cubic spline
# R2BSplinBasis are LEBSplineBasis with specific Matrices and number of basis
# 		linear bases are first and last bases components if firstlinear
#		B2[Kmin] = B_(N-1)[Kmin] = 1
# same as LTBSplineBasis but constrained on order 2 derivative continuity
R2bBSplineBasis<-function(knots,
		degree=3,
		keep.duplicates=FALSE,
		log=FALSE,
		firstlinear=TRUE) {
	# for restricted B-splines (linear extrapolation + 2nd derivative at boundaries == 0 )
	# when degree > 2 , 2 bases less than the corresponding BSplineBasis  
	# knots are 1 boundary min knots, interior knots, 1 boundary max knots
	# dim(R2bB) length(knots) + degre - 3
	# boundary knots are not tested to be duplicated
	# scaled for  B1(Kmin)' = B1(Kmin)' =1
	
	# build from linearly extended BSplinebasis
	# inner bases compenents are left unchanged
	# the degree left  basis components are replaced by 2 basis components
	# the degree right basis components are replaced by 2 basis components
	
	
# if firstlinear=TRUE, the first basis is the linear basis, the second is the constant basis (the reverse at the right bounday)
# if firstlinear=FALSE,the second basis is the linear basis, the first is the constant basis (the reverse at the right bounday)
	
	Nk <- length(knots)
	Nb <- Nk + degree -1 
	if(degree < 1){
		stop("unable to build R2bBSplineBasis with degree 0")
	}
	if(Nb < 3){
		stop(paste0("unable to build R2bBSplineBasis with ", Nk, " knots and degree ", degree))
	}
# build the complete LEBSplineBasis 
	tmpobj <- LEBSplineBasis(knots=knots,
			degree=degree,
			keep.duplicates=keep.duplicates,
			log=FALSE) 
	
	if(degree <= 2){
		# LEBS is a R2bB, but we need to change basis, # newNb = Nb
		newNb <- Nb
		
		derivboundaries <- diag((evaluate(deriv(tmpobj), knots[c(1, Nk)]))[, c(1, Nb)])
		if(Nb>3){
			# at left
			cleft <- cbind(
					c(1/derivboundaries[1], 0), 
					c(1, 1)
			)
			
			# at right
			cright <- cbind(
					c(1, 1),
					c(0, 1/derivboundaries[2])
			)
			
			# change basis matrix	
			mm <- diag(Nb)	
			mm[1:2, 1:2] <- cleft
			mm[(1:2)+(Nk+degree-1-2), 1:2 + Nk+degree-1-2] <- cright
		}
		else {
			# Nb=Nk=3
			mm <- diag(3)
			mm[1,1]<- 1/derivboundaries[1]
			mm[3,3]<- 1/derivboundaries[2]
			mm[1,2]<- 1
			mm[3,2]<- 1		
		}
	}
	else {
		# number of "interior" bases (appart from degree components at left and degree components at right)  
		# the interior nbinside basis component are left unchanged
		nbinside <- Nb - 2 * 3 
		if( nbinside < -1 ){
			stop("unable to build R2bBSplineBasis")
		} else {
			# newNb = Nb + degree - 1 -2 = Nk + degree -3
			newNb <- Nk + degree -3 
			
			# if nbinside >= 0
			if(nbinside > -1){
				marsden <- MarsdenIdentity(tmpobj)
				
				# at left
				
				cleft <- cbind(
						c(marsden[1:3]-marsden[3]),
						rep(1, 3)
				)			
				
				
				# at right
				
				cright <- cbind(
						rep(1, 3),
						c(marsden[1:3+Nb-3]-marsden[1+Nb-3])
				)
				
				
				# change basis matrix	
				
				mm <- diag(Nb)[, 2:(Nb-1)]				
				mm[1:3, 1:2] <- cleft
				mm[(1:3)+Nb-3, 1:2 + Nb-2-2] <- cright
				
			}
			if(nbinside == -1){
				# nbinside == -1
				# Nk + degre = 2*3 
				# Nb = 5
				# newNb = 3 : left linear asymptopte, right linear asymptot, constant basis componenet 
				# the constant basis component is the um of the B-basis components
				
				# using Marsden's identity
				marsden <- MarsdenIdentity(tmpobj)
				
				mm <- diag(Nb)[, (3-1):(Nb-1)]
				mm[1:3, 1] <- c(marsden[1:3]-marsden[3])
				mm[,2] <- rep(1, Nb)
				mm[(1:3)+Nb-3, 3] <- c(marsden[1:3+Nb-3]-marsden[1+Nb-3])
				
			}
		}
	}
	
	change_basis <- t(mm)
	newobj <- change_basis %*% tmpobj
	
	if(!firstlinear){
		mperm <- diag(newNb)
		mperm[1,1]<- 0
		mperm[2,2]<- 0
		mperm[1,2]<- 1
		mperm[2,1]<- 1
		if(newNb > 3){
			mperm[newNb  ,newNb]  <- 0
			mperm[newNb-1,newNb-1]<- 0
			mperm[newNb-1,newNb]  <- 1
			mperm[newNb  ,newNb-1]<- 1
			
		}
		newobj <- mperm %*% newobj
		
	}		
	
	class(newobj) <- "R2bBSplineBasis"
	newobj
}


# Linear-Tailed Bsplinebasis : linear extrapolation + derivative(of order > 1 at boundaries) == 0 and derivative continuity condition at boundary knots
# if order = 4, restricted cubic spline
# LTBSplinBasis are LEBSplineBasis with specific Matrices and number of basis
# 		linear bases are first and last bases components if firstlinear
#		B2[Kmin] = B_(N-1)[Kmin] = 1
LTBSplineBasis<-function(knots,
		degree=3,
		keep.duplicates=FALSE,
		log=FALSE,
		firstlinear=TRUE) {
	# for restricted B-splines (linear extrapolation + 2nd derivative at boundaries == 0 + derivatives continuity condition at boundary knots)
	# 2*(degree - 2) bases less than the corresponding BSplineBasis 
	# knots are 1 boundary min knots, interior knots, 1 boundary max knots
	# dim(LTB) length(knots) - degre + 3
	# boundary knots are not tested to be duplicated
	# scaled for  B1(Kmin)' = B1(Kmin)' =1
	
	# build from linearly extended BSplinebasis
	# inner bases compenents are left unchanged
	# the degree left  basis components are replaced by 2 basis components
	# the degree right basis components are replaced by 2 basis components
	
	
# if firstlinear=TRUE, the first basis is the linear basis, the second is the constant basis (the reverse at the right bounday)
# if firstlinear=FALSE,the second basis is the linear basis, the first is the constant basis (the reverse at the right bounday)
	
	Nk <- length(knots)
	Nb <- Nk + degree -1 
	if(degree < 1){
		stop("unable to build LTBSplineBasis with degree 0")
	}
	if(Nb < 3){
		stop(paste0("unable to build LTBSplineBasis with ", Nk, " knots and degree ", degree))
	}
# build the complete LEBSplineBasis 
	tmpobj <- LEBSplineBasis(knots=knots,
			degree=degree,
			keep.duplicates=keep.duplicates,
			log=FALSE) 
	
	if(degree <= 2){
		# LEBS is a LTBS, but we need to change basis, # newNb = Nb
		newNb <- Nb
		
		derivboundaries <- diag((evaluate(deriv(tmpobj), knots[c(1, Nk)]))[, c(1, Nb)])
		if(Nb>3){
			# at left
			cleft <- cbind(
					c(1/derivboundaries[1], 0), 
					c(1, 1)
			)
			
			# at right
			cright <- cbind(
					c(1, 1),
					c(0, 1/derivboundaries[2])
			)
			
			# change basis matrix	
			mm <- diag(Nb)	
			mm[1:2, 1:2] <- cleft
			mm[(1:2)+(Nk+degree-1-2), 1:2 + Nk+degree-1-2] <- cright
		}
		else {
			# Nb=Nk=3
			mm <- diag(3)
			mm[1,1]<- 1/derivboundaries[1]
			mm[3,3]<- 1/derivboundaries[2]
			mm[1,2]<- 1
			mm[3,2]<- 1		
		}
	}
	else {
		# number of "interior" bases (appart from degree components at left and degree components at right)  
		# the interior nbinside basis component are left unchanged
		nbinside <- Nb - 2 * (degree)
		if( nbinside < -1 ){
			stop("unable to build LTBSplineBasis")
		} else {
			# newNb = Nb - 2*(degree -2) = Nk - degree + 3
			newNb <- Nk - degree + 3
			
			# if nbinside > 0
			if(Nk > degree){
				marsden <- MarsdenIdentity(tmpobj)
				
				# at left
				
				cleft <- cbind(
						c(marsden[1:degree]-marsden[degree]),
						rep(1, degree)
				)			
				
				
				# at right
				
				cright <- cbind(
						rep(1, degree),
						c(marsden[1:degree+Nb-degree]-marsden[1+Nb-degree])
				)
				
				
				# change basis matrix	
				
				mm <- diag(Nb)[, (degree-1):(Nk+1)]				
				mm[1:degree, 1:2] <- cleft
				mm[(1:degree)+(Nk-1), 1:2 + Nk-(degree-2)-1] <- cright
				
			}
			if(Nk == degree){
				# nbinside == -1
				# Nb = 2 Nk -1, 
				# newNb = 3 : left linear asymptopte, right linear asymptot, constant basis componenet 
				# the constant basis component is the um of the B-basis components
				
				# using Marsden's identity
				marsden <- MarsdenIdentity(tmpobj)
				
				mm <- diag(Nb)[, (degree-1):(Nk+1)]
				mm[1:degree, 1] <- c(marsden[1:degree]-marsden[degree])
				mm[,2] <- rep(1, Nb)
				mm[(1:degree)+(Nk-1), 3] <- c(marsden[1:degree+Nb-degree]-marsden[1+Nb-degree])
				
			}
		}
	}
	
	change_basis <- t(mm)
	newobj <- change_basis %*% tmpobj
	
	if(!firstlinear){
		mperm <- diag(newNb)
		mperm[1,1]<- 0
		mperm[2,2]<- 0
		mperm[1,2]<- 1
		mperm[2,1]<- 1
		if(newNb > 3){
			mperm[newNb  ,newNb]  <- 0
			mperm[newNb-1,newNb-1]<- 0
			mperm[newNb-1,newNb]  <- 1
			mperm[newNb  ,newNb-1]<- 1
			
		}
		newobj <- mperm %*% newobj
		
	}		
	
	class(newobj) <- "LTBSplineBasis"
	newobj
}



BSplineBasis2<-function(allknots, degree=3, keep.duplicates=FALSE, log=FALSE) {
	# for B-splines, allknots are degree+1 boundary min knots, interior knots, degere+1 boundary max knots
	# the 2x4 boundary knots are equal
	# code from orthogonalsplinebasis
	order<-degree+1
	n<-length(allknots)
	knots<-allknots  # !! all the boundary knots are given
	
	if( n ==2*order ){
# no interior knots
		knots<-allknots  # !! all the boundary knots are given
	}
	else {
		if(any(table(allknots[order:(n-order+1)])>1)&&!keep.duplicates){
			warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
			# modif MGk 06/06/2011 to keep the (order-1) first and last knots 
			oknots <- allknots
			iknots<-unique(oknots[(order+1):(n-order)])
			knots<-c(rep(oknots[order], order), iknots, rep(oknots[n-order+1], order))
		}
	}
	q<-n-order
	SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=order, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	new("BSplineBasis", knots=knots, min=knots[1], max=knots[length(knots)],
			degree=as.integer(degree), nbases=as.integer(q), Matrices=M, SplineBasis=SB, log=log)
}

LEBSplineBasis2<-function(allknots, degree=3, keep.duplicates=FALSE, log=FALSE) {
	# for linearly extended M-splines,
	# for B-splines, allknots are degree+1 boundary min knots, interior knots, degere+1 boundary max knots
	# the 2x4 boundary knots are equal
	# code from orthogonalsplinebasis
	order<-degree+1
	n<-length(allknots)
	knots<-allknots  # !! all the boundary knots are given
	if( n ==2*order ){
# no interior knots
		knots<-allknots  # !! all the boundary knots are given
	}
	else {
		if(any(table(allknots[order:(n-order+1)])>1)&&!keep.duplicates){
			warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
			# modif MGk 06/06/2011 to keep the (order-1) first and last knots 
			oknots <- allknots
			iknots<-unique(oknots[(order+1):(n-order)])
			knots<-c(rep(oknots[order], order), iknots, rep(oknots[n-order+1], order))
		}
	}
	q<-n-order
	SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=order, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	
	# coef of the linear extrapolation
	linexinf <- matrix(0, ncol=q, nrow=order)
	linexinf[1,] <- evaluate(SB, knots[1])
	linexinf[2,] <- evaluate(deriv(SB), knots[1])
	
	linexsup <-  matrix(0, ncol=q, nrow=order)
	linexsup[1,] <- evaluate(SB, knots[n])
	linexsup[2,] <- evaluate(deriv(SB), knots[n])
	
	
	new("LEBSplineBasis", knots=knots, min=knots[1], max=knots[length(knots)], linexinf=linexinf, linexsup=linexsup,
			degree=as.integer(degree), nbases=as.integer(q), Matrices=M, SplineBasis=SB, log=log)
}

BSplineBasis2<-function(knots, degree=3,  Boundary.knots , keep.duplicates=FALSE, log=FALSE ) {
	# for B-splines mmimic of splines:bs()
	# knots are interior knots; uplicated knots are allowed for discontinuity conditions
	# boundary knots are not tested to be duplicated
	# code using thogonalsplinebasis
	
	
	order<-degree+1
	
	n<-length(knots)
	if( n ==0 ){
# no interior knots
	}
	else {
		if(any(table(knots)>1)&&!keep.duplicates){
			warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
			# modif MGk 06/06/2011 to keep the (order-1) first and last knots 
			oknots <- knots
			knots<-unique(oknots[order:(n-order+1)])
		}
	}
	n<-length(knots)
# number of bases	
	q<-n+order
	Aknots <- sort(c(rep(Boundary.knots, order), knots))
	SB <- orthogonalsplinebasis::SplineBasis(knots=Aknots, order=order, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	new("BSplineBasis", knots=knots, degree=as.integer(degree), nbases=as.integer(q), Matrices=M, SplineBasis=SB, log=log)
}

oldBSplineBasis<-function(knots, degree=3, keep.duplicates=FALSE, log=FALSE) {
	# for B-splines, knots are degree+1 boundary min knots, interior knots, degree+1 boundary max knots
	# boundary knots are not tested to be duplicated
	# code from orthogonalsplinebasis
	order<-degree+1
	n<-length(knots)
	if( n ==2*order ){
# no interior knots
		knots<-knots  # !! all the boundary knots are given
	}
	else {
		if(any(table(knots[order:(n-order+1)])>1)&&!keep.duplicates){
			warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
			# modif MGk 06/06/2011 to keep the (order-1) first and last knots 
			oknots <- knots
			iknots<-unique(oknots[order:(n-order+1)])
			knots<-c(oknots[1:(order-1)], iknots, oknots[(n-order+2):length(oknots)])
		}
	}
	q<-n-order
	SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=order, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	new("BSplineBasis", knots=knots, degree=as.integer(degree), nbases=as.integer(q), Matrices=M, SplineBasis=SB, log=log)
}

maketpdegrees <- function(knots, order){
	order - unlist(lapply(table(knots), function(x) x:1))
}

TPSplineBasis<-function(knots, degree=3, min, max, type=c("standard", "increasing"), coef=NULL, log=FALSE, keep.duplicates=FALSE) {
	# knots are interior knots
	type  <- match.arg(type)       # type of tp-spline
	order<-degree+1
	n <- length(knots)
	if( n ==0 ){
		# no interior knots
		knots <- as.numeric(knots)
	}
	else {
		knots <- sort(knots)
		if (!keep.duplicates & any(table(knots) > 1) ) {
			warning("Duplicate interior knots. Removing duplicates.\n")
			knots <- unique(knots)
			degrees <- rep(degree, length(knots))
		}
		else {
			degrees <- maketpdegrees(knots, order)
		}
	}
	n <- length(knots)
	if( is.null(coef)) {
		coef <- rep(1, n+order)
	}
	new("TPSplineBasis", knots=knots, degree=as.integer(degree), nbases=as.integer(n+order),
			min=min, max=max , coef=coef, degrees=as.integer(degrees),
			log=log, type=type)
}

TPRSplineBasis<-function(knots, degree=3, ref=0, min, max, type=c("standard", "increasing"), coef=NULL, log=FALSE, keep.duplicates=FALSE) {
	# knots are interior knots
	type  <- match.arg(type)       # type of tp-spline
	order<-degree+1
	n <- length(knots)
	if( n ==0 ){
		# no interior knots
		knots <- as.numeric(knots)
	}
	else {
		knots <- sort(knots)
		if (!keep.duplicates & any(table(knots) > 1) ) {
			warning("Duplicate interior knots. Removing duplicates.\n")
			knots <- unique(knots)
			degrees <- rep(degree, length(knots))
		}
		else {
			degrees <- maketpdegrees(knots, order)
		}
	}
	n <- length(knots)
	if( is.null(coef)) {
		coef <- rep(1, n+order)
	}
	new("TPSplineBasis", knots=knots, degree=as.integer(degree), nbases=as.integer(n+order),
			min=min, max=max , coef=coef, degrees=as.integer(degrees), ref=ref,
			log=log, type=type)
	
	new("TPSplineBasis", knots=knots, degree=as.integer(degree), nbases=as.integer(n+degree+1), ref=ref, min=min, max=max, log=log )
}

C0BSplineBasis<-function(knots, degree=3, keep.duplicates=FALSE, log=FALSE) {
	# B-splines, for constraints optimisation : sum(beta_i)_(i<= degree)b_i(0) = 1
	# for B-splines, knots are degree+1 boundary min knots, interior knots, degree+1 boundary max knots
	# code from orthogonalsplinebasis
	order<-degree+1
	n<-length(knots)
	if(any(table(knots[order:(n-order+1)])>1)&&!keep.duplicates){
		warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
		# modif MGk 06/06/2011 to keep the (order-1) first and last knots 
		oknots <- knots
		iknots<-unique(oknots[order:(n-order+1)])
		knots<-c(oknots[1:(order-1)], iknots, oknots[(n-order+2):length(oknots)])
	}
	q<-n-order
	SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=degree+1, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	new("C0BSplineBasis", knots=knots, degree=as.integer(degree), nbases=as.integer(q) , Matrices=M, log=log)
}



C0TPSplineBasis<-function(knots, degree=3, log=FALSE) {
	# duplicated notes are suppressed
	# for natural splines, knots are interior knots
	order<-degree+1
	n <- length(knots)
	if (any(table(knots[2:(n - 1)]) > 1) ) {
		warning("Duplicate interior knots. Removing duplicates.\n")
		knots <- unique(knots)
	}
	new("TPSplineBasis", knots=knots, degree=as.integer(degree), log=log)
}

###################################################################################
####### fin          createurs
###################################################################################


###################################################################################
#######         GETEUR
###################################################################################

setGeneric("getKnots",function(object)standardGeneric("getKnots"))
setMethod("getKnots",signature(".SplineBasis"),function(object)object@knots)

setGeneric("getInteriorKnots",function(object)standardGeneric("getInteriorKnots"))
setMethod("getInteriorKnots",signature(".SplineBasis"),function(object){
#			object@knots[-c(1:(object@degree+1), length(object@knots)-(0:object@degree))]
			object@knots[(object@degree+2):(length(object@knots)-object@degree-1)]
		}
)

setGeneric("getDegree",function(object)standardGeneric("getDegree"))
setMethod("getDegree",signature(".SplineBasis"),function(object)object@degree)

setGeneric("getOrder",function(object)standardGeneric("getOrder"))
setMethod("getOrder",signature(".SplineBasis"),function(object)object@degree+1)

setGeneric("getNBases",function(object)standardGeneric("getNBases"))
setMethod("getNBases",signature(".SplineBasis"),function(object)object@nbases)

dim.BSplineBasis<-function(x)dim(x@Matrices)
setMethod("dim","BSplineBasis", dim.BSplineBasis)
setMethod("dim","EBSplineBasis", dim.BSplineBasis)

setGeneric("getLog",function(object)standardGeneric("getLog"))
setMethod("getLog",signature(".SplineBasis"),function(object)object@log)

setGeneric("getSplineMatrix",function(object)standardGeneric("getSplineMatrix"))
setMethod("getSplineMatrix",signature("BSplineBasis"),function(object)object@Matrix)

setGeneric("getRef",function(object)standardGeneric("getRef"))
setMethod("getRef",signature("TPSplineBasis"),function(object)object@ref)

setGeneric("getMin",function(object)standardGeneric("getMin"))
setMethod("getMin",signature("TPSplineBasis"),function(object)object@min)
setMethod("getMin",signature("LEBSplineBasis"),function(object)object@min)

setGeneric("getMax",function(object)standardGeneric("getMax"))
setMethod("getMax",signature("TPSplineBasis"),function(object)object@max)
setMethod("getMax",signature("LEBSplineBasis"),function(object)object@max)

setGeneric("getInteriorKnots",function(object)standardGeneric("getInteriorKnots"))
setMethod("getInteriorKnots",signature("BSplineBasis"),
		function(object)
			object@knots[(1:(length(object@knots)-2*getOrder(object)))+getOrder(object)]
)
setMethod("getInteriorKnots",signature("BSplineBasis"),
		function(object)
			object@knots[(1:(length(object@knots)-2*getOrder(object)))+getOrder(object)]
)
setMethod("getInteriorKnots",signature("TPSplineBasis"),
		function(object) object@knots
)



###################################################################################
####### fin          GETEUR
###################################################################################

###################################################################################
#######         shower
###################################################################################


print.SplineBasis<-function(object) { 
	cat(class(object),"\n")
	cat("Order: ",object@degree+1,"\n",
			"Degree: ",object@degree,"\n",
			"Knots: ", paste(object@knots,collapse=" "),"\n",
			"Number of bases: ", object@nbases,"\n",
			"Range: ", paste(c(object@min, object@max),collapse=" ; "),"\n",
			sep="")
	invisible(object)
}

print.TPSplineBasis<-function(object) { 
	cat(class(object),"\n")
	cat("Order: ",object@degree+1,"\n",
			"Degree: ",object@degree+1,"\n",
			"Knots: ",paste(object@knots,collapse=" "),"\n",
			"Degrees: ",paste(object@degrees,collapse=" "),"\n",
			"Range: ", paste(c(object@min, object@max),collapse=" ; "),"\n",
			"Number of bases: ", object@nbases,"\n",
			"Type: ",object@type,"\n", sep="")
	invisible(object)
}

setMethod("show","BSplineBasis",  print.SplineBasis)
setMethod("show","BSplineBasis",  print.SplineBasis)
setMethod("show","EBSplineBasis",  print.SplineBasis)
setMethod("show","LEBSplineBasis",print.SplineBasis)
setMethod("show","R2bBSplineBasis",print.SplineBasis)
setMethod("show","TPSplineBasis", print.TPSplineBasis)


###################################################################################
#######   end      shower
###################################################################################


######################################################################
# using spline.des()
EvaluateBBasis0<-function(object,x,intercept=TRUE, xname=NULL,  ...) { 
	stopifnot(is.numeric(x))
	dots<-list(...)
	nx <- names(x)
	nax <- is.na(x)
	if (nas <- any(nax)){ 
		x <- x[!nax]
	}
	Aknots<-object@knots
	degree<-object@degree
	ol <- x < Aknots[degree+1]
	or <- x > Aknots[length(Aknots)-degree]
	outside <- ( or | ol)
	
	if (any(outside)) {
		basis <- array(, dim= c(length(x), object@nbases))
		if (any(inside <- !outside)){ 
			basis[inside, ] <- spline.des(knots=Aknots, x=x[inside], ord=degree+1)$design
		}
	}
	else {
		basis <- spline.des(knots=Aknots, x=x, ord=degree+1)$design
	}
	if (!intercept) {
		basis <- basis[, -1L, drop = FALSE]
	}
	
	if (object@log) {
		basis <- cbind(basis, log(x))
	}
	
	
	
#add na x
	n.col <- ncol(basis)
	if (nas) {
		nmat <- matrix(NA, length(nax), n.col)
		nmat[!nax, ] <- basis
		basis <- nmat
	}
	
	#add dimnames
	if(!is.null(xname)){
		if (intercept) {
			dbs <- paste("BS-", xname, 1:(getNBases(object)+getLog(object)), sep="")
		}
		else {
			dbs <- paste("BS-", xname, 2:(getNBases(object)+getLog(object)), sep="")
		}
		dimnames(basis) <- list(nx, dbs)
	}
	
	
	
	a <- list(degree = degree, knots =  Aknots, 
			intercept = intercept, log=getLog(object))
	attributes(basis) <- c(attributes(basis), a)
	class(basis) <- c("bs", "basis", "matrix")
	basis
}


FEvaluateBBasis<-function(object, x, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(x))
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	basis <- .Call(C_eval_spline_basis, as.double(knots), as.integer(order), M, 
			as.integer(intercept), as.double(x), as.integer(outer.ok))
	
	if (object@log) {
		return(cbind(basis, log(x)))
	}
	else {
		return(basis)
	}
	
}

EvaluateBBasis<-function(object,x,intercept=TRUE, xname=NULL, outer.ok=TRUE, namespline= "B", ...) {
	
	basis <- FEvaluateBBasis(object=object, x=x,
			intercept=intercept, 
			xname=xname,
			outer.ok= outer.ok, ...) 
	
	Aknots<-object@knots
	degree<-object@degree
	nbases<-object@nbases
	
	
	#add dimnames
	if(!is.null(xname)){
		if (intercept) {
			dbs <- paste(namespline, "-", xname, ":", 1:(getNBases(object)+getLog(object)), sep="")
		}
		else {
			dbs <- paste(namespline, "-", xname, ":", 2:(getNBases(object)+getLog(object)), sep="")
		}
		dimnames(basis)[[2]] <- dbs
	}
	
	
#  dimnames(basis) <- list(nx, 1L:n.col)
	a <- list(degree = degree, knots =  Aknots, 
			intercept = intercept, log=getLog(object))
	attributes(basis) <- c(attributes(basis), a)
	class(basis) <- c("bs", "basis", "matrix")
	basis
}


FEvaluateEBBasis<-function(object, x, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(x))
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	basis <- .Call(C_eval_spline_basis, as.double(knots), as.integer(order), M, 
			as.integer(intercept), as.double(x), as.integer(outer.ok))
	
	if (object@log) {
		return(cbind(basis, log(x)))
	}
	else {
		return(basis)
	}
	
}


EvaluateEBBasis<-function(object,x,intercept=TRUE, xname=NULL, outer.ok=TRUE, namespline= "B", ...) {
	
	basis <- FEvaluateEBBasis(object=object, x=x,
			intercept=intercept, 
			xname=xname,
			outer.ok= outer.ok, ...) 
	
	Aknots<-object@knots
	degree<-object@degree
	nbases<-object@nbases
	
	
	#add dimnames
	if(!is.null(xname)){
		if (intercept) {
			dbs <- paste(namespline, "-", xname, ":", 1:(getNBases(object)+getLog(object)), sep="")
		}
		else {
			dbs <- paste(namespline, "-", xname, ":", 2:(getNBases(object)+getLog(object)), sep="")
		}
		dimnames(basis)[[2]] <- dbs
	}
	
	
#  dimnames(basis) <- list(nx, 1L:n.col)
	a <- list(degree = degree, knots =  Aknots, 
			intercept = intercept, log=getLog(object))
	attributes(basis) <- c(attributes(basis), a)
	class(basis) <- c("bs", "basis", "matrix")
	basis
}

# evaluate linearly extended Bspline
FEvaluateLEBBasis<-function(object, x, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(x))
	dots<-list(...)
	M<-object@Matrices
	linexinf<-object@linexinf
	linexsup<-object@linexsup
	knots<-object@knots
	order<-object@degree+1L
	orderextrapol <- object@orderextrapol
	basis <- .Call(C_eval_linex_spline_basis, as.double(knots), as.integer(order), M, linexinf, linexsup,
			as.integer(orderextrapol), as.integer(intercept), as.double(x), as.integer(outer.ok))
	
	if (object@log) {
		return(cbind(basis, log(x)))
	}
	else {
		return(basis)
	}
	
}

EvaluateLEBBasis<-function(object, x, intercept=TRUE, xname=NULL, outer.ok=TRUE, namespline= "B", ...) {
	basis <- FEvaluateLEBBasis(object=object, x=x,
			intercept=intercept, 
			outer.ok=outer.ok) 
	
	Aknots<-object@knots
	degree<-object@degree
	nbases<-object@nbases
	
	
	#add dimnames
	if(!is.null(xname)){
		if (intercept) {
			dbs <- paste(namespline, "-", xname, ":", 1:(getNBases(object)+getLog(object)), sep="")
		}
		else {
			dbs <- paste(namespline, "-", xname, ":", 2:(getNBases(object)+getLog(object)), sep="")
		}
		dimnames(basis)[[2]] <- dbs
	}
	
	
#  dimnames(basis) <- list(nx, 1L:n.col)
	a <- list(degree = degree, knots =  Aknots, 
			intercept = intercept, log=getLog(object))
	attributes(basis) <- c(attributes(basis), a)
	class(basis) <- c("bs", "basis", "matrix")
	basis
}



# evaluate tp Basis, no ndimnames
FEvaluateTPBasisPosOld<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
	stopifnot(is.numeric(x))
	nax <- is.na(x)
	if (nas <- any(nax)){ 
		x <- x[!nax]
	}
	
	Aknots<-object@knots
	degree<-object@degree
	if (xmin != -Inf){
		ol <- x < xmin
		if (xmax != +Inf){
			or <- x > xmax
			outside <- ( or | ol)
		}
		else {
			outside <- ol
		}
	}
	else {
		if (xmax != +Inf){
			outside <- x > xmax
		}
		else {
			outside <- FALSE
		}
	}
	
	Aknotref <- Aknots 
	if( !is.null(ref)){
		x <- x - ref
		if(length(Aknots)>0){
			Aknotref <- Aknots - ref
		}
	}
	
	
	if (intercept) {
		nb <- degree+ 1 + length(Aknots) 
		last <- degree+1
	}
	else{
		nb <- degree  + length(Aknots) 
		last <- degree
	}    
	if( outer.ok) {
		outervalue <- 0.0
	}
	else {
		outervalue <- NA
	}
	
	basis <- matrix(1, nrow=length(x), ncol=nb)
	if (any(!outside)){ 
		if (intercept) {
			for(i in 1:degree){
				basis[,i+1]<-x^i
			}
		}
		else{
			for(i in 1:degree){
				basis[,i]<-x^i
			}
		}
		if(length(Aknots)>0){
			for (k in 1:(length(Aknots))){
				basis[,last+k] <- ifelse(x>Aknotref[k], (x-Aknotref[k])^degree, 0)
			}
		}
		if (object@log) {
			basis <- cbind(basis, log(x))
		}
	}    
	
#add outside values
	if(any(outside)){
		basis[outside,]<-outervalue
	}
	
#add na x
	if (nas) {
		nabasis <- matrix(NA, nrow=length(nax), ncol=nb)
		nabasis[!nax, ] <- basis
		nabasis
	}
	else {
		basis
	}
}

FEvaluateTPBasisPos<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
	stopifnot(is.numeric(x))
	
	min<-object@min
	max<-object@max
	allknots<-object@knots
	order<-object@degree+1
	knots <- unique(allknots)
	
	if( !is.null(ref)){
		x <- x - ref
		min <- min- ref
		max <- max - ref
		if(length(knots)>0){
			knot <- knots - ref
		}
	}
	
	replicates <- table(allknots)
	degrees <- object@degrees
	coef <- object@coef
	basis <- .Call(C_eval_trunc_power_basis, as.double(knots), as.double(replicates),
			as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
			as.double(x), as.integer(outer.ok))
	
	return(basis)
	
}

# evaluate TP basis, add dimnames
EvaluateTPBasisPos<-function(object,x,intercept=TRUE, ref=NULL, xname=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) {
	
	basis <- FEvaluateTPBasisPos(object=object, x=x,
			intercept=intercept, ref=ref,
			xname=xname,
			outer.ok=outer.ok, ...) 
	
	
	Aknots<-object@knots
	degree<-object@degree
	
	#add dimnames
	if(!is.null(xname)){
		if(degree>2 ){
			if( !is.null(ref)){
				dnva <- c(paste("(",xname, "-", ref, ")", sep=""),
						paste("(",xname, "-", ref, ")^", 2:degree, sep=""))
			}
			else {
				dnva <- c(xname, paste(xname, "^", 2:degree, sep=""))
			}
		}
		else {
			if( !is.null(ref)){
				dnva <- c(paste("(",xname, "-", ref, ")", sep=""))
			}
			else {
				dnva <- c(xname)
			}
		}
		
		if (intercept) {
			dnva <- c("Intercept", dnva)
		}
		
		if(length(Aknots)>0){
			dnvaplus <- paste("(",xname, "-", Aknots[1:(length(Aknots))], ")_+^", degree, sep="")
		}
		if (getLog(object)) {
			dlog  <- paste("log(",xname, ")", sep="")
		}
		else {
			dlog  <- NULL
		}
		dimnames(basis)[[2]]<-c(dnva, dnvaplus, dlog)
	}
	
	a <- list(degree = degree, knots =  Aknots, 
			intercept = intercept, ref=ref, log=getLog(object))
	attributes(basis) <- c(attributes(basis), a)
	class(basis) <- c("tp", "basis", "matrix")
	basis
}


# similar to FEvaluateTPBasisPos but if knot < ref, -(k-x)+^d if x<ref
FEvaluateTPBasisOld<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
	stopifnot(is.numeric(x))
	nax <- is.na(x)
	if (nas <- any(nax)){ 
		x <- x[!nax]
	}
	
	Aknots<-object@knots
	degree<-object@degree
	if (xmin != -Inf){
		ol <- x < xmin
		if (xmax != +Inf){
			or <- x > xmax
			outside <- ( or | ol)
		}
		else {
			outside <- ol
		}
	}
	else {
		if (xmax != +Inf){
			outside <- x > xmax
		}
		else {
			outside <- FALSE
		}
	}
	
}

FEvaluateTPBasisstd<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
	stopifnot(is.numeric(x))
	
	min<-object@min
	max<-object@max
	allknots<-object@knots
	order<-object@degree+1
	knots <- unique(allknots)
	
	if( !is.null(ref)){
		x <- x - ref
		min <- min- ref
		max <- max - ref
		if(length(knots)>0){
			knot <- knots - ref
		}
	}
	
	replicates <- table(allknots)
	degrees <- object@degrees
	coef <- object@coef
	basis <- .Call(C_eval_trunc_power_increasing_basis, as.double(knots), as.double(replicates),
			as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
			as.double(x), as.integer(outer.ok))
	
	return(basis)
}


# truncated powezr basis
FEvaluateTPBasis<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
	stopifnot(is.numeric(x))
	
	min<-object@min
	max<-object@max
	allknots<-object@knots
	order<-object@degree+1
	knots <- unique(allknots)
	
	if( !is.null(ref)){
		x <- x - ref
		min <- min- ref
		max <- max - ref
		if(length(knots)>0){
			knot <- knots - ref
		}
	}
	
	replicates <- table(allknots)
	degrees <- object@degrees
	coef <- object@coef
	
	if(object@type == "increasing"){
		basis <- .Call(C_eval_trunc_power_increasing_basis, as.double(knots), as.double(replicates),
				as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
				as.double(x), as.integer(outer.ok))
	}
	else {
		basis <- .Call(C_eval_trunc_power_basis, as.double(knots), as.double(replicates),
				as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
				as.double(x), as.integer(outer.ok))
	}
	if (object@log) {
		basis <- cbind(basis, log(x))
	}
	return(basis)
}



# same as FEvaluateTPBasis, but with dimnames, attribures, ...
EvaluateTPBasis<-function(object,x,intercept=TRUE, ref=NULL, xname=NULL, 
		outer.ok=TRUE, ...) { 
	
	
	basis <- FEvaluateTPBasis(object=object, x=x,
			intercept=intercept, ref=ref,
			xname=xname,
			outer.ok= outer.ok, ...)  
	
	Aknots<-object@knots
	degree<-object@degree
	
	#add dimnames
	if(!is.null(xname)){
		if(degree>2 ){
			if( !is.null(ref)){
				dnva <- c(paste("(",xname, "-", ref, ")", sep=""),
						paste("(",xname, "-", ref, ")^", 2:degree, sep=""))
			}
			else {
				dnva <- c(xname, paste(xname, "^", 2:degree, sep=""))
			}
		}
		else {
			if( !is.null(ref)){
				dnva <- c(paste("(",xname, "-", ref, ")", sep=""))
			}
			else {
				dnva <- c(xname)
			}
		}
		
		if (intercept) {
			dnva <- c("Intercept", dnva)
		}
		
		if(length(Aknots)>0){
			dnvaplus <- paste("(",xname, "-", Aknots[1:(length(Aknots))], ")_+^", degree, sep="")
		}
		if (getLog(object)) {
			dlog  <- paste("log(",xname, ")", sep="")
		}
		else {
			dlog  <- NULL
		}
		dimnames(basis)[[2]]<-c(dnva, dnvaplus, dlog)
	}
	
	a <- list(degree = degree, knots =  Aknots, 
			intercept = intercept, ref=ref, log=getLog(object))
	attributes(basis) <- c(attributes(basis), a)
	class(basis) <- c("tp", "basis", "matrix")
	basis
}




######################################################################
EvaluateC0SBasis<-function(object,x,intercept=TRUE, ...) { 
# output : b_1(t) , b_i(t) - b_1(t) if i in 2:(degree), b_i(t) if i > degree
	stopifnot(is.numeric(x))
	dots<-list(...)
	nx <- names(x)
	nax <- is.na(x)
	if (nas <- any(nax)){ 
		x <- x[!nax]
	}
	
	Aknots<-object@knots
	degree<-object@degree
	ol <- x < Aknots[degree+1]
	or <- x > Aknots[length(Aknots)-degree]
	outside <- ( or | ol)
	
	
	if (any(outside)) {
		basis <- array(NA, dim= c(length(x), length(Aknots) - degree - 1L))
		if (any(inside <- !outside)){ 
			thebasis <- spline.des(knots=Aknots, x=x[inside], ord=degree+1)$design
			for( i in 2:degree){
				thebasis[,i]<-thebasis[,i] - thebasis[,1]
			}
			basis[inside, ] <- thebasis
		}
	}
	else {
		basis <- spline.des(knots=Aknots, x=x, ord=degree+1)$design
		for( i in 2:degree){
			basis[,i]<-basis[,i] - basis[,1]
		}
	}
	if (!intercept) {
		basis <- basis[, -1L, drop = FALSE]
	}
	if (object@log) {
		basis <- cbind(basis, log(x))
	}
	
#add na x
	if (nas) {
		n.col <- ncol(basis)
		nmat <- matrix(NA, length(nax), n.col)
		nmat[!nax, ] <- basis
		basis <- nmat
	}
	dimnames(basis) <- list(nx, 1L:n.col)
	a <- list(degree = degree, knots =  Aknots, 
			intercept = intercept, log=getLog(object))
	attributes(basis) <- c(attributes(basis), a)
	class(basis) <- c("bs", "basis", "matrix")
	basis
}

# here, evaluate for tp-spline is EvaluateTPBasis
# thus, evaluate(tpspline, 0 , intercept = TRUE) == c(1, rep(0, nn)) 

#setMethod("evaluate",signature("BSplineBasis","numeric"),function(object, x, ...)EvaluateBBasis(object=object, x=x, ...))
setMethod("evaluate",signature("BSplineBasis","numeric"),function(object, x, ...)EvaluateBBasis(object=object, x=x, ...))
setMethod("evaluate",signature("EBSplineBasis","numeric"),function(object, x, ...)EvaluateEBBasis(object=object, x=x, ...))
setMethod("evaluate",signature("LEBSplineBasis","numeric"),function(object, x, ...)EvaluateLEBBasis(object=object, x=x, ...))
setMethod("evaluate",signature("R2bBSplineBasis","numeric"),function(object, x, ...)EvaluateLEBBasis(object=object, x=x, ...))
setMethod("evaluate",signature("TPSplineBasis","numeric"),function(object, x, ...)EvaluateTPBasis(object=object, x=x, ...))
setMethod("evaluate",signature("TPRSplineBasis","numeric"),function(object, x, ...)EvaluateTPBasis(object=object, x=x, ref=object@ref, ...))
setMethod("evaluate",signature("C0BSplineBasis","numeric"),function(object, x, ...)EvaluateC0SBasis(object=object, x=x, ...))

setMethod("fevaluate",signature("BSplineBasis","numeric"),function(object, x, ...)FEvaluateBBasis(object=object, x=x, ...))
setMethod("fevaluate",signature("EBSplineBasis","numeric"),function(object, x, ...)FEvaluateEBBasis(object=object, x=x, ...))
setMethod("fevaluate",signature("LEBSplineBasis","numeric"),function(object, x, ...)FEvaluateLEBBasis(object=object, x=x, ...))
setMethod("fevaluate",signature("R2bBSplineBasis","numeric"),function(object, x, ...)FEvaluateLEBBasis(object=object, x=x, ...))
setMethod("fevaluate",signature("TPSplineBasis","numeric"),function(object, x, ...)FEvaluateTPBasis(object=object, x=x, ...))
setMethod("fevaluate",signature("TPRSplineBasis","numeric"),function(object, x, ...)EvaluateTPBasis(object=object, x=x, ref=object@ref, ...))
setMethod("fevaluate",signature("C0BSplineBasis","numeric"),function(object, x, ...)EvaluateC0SBasis(object=object, x=x, ...))
setMethod("fevaluate",signature("NULL","numeric"),function(object, x, ...) NULL)


######################################################################
#  evaluate linear combination of spline basis

# M-spline
EvaluateLCBBasis<-function(object, x, beta, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(x))
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	cl <- .Call(C_eval_lc_spline_basis, as.double(knots), as.integer(order), M,
			as.integer(intercept), as.double(x), as.double(beta), as.integer(outer.ok))
	
	if (object@log) {
		return(cl + log(x) * beta[length(beta)])
	}
	else {
		return(cl)
	}
	
}

# linearly extended  M-spline
EvaluateLCLEBBasis<-function(object, x, beta, intercept=TRUE, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(x))
	dots<-list(...)
	M<-object@Matrices
	linexinf<-object@linexinf
	linexsup<-object@linexsup
	knots<-object@knots
	order<-object@degree+1
	cl <- .Call(C_eval_lc_linex_spline_basis, as.double(knots), as.integer(order), M, linexinf, linexsup,
			as.integer(intercept), as.double(x), as.double(beta), as.integer(outer.ok))
	
	if (object@log) {
		return(cl + log(x) * beta[length(beta)])
	}
	else {
		return(cl)
	}
	
}

EvaluateLCTPBasis<-function(object, x, beta, intercept=TRUE, ref=NULL, outer.ok=TRUE, ...) { 
	stopifnot(is.numeric(x))
	
	min<-object@min
	max<-object@max
	allknots<-object@knots
	order<-object@degree+1
	knots <- unique(allknots)
	
	if( !is.null(ref)){
		x <- x - ref
		min <- min- ref
		max <- max - ref
		if(length(knots)>0){
			knot <- knots - ref
		}
	}
	
	replicates <- table(allknots)
	degrees <- object@degrees
	coef <- object@coef
	if(object@type == "increasing"){
		cl <- .Call(C_eval_lc_trunc_power_increasing_basis, as.double(knots), as.double(replicates),
				as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
				as.double(x), as.double(beta), as.integer(outer.ok))
	}
	else {
		cl <- .Call(C_eval_lc_trunc_power_basis, as.double(knots), as.double(replicates),
				as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
				as.double(x), as.double(beta), as.integer(outer.ok))
	}
	if (object@log) {
		return(cl + log(x) * beta[length(beta)])
	}
	else {
		return(cl)
	}
}





setMethod("evaluatelc",signature("BSplineBasis","numeric","numeric"),function(object, x, beta, ...)EvaluateLCBBasis(object=object, x=x, beta=beta, ...))
setMethod("evaluatelc",signature("LEBSplineBasis","numeric","numeric"),function(object, x, beta, ...)EvaluateLCLEBBasis(object=object, x=x, beta=beta, ...))
setMethod("evaluatelc",signature("R2bBSplineBasis","numeric","numeric"),function(object, x, beta, ...)EvaluateLCLEBBasis(object=object, x=x, beta=beta, ...))
setMethod("evaluatelc",signature("TPSplineBasis","numeric","numeric"),function(object, x, beta, ...)EvaluateLCTPBasis(object=object, x=x, beta=beta, ...))


######################################################################
#  Method predic for SplineParam

# MSplineParam
# computes f(x) = sum_i beta_i b_i(x)
PredictBBasisBeta <- function(object=object, x=x, beta=beta, intercept=TRUE, outer.ok=TRUE, ...){
	if(intercept){
		predict(object * beta, x, intercept=TRUE, ...)
	}
	else {
		predict(object * c(0, beta), x,  intercept=FALSE, ...)
	}
}

# EMSplineParam
# computes f(x) = sum_i beta_i b_i(x)
PredictEBBasisBeta <- function(object=object, x=x, beta=beta, intercept=TRUE, outer.ok=TRUE, ...){
	if(intercept){
		predict(object * beta, x, intercept=TRUE, ...)
	}
	else {
		predict(object * c(0, beta), x,  intercept=FALSE, ...)
	}
}

# LEMSplineParam
# computes f(x) = sum_i beta_i b_i(x)
PredictLEBBasisBeta <- function(object=object, x=x, beta=beta, intercept=TRUE, outer.ok=TRUE, ...){
	if(intercept){
		predict(object * beta, x, intercept=TRUE, ...)
	}
	else {
		predict(object * c(0, beta), x,  intercept=FALSE, ...)
	}
}

# computes f(x) = sum_i B_i(x)
# assuming B_i = beta_i * b_i
predict.BSplineBasis <- function(object=object, x=x, intercept=TRUE, outer.ok=TRUE, ...){
	stopifnot(is.numeric(x))
	
#  if(!is.null(beta)){
#    if(intercept){
#      object <- object * beta
#    }
#    else {
#      object <- object * c(0, beta)
#    }
#  }
	
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	cl <- .Call(C_predict_spline_basis, as.double(knots), as.integer(order), M,
			as.integer(intercept), as.double(x), as.integer(outer.ok))
	
	if (object@log) {
		return(cl + log(x) * beta[length(beta)])
	}
	else {
		return(cl)
	}
	
}

# computes f(x) = sum_i B_i(x)
# assuming B_i = beta_i * b_i
slowpredict.BSplineBasis <- function(object=object, x=x, intercept=TRUE, outer.ok=TRUE, ...){
	stopifnot(is.numeric(x))
	
#  if(!is.null(beta)){
#    if(intercept){
#      object <- object * beta
#    }
#    else {
#      object <- object * c(0, beta)
#    }
#  }
	
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	cl <- .Call(C_slow_predict_spline_basis, as.double(knots), as.integer(order), M,
			as.integer(intercept), as.double(x), as.integer(outer.ok))
	
	if (object@log) {
		return(cl + log(x) * beta[length(beta)])
	}
	else {
		return(cl)
	}
	
}

# EBSplineBasis
# computes f(x) = sum_i B_i(x)
# assuming B_i = beta_i * b_i
predict.EBSplineBasis <- function(object=object, x=x, intercept=TRUE, outer.ok=TRUE, ...){
	stopifnot(is.numeric(x))
	
	dots<-list(...)
	M<-object@Matrices
	knots<-object@knots
	order<-object@degree+1
	cl <- .Call(C_predict_spline_basis, as.double(knots), as.integer(order), M,
			as.integer(intercept), as.double(x), as.integer(outer.ok))
	
	if (object@log) {
		return(cl + log(x) * beta[length(beta)])
	}
	else {
		return(cl)
	}
	
}

# computes f(x) = sum_i B_i(x)
# assuming B_i = beta_i * b_i
predict.LEBSplineBasis <- function(object=object, x=x, intercept=TRUE, outer.ok=TRUE, ...){
	stopifnot(is.numeric(x))
	
#  if(!is.null(beta)){
#    if(intercept){
#      object <- object * beta
#    }
#    else {
#      object <- object * c(0, beta)
#    }
#  }
	
	dots<-list(...)
	M<-object@Matrices
	linexinf<-object@linexinf
	linexsup<-object@linexsup
	orderextrapol <- object@orderextrapol
	knots<-object@knots
	order<-object@degree+1
	cl <- .Call(C_predict_linex_spline_basis, as.double(knots), as.integer(order), M, linexinf, linexsup,
			as.integer(orderextrapol), as.integer(intercept), as.double(x), as.integer(outer.ok))
	
	if (object@log) {
		return(cl + log(x) * beta[length(beta)])
	}
	else {
		return(cl)
	}
	
}


slowpredict.LEBSplineBasis <- function(object=object, x=x, intercept=TRUE, outer.ok=TRUE, ...){
	stopifnot(is.numeric(x))
	
#  if(!is.null(beta)){
#    if(intercept){
#      object <- object * beta
#    }
#    else {
#      object <- object * c(0, beta)
#    }
#  }
	
	dots<-list(...)
	M<-object@Matrices
	linexinf<-object@linexinf
	linexsup<-object@linexsup
	knots<-object@knots
	order<-object@degree+1
	orderextrapol <- object@orderextrapol
	cl <- .Call(C_slow_predict_linex_spline_basis, as.double(knots), as.integer(order), M, linexinf, linexsup,
			as.integer(orderextrapol), as.integer(intercept), as.double(x), as.integer(outer.ok))
	
	if (object@log) {
		return(cl + log(x) * beta[length(beta)])
	}
	else {
		return(cl)
	}
	
}


# TPSplineParam
# computes f(x) = sum_i beta_i b_i(x)
PredictTPBasisBeta <- function(object=object, x=x, beta=beta, intercept=TRUE, ref=NULL, outer.ok=TRUE, ...){
	if(intercept){
		predict(object * beta, x, intercept=TRUE, ...)
	}
	else {
		predict(object * c(0, beta), x,  intercept=FALSE, ...)
	}
}

# computes f(x) = sum_i beta_i * B_i(x)
# if beta == NULL, assuming beta_i = 1
predict.TPSplineBasis <- function(object=object, x=x, beta=NULL, intercept=TRUE, ref=NULL, outer.ok=TRUE, ...){
	stopifnot(is.numeric(x))
	
	if(!is.null(beta)){
		if(intercept){
			object <- object * beta
		}
		else {
			object <- object * c(0, beta)
		}
	}
	min<-object@min
	max<-object@max
	allknots<-object@knots
	order<-object@degree+1
	knots <- unique(allknots)
	
	if( !is.null(ref)){
		x <- x - ref
		min <- min- ref
		max <- max - ref
		if(length(knots)>0){
			knot <- knots - ref
		}
	}
	
	replicates <- table(allknots)
	degrees <- object@degrees
	coef <- object@coef
	if(object@type == "increasing"){
		cl <- .Call(C_predict_trunc_power_increasing_basis, as.double(knots), as.double(replicates),
				as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
				as.double(x), as.integer(outer.ok))
	}
	else {
		cl <- .Call(C_predict_trunc_power_basis, as.double(knots), as.double(replicates),
				as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
				as.double(x), as.integer(outer.ok))
	}
	if (object@log) {
		return(cl + log(x) * beta[length(beta)])
	}
	else {
		return(cl)
	}
	
}


setMethod("predictSpline",signature(object="BSplineBasis",x="numeric"),function(object, x, ...)predict.BSplineBasis(object=object, x=x,  ...))
setMethod("predictSpline",signature(object="EBSplineBasis",x="numeric"),function(object, x, ...)predict.EBSplineBasis(object=object, x=x,  ...))
setMethod("predictSpline",signature(object="LEBSplineBasis",x="numeric"),function(object, x, ...)predict.LEBSplineBasis(object=object, x=x,  ...))
setMethod("predictSpline",signature(object="R2bBSplineBasis",x="numeric"),function(object, x, ...)predict.LEBSplineBasis(object=object, x=x,  ...))
setMethod("predictSpline",signature(object="TPSplineBasis",x="numeric"),function(object, x, ...)predict.TPSplineBasis(object=object, x=x,  ...))

#setMethod("predict",signature("BSplineBasis","numeric","numeric"),function(object, x, beta, ...)PredictBBasisBeta(object=object, x=x, beta=beta, ...))
#setMethod("predict",signature("TPSplineBasis","numeric","numeric"),function(object, x, beta, ...)PredictTPBasisBeta(object=object, x=x, beta=beta, ...))

setMethod("slowpredictSpline",signature(object="BSplineBasis",x="numeric", beta="missing"),function(object, x, ...)slowpredict.BSplineBasis(object=object, x=x,  ...))
setMethod("slowpredictSpline",signature(object="LEBSplineBasis",x="numeric", beta="missing"),function(object, x, ...)slowpredict.LEBSplineBasis(object=object, x=x,  ...))

######################################################################
#  integrate 

# define parameters for integrated Spline Basis
integrate.TPSplineBasis<-function(object){
	
	min<-object@min
	max<-object@max
	allknots<-object@knots
	order<-object@degree+1
	knots <- unique(allknots)
	replicates <- table(allknots)
	
	
	
	idegrees <- object@degrees + 1
	coef <- object@coef
	nbases <- length(coef)
	
	if(object@type == "standard"){
		icoef <- c(0, coef[1:order] / (1:order), coef[(order+1):nbases] / idegrees)
	}
	else {
		icoef <- c(0, coef[1:order] / (1:order), coef[(order+1):nbases] / ifelse(allknots<0, -  idegrees,   idegrees))
	}
	new("TPSplineBasis", knots=object@knots, min=object@min, max=object@max, 
			degree=object@degree+1, nbases=nbases, coef=icoef, degrees=idegrees, log=FALSE, type=object@type)
}

# define parameters for integrated Spline Basis
integrate.BSplineBasis<-function(object){
	SB <- orthogonalsplinebasis::integrate(object@SplineBasis)
	new("BSplineBasis", knots=SB@knots, min=object@min, max=object@max,
			degree=as.integer(SB@order-1), nbases=as.integer(dim(SB@Matrices)[2]), Matrices=SB@Matrices, SplineBasis=SB, log=FALSE)
}

# define parameters for integrated Spline Basis
integrate.MSplineBasis<-function(object){
	SB <- orthogonalsplinebasis::integrate(object@SplineBasis)
	new("MSplineBasis", knots=SB@knots, min=object@min, max=object@max,
			degree=as.integer(SB@order-1), nbases=as.integer(dim(SB@Matrices)[2]), Matrices=SB@Matrices, SplineBasis=SB, log=FALSE)
}

# compute values of integrated basis
IntegrateBBasis<-function(object,x,intercept=TRUE, xname=NULL, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(x))
	if (object@log) {
		stop("no method 'integrate' for BSplineBasis with additional log basis" )
	}
	else {
		evaluate(integrate(object), x=x, intercept=intercept, outer.ok=outer.ok, ...)
	}
}

# EBSplineBasis
# define parameters for integrated Spline Basis
integrate.EBSplineBasis<-function(object, ...){
	dnew<-d<-dim(object)
	M<-object@Matrices
	order<-getOrder(order)
	knots<- getKnots(knots)
	nknots<-length(knots)
	
	dnew[1]=d[1]+1L
	w<-c(diff(knots)[order:(nknots-order)], 1)
	
	N<-array(0,dim=dnew)
	for(i in 1:d[3]){
		N[,,i]<-rbind( if(i>1)(rep(1,order+1))%*%N[,,i-1] else 0,w[i]*t(diag(1/(1:order)))%*%M[,,i] )
	}
	
	newknots<-knots[c(1,seq(nknots),nknots)]
	neworder <- order+1L
	new("EBSplineBasis", knots=newknots, min=object@min, max=object@max,
			degree=as.integer(neworder-1), nbases=getNBases(object), Matrices=N, log=FALSE)
}

# compute values of integrated basis
IntegrateEBBasis<-function(object,x,intercept=TRUE, xname=NULL, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(x))
	if (object@log) {
		stop("no method 'integrate' for EBSplineBasis with additional log basis" )
	}
	else {
		evaluate(integrate(object), x=x, intercept=intercept, outer.ok=outer.ok, ...)
	}
}

# LEMSplinebasis
# define parameters for integrated Spline Basis
integrate.LEBSplineBasis<-function(object){
	SB <- orthogonalsplinebasis::integrate(object@SplineBasis)
	
	no <- new("LEBSplineBasis", knots=SB@knots, min=object@min, max=object@max,
			degree=as.integer(SB@order-1), nbases=as.integer(dim(SB@Matrices)[2]), Matrices=SB@Matrices, SplineBasis=SB, log=FALSE)
	linexinf <- rbind(rep(0, object@nbases), diag(1/(1:(object@degree+1)))%*%object@linexinf)
	linexsup <- rbind(evaluate(no, object@knots[length(object@knots)]), diag(1/(1:(object@degree+1)))%*%object@linexsup)
	no@linexinf <- linexinf
	no@linexsup <- linexsup
	no
}

# compute values of integrated basis
IntegrateLEBBasis<-function(object,x,intercept=TRUE, xname=NULL, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(x))
	if (object@log) {
		stop("no method 'integrate' for BSplineBasis with additional log basis" )
	}
	else {
		evaluate(integrate(object), x=x, intercept=intercept, outer.ok=outer.ok, ...)
	}
}



# compute values of integrated basis
IntegrateBBasisOld<-function(object,x,intercept=TRUE, xname=NULL, outer.ok=TRUE, ...) {
	stopifnot(is.numeric(x))
	nax <- is.na(x)
	if (nas <- any(nax)){ 
		x <- x[!nax]
	}
	
	Aknots<-object@knots
	degree<-object@degree
	nbases<-object@nbases
	ol <- x < Aknots[degree+1]
	or <- x > Aknots[length(Aknots)-degree]
	outside <- ( or | ol)
	
	if (any(outside)) {
		if(outer.ok) {
			basis <- array(0, dim= c(length(x), nbases))
		}
		else {
			basis <- array(NA, dim= c(length(x), nbases))
		}
		if (any(inside <- !outside)){ 
			basis[inside, ] <- orthogonalsplinebasis::evaluate(orthogonalsplinebasis::integrate(object@SplineBasis), x[inside] )
		}
	}
	else {
		basis <- orthogonalsplinebasis::evaluate(orthogonalsplinebasis::integrate(object@SplineBasis), x )
	}
	if (!intercept) {
		basis <- basis[, -1L, drop = FALSE]
	}
	
	if (object@log) {
		basis <- cbind(basis, 1/x)
	}
	
#add na x
	if (nas) {
		nabasis <- matrix(NA, length(nax), nbases)
		nabasis[!nax, ] <- basis
		nabasis
	}
	else {
		basis
	}
	
	
}

# compute values of integrated basis
# truncated powezr basis
integrateTPBasis<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
	stopifnot(is.numeric(x))
	if (object@log) {
		stop("no method 'integrate' for BSplineBasis with additional log basis" )
	}
	else {
		evaluate(integrate(object), x=x, intercept=intercept, ref=ref, xmin=xmin, xmax=xmax, outer.ok=outer.ok, ...)
	}
}

FIntegrateTPBasis0<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, ...) { 
# similar to FEvaluateTPBasis but if knots<ref, non 0 if x<0
# thus, active bases around ref are the full monomials (x-ref)^k
	stopifnot(is.numeric(x))
	nax <- is.na(x)
	if (nas <- any(nax)){ 
		x <- x[!nax]
	}
	
	Aknots<-object@knots
	degree<-object@degree
	degreep1 <- degree+1
	if (xmin != -Inf){
		ol <- x < xmin
		if (xmax != +Inf){
			or <- x > xmax
			outside <- ( or | ol)
		}
		else {
			outside <- ol
		}
	}
	else {
		if (xmax != +Inf){
			outside <- x > xmax
		}
		else {
			outside <- FALSE
		}
	}
	
	Aknotref <- Aknots 
	if( !is.null(ref)){
		x <- x - ref
		if(length(Aknots)>0){
			Aknotref <- Aknots - ref
		}
	}
	
	
	if (intercept) {
		nb <- degree+ 1 + length(Aknots) 
		last <- degree+1
	}
	else{
		nb <- degree  + length(Aknots) 
		last <- degree
	}    
	basis <- matrix(1, nrow=length(x), ncol=nb)
	
	
	
	if (any(!outside)){ 
		if (intercept) {
			for(i in 1:(degree+1)){
				basis[,i]<-x^i/i
			}
		}
		else{
			for(i in 2:(degree+1)){
				basis[,i+1]<-x^i/i
			}
		}
		if(length(Aknots)>0){
			for (k in 1:(length(Aknots))){
				if(Aknots[k]>0){
					basis[,last+k] <- ifelse(x>Aknotref[k], (x-Aknotref[k])^degreep1/degreep1, 0)
				} else {
					basis[,last+k] <- ifelse(x>Aknotref[k], 0, (x-Aknotref[k])^degreep1/degreep1)
				}
			}
		}
		if (object@log) {
			basis <- cbind(basis, log(x))
		}
	}    
	if (any(outside)) {
		basis[outside,]<-rep(NA, dim(basis)[2])
		if (object@log) {
			basis <- cbind(basis, rep(NA, dim(basis)[1]))
		}
	}
	
#add na x
	if (nas) {
		nabasis <- matrix(NA, nrow=length(nax), ncol=nb)
		nabasis[!nax, ] <- basis
		nabasis
	}
	else {
		basis
	}
}




FIntegrateTPBasis<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, ...) { 
	stopifnot(is.numeric(x))
	nax <- is.na(x)
	if (nas <- any(nax)){ 
		x <- x[!nax]
	}
	Aknots<-object@knots
	degree<-object@degree
	degreep1 <- degree+1
	if (xmin != -Inf){
		ol <- x < xmin
		if (xmax != +Inf){
			or <- x > xmax
			outside <- ( or | ol)
		}
		else {
			outside <- ol
		}
	}
	else {
		if (xmax != +Inf){
			outside <- x > xmax
		}
		else {
			outside <- FALSE
		}
	}
	
	Aknotref <- Aknots 
	if( !is.null(ref)){
		x <- x - ref
		if(length(Aknots)>0){
			Aknotref <- Aknots - ref
		}
	}
	
	
	if (intercept) {
		nb <- degree+ 1 + length(Aknots) 
		last <- degree+1
	}
	else{
		nb <- degree  + length(Aknots) 
		last <- degree
	}    
	basis <- matrix(1, nrow=length(x), ncol=nb)
	
	
	
	if (any(!outside)){ 
		if (intercept) {
			for(i in 1:(degree+1)){
				basis[,i]<-x^i/i
			}
		}
		else{
			for(i in 2:(degree+1)){
				basis[,i+1]<-x^i/i
			}
		}
		if(length(Aknots)>0){
			for (k in 1:(length(Aknots))){
				basis[,last+k] <- ifelse(x>Aknotref[k], (x-Aknotref[k])^degreep1/degreep1, 0)
			}
		}
		if (object@log) {
			basis <- cbind(basis, 1/x)
		}
	}    
	if (any(outside)) {
		basis[outside,]<-rep(NA, dim(basis)[2])
		if (object@log) {
			basis <- cbind(basis, rep(NA, dim(basis)[1]))
		}
	}
	
#add na x
	if (nas) {
		nabasis <- matrix(NA, nrow=length(nax), ncol=nb)
		nabasis[!nax, ] <- basis
		nabasis
	}
	else {
		basis
	}
}



setMethod("integrate",signature("BSplineBasis", "missing"),integrate.BSplineBasis)
setMethod("integrate",signature("MSplineBasis", "missing"),integrate.MSplineBasis)
setMethod("integrate",signature("EBSplineBasis", "missing"),integrate.EBSplineBasis)
setMethod("integrate",signature("LEBSplineBasis", "missing"),integrate.LEBSplineBasis)
setMethod("integrate",signature("TPSplineBasis", "missing"),integrate.TPSplineBasis)

setMethod("integrate",signature("BSplineBasis","numeric"),function(object, x, ...)IntegrateBBasis(object=object, x=x, ...))
setMethod("integrate",signature("EBSplineBasis","numeric"),function(object, x, ...)IntegrateEBBasis(object=object, x=x, ...))
setMethod("integrate",signature("LEBSplineBasis","numeric"),function(object, x, ...)IntegrateLEBBasis(object=object, x=x, ...))
setMethod("integrate",signature("TPSplineBasis","numeric"),function(object, x, ...)FIntegrateTPBasis(object=object, x=x, ...))
#setMethod("integrate",signature("TPRSplineBasis","numeric"),function(object, x, ...)IntegrateTPBasis(object=object, x=x, ref=object@ref, ...))
#setMethod("integrate",signature("C0BSplineBasis","numeric"),function(object, x, ...)IntegrateC0SBasis(object=object, x=x, ...))



######################################################################
# deriv method
# define parameters for derived Spline Basis
# MSplinebasis
deriv.BSplineBasis<-function(expr){
	SB <- orthogonalsplinebasis::deriv(expr@SplineBasis)
	
	no <- new("BSplineBasis", knots=SB@knots, min=expr@min, max=expr@max,
			degree=as.integer(SB@order-1), nbases=as.integer(dim(SB@Matrices)[2]), Matrices=SB@Matrices, SplineBasis=SB, log=FALSE)
	
	no
}


# EBSplinebasis
deriv.EBSplineBasis<-function(expr){
	# 
	dnew<-d<-dim(expr)
#  DM<-orthogonalsplinebasis:::DerivativeMatrix(d[1])
	DM <- rbind(0, diag(x = 1:(d[1]-1), nrow = d[1] - 1, ncol = d[1]))
	
	M<-expr@Matrices
	order<-getOrder(order)
	knots<- getKnots(knots)
	nknots<-length(knots)
	
	dnew[1]=d[1]-1L
	N<-array(0,dim=dnew)
#ori	w<-(diff(knots)[k:(n-k)])^-1
	w<-c(1/(diff(knots)[order:(nknots-order)]),1)
	for( i in 1:d[3]) { 
		N[,,i]<-(w[i]*t(DM)%*%M[,,i])[seq(d[1]-1),]
	}
	
	newknots<-knots[seq(2,(nknots-1))]
	neworder <- order-1
	no <- new("EBSplineBasis", knots=newknots, min=expr@min, max=expr@max,
			degree=as.integer(neworder-1), nbases=getNBases(expr), Matrices=N, log=FALSE)
	no
}


# LEMSplinebasis
deriv.LEBSplineBasis<-function(expr){
	SB <- orthogonalsplinebasis::deriv(expr@SplineBasis)
	
	no <- new("LEBSplineBasis", knots=SB@knots, min=expr@min, max=expr@max,
			degree=as.integer(SB@order-1), nbases=as.integer(dim(SB@Matrices)[2]), Matrices=SB@Matrices, SplineBasis=SB,
			orderextrapol = ifelse(expr@orderextrapol > 0, expr@orderextrapol-1L, 0L), log=FALSE)
	
	linexinf <- ((diag(0:expr@degree))%*%expr@linexinf)[-1, , drop=FALSE]
	linexsup <- ((diag(0:expr@degree))%*%expr@linexsup)[-1, , drop=FALSE]
	no@linexinf <- linexinf
	no@linexsup <- linexsup
	no
}


setMethod("deriv",signature("BSplineBasis"),deriv.BSplineBasis)
setMethod("deriv",signature("EBSplineBasis"),deriv.EBSplineBasis)
setMethod("deriv",signature("LEBSplineBasis"),deriv.LEBSplineBasis)






InitCoefSBasis<-function(object,ncol, init=1, intercept=TRUE, xname=NULL, ...) { 
	# knots are all knots (first and last replicated
	# output matrix with all "init" 
	stopifnot(is.integer(ncol))
	nb <- object@nbases + 1 - intercept + object@log
	matrix(init, ncol=ncol , nrow=nb)
	
}

InitCoefBBasis<-function(object,ncol, init=1, intercept=TRUE, xname=NULL, ...) { 
	# output matrix with all "init" 
	# knots are all knots (first and last replicated
	stopifnot(is.integer(ncol))
	nb <- object@nbases + 1 - intercept + object@log
	
	matrix(init, ncol=ncol , nrow=nb)
	
}


InitCoefTPBasis<-function(object,ncol, init=1, intercept=TRUE, xname=NULL, ...) { 
	stopifnot(is.integer(ncol))
	# knots are ordered  interior knots
	# init vectoe to (1, 0, 0, ...)
	nb <- object@nbases + 1 - intercept + object@log
	
	vv <- rep(0,nb)
	vv[1] <- 1
	
	matrix(vv, ncol=ncol , nrow=nb)
	
}

setMethod("initcoef",
		signature("BSplineBasis","integer"),
		function(object, ncol, ...)InitCoefSBasis(object=object, ncol=ncol, ...)
)
setMethod("initcoef",
		signature("BSplineBasis","integer"),
		function(object, ncol, ...)InitCoefBBasis(object=object, ncol=ncol, ...)
)
setMethod("initcoef",
		signature("LEBSplineBasis","integer"),
		function(object, ncol, ...)InitCoefBBasis(object=object, ncol=ncol, ...)
)
setMethod("initcoef",
		signature("TPSplineBasis","integer"),
		function(object, ncol, ...)InitCoefTPBasis(object=object, ncol=ncol, ...)
)


# idem but for constraints parameters
# for B-spline, all coefs are 1
InitCoefCSBasis<-function(object,ncol,intercept=TRUE, xname=NULL, ...) { 
	# knots are all knots (first and last replicated
	stopifnot(is.integer(ncol))
	nb <- length(object@knots) - object@degree - 3 + intercept  + object@log
	
	matrix(1, ncol=ncol , nrow=nb)
}



# for natural spline, all coefs are null for non intercvept term
InitCoefCTPBasis<-function(object,ncol,intercept=TRUE, xname=NULL, ...) {
	# knots are knot_min and interior knots
	stopifnot(is.integer(ncol))
	nb <- object@degree+ length(object@knots) -3 + intercept + object@log
	
	vv <- rep(0,nb)
	matrix(vv, ncol=ncol , nrow=nb)
	
}

setMethod("initcoefC",
		signature("BSplineBasis","numeric"),
		function(object, ncol, ...)InitCoefCSBasis(object=object, ncol=ncol, ...)
)

setMethod("initcoefC",
		signature("TPSplineBasis","numeric"),
		function(object, ncol, ...)InitCoefCTPBasis(object=object, ncol=ncol, ...))


ExpandMCoefSBasis0<-function(object,ncol,coef, intercept=TRUE, xname=NULL, ...) { 
# add a first row of nrow - sum(matcoef) 
	stopifnot(is.integer(ncol))
	nb <- dim(coef)[1]
	
	rbind(intercept - rep(1, nb) %*% coef , coef)
}

ExpandVCoefSBasis0<-function(object,ncol,coef, intercept=TRUE, xname=NULL, ...) { 
# add a first row of nrow - sum(matcoef) 
	stopifnot(is.integer(ncol))
	if(!is.matrix(coef)){
		coef <- matrix(coef, ncol=ncol)
	}
	
	nb <- dim(coef)[1]+1
	
	rbind(intercept - rep(1, nb-1) %*% coef , coef)
}



######################################################################
# operators


# SplineBasis, from package orthogonalsplinebasis
Prod.S.n <- function(e1, e2) {
	d <- dim(e1@Matrices)
	le <- length(e2)
	if(le == 1 ) { # matching dim
		e1@Matrices <- e1@Matrices * e2
		e1
	} else if (le >= d[2]){
# on multiplie chaque Matrices[,,i] par diag(e2)
		for(i in 1:d[3]){
			e1@Matrices[,,i] <- e1@Matrices[,,i] %*%diag(e2[1:d[2]])
		}
		e1
	}
	else stop ("length of args does not")
}
setMethod("*", signature(e1 = "SplineBasis", e2 = "numeric"), Prod.S.n)
Prod.n.S <- function(e1, e2) {
	e2 * e1
}
setMethod("*", signature(e1 ="numeric" , e2 = "SplineBasis"), Prod.n.S)


Sum.S.S <- function(e1, e2) {
	if(e1@knots== e2@knots && e1@order == e2@order){
		e1@Matrices <- e1@Matrices + e2@Matrices
		e1
	}
	else stop("args cannot be added")
}
setMethod("+", signature(e1 = "SplineBasis", e2 = "SplineBasis"), Sum.S.S)


Sum.S.n <- function(e1, e2) {
	e1 + e2 * SplineBasis(knots=e1@knots, order=e1@order, keep.duplicates=TRUE)
}
setMethod("+", signature(e1 = "SplineBasis", e2 = "numeric"), Sum.S.n)

Sum.n.S <- function(e1, e2) {
	e2 + e1
}
setMethod("+", signature(e1 ="numeric" , e2 = "SplineBasis"), Sum.n.S)



Dif.S.S <- function(e1, e2) {
	if(e1@knots== e2@knots && e1@order == e2@order){
		e1@Matrices <- e1@Matrices - e2@Matrices
		e1
	}
	else stop("args cannot be added")
}
setMethod("-", signature(e1 = "SplineBasis", e2 = "SplineBasis"), Dif.S.S)

Dif.S.n <- function(e1, e2) {
	e1 + ((-1) * e2)
}
setMethod("-", signature(e1 = "SplineBasis", e2 = "numeric"), Dif.S.n)


Dif.n.S <- function(e1, e2) {
	( e2 * (-1)) + e1
}
setMethod("-", signature(e1 = "numeric", e2 = "SplineBasis"), Dif.n.S)


# product matrix %*% SplineBasis makes change basis
# if dim(matrix)[1] < nbases, we get a projection in the subspace
# by extension/convention SplineBasis %*% matrix too


Prod.m.S <- function(x, y) {
	d1 <- dim(x)
	d2 <- dim(y)
	if(d1[2] == d2[2] ) { # matching dim
		newmatrices <- array(0, dim=c(d2[1], d1[1],d2[3])) 
		for( i in 1:d2[3]){
			# because of the order of dims in Matrices
			newmatrices[,,i] <- y@Matrices[,,i] %*% t(x)
		}
		y@Matrices <- newmatrices
	} 
	else stop ("dim of args do not match")
	return(y)
}
setMethod("%*%", signature(x ="matrix" , y = "SplineBasis"), Prod.m.S)


Prod.S.m <- function(x, y) {
	y %*% x
}
setMethod("%*%", signature(x = "SplineBasis", y = "matrix"), Prod.S.m)


######################################################################
# MSplineBasis
Prod.MS.n <- function(e1, e2) {
	e1@SplineBasis <- e1@SplineBasis * e2
	e1@Matrices <- e1@SplineBasis@Matrices
	e1
}
setMethod("*", signature(e1 = "BSplineBasis", e2 = "numeric"), Prod.MS.n)

Prod.n.MS <- function(e1, e2) {
	e2 * e1
}
setMethod("*", signature(e1 ="numeric" , e2 = "BSplineBasis"), Prod.n.MS)


Sum.MS.MS <- function(e1, e2) {
	e1@SplineBasis <- e1@SplineBasis + e2@SplineBasis 
	e1@Matrices <- e1@SplineBasis@Matrices
	e1
}
setMethod("+", signature(e1 = "BSplineBasis", e2 = "BSplineBasis"), Sum.MS.MS)


Sum.MS.n <- function(e1, e2) {
	e1@SplineBasis <- e1@SplineBasis + e2 
	e1@Matrices <- e1@SplineBasis@Matrices
	e1
}
setMethod("+", signature(e1 = "BSplineBasis", e2 = "numeric"), Sum.MS.n)

Sum.n.MS <- function(e1, e2) {
	e2 + e1
}
setMethod("+", signature(e1 ="numeric" , e2 = "BSplineBasis"), Sum.n.MS)



Dif.MS.MS <- function(e1, e2) {
	e1@SplineBasis <- e1@SplineBasis - e2@SplineBasis 
	e1@Matrices <- e1@SplineBasis@Matrices
	e1
}
setMethod("-", signature(e1 = "BSplineBasis", e2 = "BSplineBasis"), Dif.MS.MS)

Dif.MS.n <- function(e1, e2) {
	e1 + ((-1) * e2)
}
setMethod("-", signature(e1 = "BSplineBasis", e2 = "numeric"), Dif.MS.n)

Dif.n.MS <- function(e1, e2) {
	( e2 * (-1)) + e1
}
setMethod("-", signature(e1 = "numeric", e2 = "BSplineBasis"), Dif.n.MS)

# product matrix %*% SplineBasis makes change basis
# if dim(matrix)[1] < nbases, we get a projection in the subspace
# by extension/convention SplineBasis %*% matrix too


Prod.m.MS <- function(x, y) {
	y@SplineBasis <- x %*% y@SplineBasis
	y@Matrices <- y@SplineBasis@Matrices
	y@nbases <- dim(y@SplineBasis)[2]
	return(y)
}
setMethod("%*%", signature(x ="matrix" , y = "BSplineBasis"), Prod.m.MS)


Prod.MS.m <- function(x, y) {
	y %*% x
}
setMethod("%*%", signature(x = "BSplineBasis", y = "matrix"), Prod.MS.m)



######################################################################
# LEBSplineBasis
Prod.LEMS.n <- function(e1, e2) {
	e1@SplineBasis <- e1@SplineBasis * e2
	e1@Matrices <- e1@SplineBasis@Matrices
	e1@linexinf <- e1@linexinf %*% diag(e2, ncol=e1@nbases, nrow=e1@nbases)
	e1@linexsup <- e1@linexsup %*% diag(e2, ncol=e1@nbases, nrow=e1@nbases)
	e1
}
setMethod("*", signature(e1 = "LEBSplineBasis", e2 = "numeric"), Prod.LEMS.n)

Prod.n.LEMS <- function(e1, e2) {
	e2 * e1
}
setMethod("*", signature(e1 ="numeric" , e2 = "LEBSplineBasis"), Prod.n.LEMS)


Sum.LEMS.LEMS <- function(e1, e2) {
	e1@SplineBasis <- e1@SplineBasis + e2@SplineBasis 
	e1@Matrices <- e1@SplineBasis@Matrices
	e1@linexinf <- e1@linexinf + e2@linexinf 
	e1@linexsup <- e1@linexsup + e2@linexsup 
	
	e1
}
setMethod("+", signature(e1 = "LEBSplineBasis", e2 = "LEBSplineBasis"), Sum.LEMS.LEMS)


Sum.LEMS.n <- function(e1, e2) {
	e1 + ( LEBSplineBasis2(allknots=e1@knots, degree=e1@degree, keep.duplicates=FALSE) * e2)
}
setMethod("+", signature(e1 = "LEBSplineBasis", e2 = "numeric"), Sum.LEMS.n)

Sum.n.LEMS <- function(e1, e2) {
	e2 + e1
}
setMethod("+", signature(e1 ="numeric" , e2 = "LEBSplineBasis"), Sum.n.LEMS)



Dif.LEMS.LEMS <- function(e1, e2) {
	e1@SplineBasis <- e1@SplineBasis - e2@SplineBasis 
	e1@Matrices <- e1@SplineBasis@Matrices
	e1@linexinf <- e1@linexinf - e2@linexinf 
	e1@linexsup <- e1@linexsup - e2@linexsup 
	e1
}
setMethod("-", signature(e1 = "LEBSplineBasis", e2 = "LEBSplineBasis"), Dif.LEMS.LEMS)

Dif.LEMS.n <- function(e1, e2) {
	e1 + ((-1) * e2)
}
setMethod("-", signature(e1 = "LEBSplineBasis", e2 = "numeric"), Dif.LEMS.n)

Dif.n.LEMS <- function(e1, e2) {
	( e2 * (-1)) + e1
}
setMethod("-", signature(e1 = "numeric", e2 = "LEBSplineBasis"), Dif.n.LEMS)


# product matrix %*% SplineBasis makes change basis
# if dim(matrix)[1] < nbases, we get a projection in the subspace
# by extension/convention SplineBasis %*% matrix too


Prod.m.LEMS <- function(x, y) {
	y@SplineBasis <- x %*% y@SplineBasis
	y@Matrices <- y@SplineBasis@Matrices
	y@linexinf <- y@linexinf %*% t(x)
	y@linexsup <- y@linexsup %*% t(x)
	y@nbases <- dim(y@SplineBasis)[2]
	return(y)
}
setMethod("%*%", signature(x ="matrix" , y = "LEBSplineBasis"), Prod.m.LEMS)


Prod.LEMS.m <- function(x, y) {
	y %*% x
}
setMethod("%*%", signature(x = "LEBSplineBasis", y = "matrix"), Prod.LEMS.m)


######################################################################
# TPSplineBasis
Prod.TPS.n <- function(e1, e2) {
	d <- length(e1@coef)
	le <- length(e2)
	if(le == 1 ) { # matching dim
		e1@coef <- e1@coef * e2
		e1
	} else if (le >= d){
		e1@coef <- (e1@coef * e2[1:d])
		e1
	}
	else stop ("length of args does not")
}
setMethod("*", signature(e1 = "TPSplineBasis", e2 = "numeric"), Prod.TPS.n)

Prod.n.TPS <- function(e1, e2) {
	e2 * e1
}
setMethod("*", signature(e1 ="numeric" , e2 = "TPSplineBasis"), Prod.n.TPS)


Sum.TPS.TPS <- function(e1, e2) {
	if(e1@knots== e2@knots && e1@degree == e2@degree && e1@min == e2@min && e1@max == e2@max && e1@degrees == e2@degrees && e1@type == e2@type){
		e1@coef <- e1@coef + e1@coef 
		e1
	}
	else stop("args cannot be added")
}
setMethod("+", signature(e1 = "TPSplineBasis", e2 = "TPSplineBasis"), Sum.TPS.TPS)


Sum.TPS.n <- function(e1, e2) {
	d <- length(e1@coef)
	le <- length(e2)
	if(le == 1 ) { # matching dim
		e1@coef[1] <- e1@coef[1]  +  e2
		e1
	} else if (le >= d){
		e1@coef <- e1@coef + e2[1:d]
		e1
	}
	else stop ("length of args does not match")
}
setMethod("+", signature(e1 = "TPSplineBasis", e2 = "numeric"), Sum.TPS.n)


Sum.n.TPS <- function(e1, e2) {
	e2 + e1
}
setMethod("+", signature(e1 ="numeric" , e2 = "TPSplineBasis"), Sum.n.TPS)


Dif.TPS.TPS <- function(e1, e2) {
	if(e1@knots== e2@knots && e1@degree == e2@degree && e1@min == e2@min && e1@max == e2@max && e1@degrees == e2@degrees && e1@type == e2@type){
		e1@coef <- e1@coef - e1@coef 
		e1
	}
	else stop("args cannot be added")
}
setMethod("-", signature(e1 = "SplineBasis", e2 = "SplineBasis"), Dif.S.S)

Dif.TPS.n <- function(e1, e2) {
	e1 + ((-1) * e2)
}
setMethod("-", signature(e1 = "TPSplineBasis", e2 = "numeric"), Dif.TPS.n)

Dif.n.TPS <- function(e1, e2) {
	( e2 * (-1)) + e1
}
setMethod("-", signature(e1 = "numeric", e2 = "TPSplineBasis"), Dif.n.TPS)



##################################################
#  Marsden's Identity

MarsdenIdentity <- function(object, order = 1){
	Knots.marsden <- getKnots(object)
	degree <- getDegree(object)
	NBasis <- getNBases(object)
	marsden.matrix <- matrix(0, ncol=length(Knots.marsden), nrow=NBasis)
	for (j in 1:NBasis) {
		marsden.matrix[j, j+(1:degree)] <- 1.0
	}
	return(marsden.matrix%*% Knots.marsden/degree) 
	
}

