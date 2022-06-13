SumDirect <- function(m1, m2){
	m1 <-as.matrix(m1)
	m2 <-as.matrix(m2)
	
	topright <- matrix(0.0, nrow(m1), ncol(m2))
#		colnames(topright) <- colnames(m2)
	lowleft <- matrix(0.0, nrow(m2), ncol(m1))
	return( rbind(cbind(m1, topright), cbind(lowleft, m2)))
}

"%sd%" <- function( x, y )
{
	return( SumDirect( x, y ) )
}


duplicMat <- function(m, n){
	dm <-NULL
	for(i in 1:n){
		dm <- cbind(dm, m)
	}
	dm
}

duplicSumDirect <- function(m, n){
	m <- as.matrix(m)
	dm <-m
	if(n> 1){
		for(i in 2:n){
			dm <- dm %sd% m
		}
	}
	dm
}


StackDiag <- function(vect, n){
	dm <-NULL
	for(i in 1:length(vect)){
		dm <- rbind(dm, diag(vect[i],n))
	}
	dm
}

