getIndexId <- function(Id){
	uid <- unique(Id)
	fId <- factor(Id, levels=as.character(uid))
	nline <- tabulate(fId)
	Nline <- cumsum(nline)
	LastId <- rep(Nline, nline)
	FirstId <- LastId - rep(nline, nline) + 1
	return(list(FirstId=FirstId, LastId=LastId))
}

getLastId <- function(FirstId){
	first <- unique(FirstId)
	nline <- c(first[-1],length(FirstId)+1)-first
	LastId <- FirstId+rep(nline, nline)-1
	
	return(LastId)
}

getFirstId <- function(LastId){
	last <- unique(LastId)
	nline <- last -  c(0, last[-length(LastId)])
	FirstId <- LastId - rep(nline, nline) + 1
	
	return(FirstId)
}

