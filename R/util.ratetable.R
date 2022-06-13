# permute the dimensions of a rate table
aperm.ratetable <- function(a, perm = NULL, ...){
	na <- aperm.default(a, perm = perm)
	if(is.null(perm)){
		perm <- length(dim(a)):1
	}
	attr(na, "cutpoints") <- attr(a, "cutpoints")[perm]
	dimid <- attr(a, "dimid")
	if (!is.null( attr(a, "dimid"))){
		attr(na, "dimid") <- attr(a, "dimid")[perm]
	}
	attr(na, "type") <- attr(a, "type")[perm]
	attr(na, "factor") <- attr(a, "factor")[perm]
	class(na) <- class(a)
	return(na)
}


# extract one matrix of rates, dimnames=c(age, year)
extract.ratetable <- function(ratetable, which=NULL, age="age", year="year", ...){
	Call <- match.call(expand.dots = TRUE)
	mrt <- match(c("which", "age", "year"), names(Call), 0L)
	if(sum(mrt)==0){
		# extract.array
		return(NextMethod("extract", ratetable, ...) )
	}
	# which : named dataframe with teh demographic variables of one subject 
	dimid <- attr(ratetable, "dimid")
	if (is.null(dimid)) dimid <- names(dimnames(ratetable))
	iage <- match(age, dimid)
	iyear <- match(year, dimid)
	
	if(is.null(which)){
		which <- as.data.frame(matrix(rep(1, length(dimid)-2), nrow=1))
		names(which) <- dimid[-c(iage, iyear)]
	}
	allnames <- unique(c(names(which),age, year))
	allnames <- names(which)
	onames <- match(dimid[-c(iage, iyear)], allnames, nomatch = 0)
	
	xtr <- as.matrix(expand.grid(c(list(age=1:dim(ratetable)[iage], year=1:dim(ratetable)[iyear]), 
							lapply(as.list(which[,onames]), unclass))))
	mr <- aperm(ratetable, perm=c(iage, iyear, (1:length(dim(ratetable)))[-c(iage, iyear)]))[xtr]
	
	attributes(mr) <- list(dim = dim(ratetable)[c(iage, iyear)],
			dimnames = dimnames(ratetable)[c(iage, iyear)],
			cutpoints = attr(ratetable, "cutpoints")[c(iage, iyear)],
			dimid = dimid[c(iage, iyear)],
			type = attr(ratetable, "type")[c(iage, iyear)],
			factor = attr(ratetable, "factor")[c(iage, iyear)],
			class=class(ratetable))
	return(mr)
}

