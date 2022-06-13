# return the rateend, cumrateend et cumrateenter

getHazardFromTable <- function(Y, startdate, startage, matchdata=NULL,
		ratetable=survival::survexp.us, 
		age=1, year=2,
		rmap,
		agemin=16, scale=365.25,
		ratename = "rateend", cumrateendname ="cumrateend", cumrateentername ="cumrateenter", 
		origin = "01/01/1970", format="%d/%m/%Y",
		left.open = FALSE,
		verbose=FALSE
){
# Y object of class Surv with start, end, status columns
# startage, startdate : age et date du T=0
	
	if (!inherits(Y, "Surv")) {
		stop("Response must be a survival object")
	}
	Survtype <- attr(Y, "type")
	
	if ((ncol(Y) ==  2) ) {
		if (Survtype != "right"){
			stop(gettextf("flexrsurv does not support %s type of censoring with (0, end] survival data", dQuote(Survtype), domaine=NA))
		} else {
			Y <- cbind(rep(0, dim(Y)[1]), Y)
		}
	} else if ((ncol(Y) ==  3) && (Survtype != "counting") ) {
		stop(gettextf("flexrsurv does not support %s type of censoring with (start, end] survival data", dQuote(Survtype), domaine=NA))
	}
	
	if (is.ratetable(ratetable)) {
		dimid <- attr(ratetable, "dimid")
		if (is.null(dimid)) dimid <- names(dimnames(ratetable))
	}
	else stop("Invalid rate table")
	iagert <- match(age, dimid)
	iyearrt <- match(year, dimid)
	
	if (!inherits(startdate, "Date")) {
		startdate <- as.Date(startdate, origin = origin, format=format)
	}
	
	
	minage <- rep(agemin * scale, length(startage) )
	mindate <- as.Date(startdate - (startage - minage), origin = origin, format=format)
	
	#extract variables from matchdata in ratetable according to rmap
	if (missing(rmap) ) {
		rmap <- NULL
	}
	if (!is.null(substitute(rmap)) ) {
		rcall0 <- substitute(rmap)		
		if(!is.call(rcall0)){
			if(is.call(rmap)){
				rcall0 <- rmap		  
			}
			else {
				stop("Invalid rmap argument")
			}
		}
		
		if (!is.call(rcall0) || rcall0[[1]] != as.name("list")){
			stop("Invalid rmap argument")
		}
		
		namesinrt <- names(as.list(rcall0)[-1])
		namesindata <- as.character(as.list(rcall0)[-1])
	}
	else {
		rcall0 <- NULL
		namesinrt <- NULL
		namesindata <- NULL
	}
	
	temp <- match(namesinrt, dimid)
	
	if (any(is.na(temp))) {
		stop("Variable not found in the ratetable:", (names(rcall))[is.na(temp)])
	}
	
	# variable in ratetable that are unlisted are assumed to be in the data with the same name 
	if (any(!(dimid %in% namesinrt))) {
		to.add <- dimid[!(dimid[-c(iagert, iyearrt)] %in% namesinrt)]
		if( !all(to.add %in% names(matchdata))) {
			stop("Some variables in the ratetable not found in the data:", to.add)
		}
		temp1 <- paste(text = paste(to.add, to.add, sep = "="), collapse = ",")
		if (is.null(rcall0)){ 
			rcall <- parse(text = paste("list(", temp1, ")"))[[1]]
			allmapvarrt <- to.add
			allmapvardata <- to.add 
		}
		else {
			temp2 <- deparse(rcall0)
			rcall <- parse(text = paste("c(", temp2, ",list(", temp1, "))"))[[1]]
			namesinrt <- c(namesinrt, to.add)
			namesindata <- c(namesindata, to.add) 
		}
	}
	
	# variables in matchdata in the same order than in dimid = namesinrt[onamesinrt]	
	onamesinrt <- match(namesinrt, dimid[-c(iagert, iyearrt)]) 
	orderednamesinrt <- namesinrt[onamesinrt]
	orderednamesindata <- namesindata[onamesinrt]
	
	newvar <- all.vars(rcall)
	#dim number of matched var in diid 
	imapmdrt <- match(orderednamesinrt, dimid)
	matchdata2 <- matchdata[,orderednamesindata, drop=FALSE]
	tmap2 <- paste(dimid[iagert], "= minage, ", dimid[iyearrt], " = mindate", sep="")
	
	if(!is.null(matchdata)){
		if (!missing(rcall0)) {
			tmap <- deparse(rcall0)
			rcall2 <- parse(text = paste( strsplit(tmap, ")"), ", ", tmap2, ")", sep=""))[[1]]
			if (!is.call(rcall2) || rcall2[[1]] != as.name("list")){
				stop("Invalid rmap argument")
			}
		}
		else {
			stop("rmap must be given")
		}
	}
	else {
		rcall2 <- parse(text = paste("list(", tmap2, ")"))[[1]]
	}
	temp <- match(names(rcall2)[-1], dimid)
	
	
	tmpdata <- data.frame(startage =startage,
			startdate=startdate,
			dstartage=startage  + Y[,1] - minage,
			dendage = startage + Y[,2] - minage,
			minage=minage, mindate=mindate)
#  tmpdata <- data.frame(date=startdate+Y[,1], age=startage+Y[,1])
	if(!is.null(matchdata)){
		tmpdata <- cbind(tmpdata, matchdata2) 
	}
	
	
	
	# compute cumulative hazard at start age/date
#   at enter Y[,1]
	callsurvexpenter <- as.call(list(survexp, formula = dstartage  ~ 1 , data=tmpdata, rmap=rcall2,
					method="individual.h",
					ratetable=ratetable))

	cumratestart <- eval(callsurvexpenter)
	
	
# cumratestart <-  cumratestart0 - cumrateminage
	
	
#   at enter Y[,2]
#  tmpdata <- data.frame(date=startdate+Y[,2], age=startage+Y[,2])
#  if(!is.null(matchdata)){
#    tmpdata <- cbind(tmpdata, matchdata) 
#  }
	
	ageend <- startage+Y[,2]
	dateend <- startdate+Y[,2]
	
#  tmpdata$ageend <- ageend
#  tmpdata$startdate <- dateend
	# compute cumulative hazard at end date
	
	callsurvexpend <- as.call(list(survexp, formula = dendage  ~ 1 , data=tmpdata, rmap=rcall2,
					method="individual.h",
					ratetable=ratetable))
	
	cumrateend <- eval(callsurvexpend)

#  cumrateend <- cumrateend0 - cumrateminage
	
# le taux à la sortie
	n <- length(startdate)
	rateend <- rep(0, n)
	if(verbose){
		lengthmessage <- floor(log10(n+1))+1
		backwd <- paste(rep('\b', lengthmessage), sep='', collapse='')
		cat(paste("compute cumulative hazard:", format(0, width=lengthmessage, trim=FALSE), sep=""))
	}
	
	
	if(!is.null(matchdata)){
		names(matchdata2) <- orderednamesinrt
	}
	for( i in 1:n) {
		if(verbose){
			cat(paste(backwd, format(i, width=lengthmessage, trim=FALSE), sep=""))
		}
		if(is.null(matchdata)){
			md <- NULL
		}
		else {
			md <- matchdata2[i,, drop=FALSE]
			nmd <- names(matchdata2)
		}
		
		# get the rate matrix, with dim1=age, dim2 = peroide 
	
		if(!inherits(ratetable, "array")){
			# in order to use extract.array 
			class(ratetable) <- c(class(ratetable), "array")
		} 
		
		ratematrix <- extract.ratetable(ratetable, which=md, age=age, year=year)
		cutdate <- c(attr(ratetable, "cutpoints")[[iyearrt]], +Inf)
		cutage  <- c(attr(ratetable, "cutpoints")[[iagert]], +Inf)
		cutdate <- c(attr(ratetable, "cutpoints")[[iyearrt]])
		cutage  <- c(attr(ratetable, "cutpoints")[[iagert]])
		cutpoints <- list(cutage=cutage, cutdate=cutdate)
		# get cell index of endage, enddate
		iendage  <- findInterval(ageend[i],  cutage, left.open=left.open)
		
		if(attr(ratetable, "type")[iyearrt] ==  3){
			# for ageXperiode rate table
			ienddate <- findInterval(dateend[i], cutdate, left.open=left.open)
		} else if(attr(ratetable, "type")[iyearrt] ==  4){
			# for us.special rate table (ie pseudo ageXgeneration rate table)
			# find the date at cutage[iendage]
			ienddate <- findInterval(dateend[i]-(ageend[i]-cutage[iendage]), cutdate, left.open=left.open)
#      ienddate <- findInterval(startdate[i]-(ageend[i]-cutage[iendage]), cutdate, left.open=left.open)
		} else {
			stop("wrong type of rate table, attr(ratetable, 'type')} sould be 3 or 4.")
		}	
		
		cell <- c(iendage, ienddate)
		rateend[i] <- ratematrix[cell[1], cell[2]]
		
	}
	
	if(verbose){
		cat("\n")
	}
	
	ret <- as.data.frame(cbind(rateend, cumratestart, cumrateend))
	names(ret) <- c( ratename, cumrateentername , cumrateendname )
	
	return(ret)
	
	
	
	
}
