# computation in a age-periode table of individuals's indices, cut by periode, 
# pseudo cumulative rates are computed staying in the same period (ie distribution of ages in a period)
# returned value :
# vector rateend, cumrateend et cumrateenter
# vector FirstHId, LastHId such that **rate**[FirstHId[i], LastHId[i]] corresponds to the 
# **rate** of row i in Y
getBrassPseudoHazardFromTable <- function(Y, startdate, startage, matchdata=NULL,
		ratetable=survival::survexp.us, 
		age=1, year=2,
		rmap,
		agemin=16, scale=365.25,
		ratename = "rateend", cumrateendname ="cumrateend", cumrateentername ="cumrateenter", idname = "Id_byperiod", 
		origin = "01/01/1970", format="%d/%m/%Y",
		left.open = FALSE,
		SplineBrass=R2bBSplineBasis(knots=c(-2.5,0,2.5), degree=3) * c(1.0, 0.0, 1.0),
		verbose=FALSE
){
# Y : object of class Surv
# startdate, startage, age and date when Y[,] == 0
	# the beginning age and date of the row i is startage[i] + Y[i,1] startdate[i] + Y[i,1]
	# the exit      age and date of the row i is startage[i] + Y[i,2] startdate[i] + Y[i,2]
# ratetable object of class ratetable 
# age=1, year=2, name or index of the age and periode variable in the rate table
# rmap,          list for mappind variables in the ratetable and the data : list(var1rt = var1data, var2rt = var2data)
# agemin=16, scale=365.25,    age at which starts the cumulative hazards
# ratename, cumrateendname , cumrateentername, idname , names of the variable created
# left.open = FALSE,  passed to findInterval
	# object (splinBasis) used to transform the logit of the survival
# verbose=FALSE       print intermediate messages
	
	
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
		if(!inherits(ratetable, "array")){
			# in order to use extraxt.array 
			class(ratetable) <- c(class(ratetable), "array")
		} 
	}
	else stop("Invalid rate table")
	iagert <- match(age, dimid)
	iyearrt <- match(year, dimid)
	
#	if (!inherits(startdate, "Date")) {
#		startdate <- as.Date(startdate, origin = origin, format=format)
#	}
	
	class.startdate <- class(startdate)
	# from now on, all dates are with class Date
	startdate <- as.Date(as.numeric(ratetableDate(startdate)))
	
	minage <- rep(agemin * scale, length(startage) )
	mindate <- startdate - (startage - minage)
	
	#extract variables from matchdata in ratetable according to rmap
	if (!missing(rmap)) {
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
		rmap <- NULL
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
	# here, the pseudo cumulative rate is computed at the period pointed bys the starting date in the split data
	tmap2 <- paste(dimid[iagert], "= minage, ", dimid[iyearrt], " = mindate", sep="")
	tmap3 <- paste(dimid[iagert], "= tagestop1, ", dimid[iyearrt], " = tdatestop1", sep="")
	
	if(!is.null(matchdata)){
		if (!missing(rmap)) {
			tmap <- deparse(rcall0)
			rcall2 <- parse(text = paste( strsplit(tmap, ")"), ", ", tmap2, ")", sep=""))[[1]]
			rcall3 <- parse(text = paste( strsplit(tmap, ")"), ", ", tmap3, ")", sep=""))[[1]]
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
		rcall3 <- parse(text = paste("list(", tmap3, ")"))[[1]]
	}
	temp <- match(names(rcall2)[-1], dimid)
	
	# dates in tmpdata are of class Date,  
	tmpdata <- data.frame(startage = startage,
			startdate= startdate,
			Tstartdate= startdate + Y[,1],
#						Tenddate =  as.Date(startdate, origin = origin, format=format) + Y[,2],
			Tenddate =  startdate + Y[,2],
			minage=minage, mindate=mindate)
# matching data
	if(!is.null(matchdata)){
		tmpdata <- cbind(tmpdata, matchdata2) 
	}
	
	# cutpoints in the ratetable, with dates converted to class Date
	
	cutdate <- ratetableDate(attr(ratetable, "cutpoints")[[iyearrt]])
	cutage  <- c(attr(ratetable, "cutpoints")[[iagert]])
	cutpoints <- list(cutage=cutage, cutdate=cutdate)
	
	# split data at the period/year cutpoints of the rate table
  # survSplit need numerics in Surv() and cut
	splitdata <- survSplit(Surv(as.numeric(Tstartdate), as.numeric(Tenddate), rep(0,dim(tmpdata)[1])) ~., data=tmpdata,
			cut=as.numeric(cutdate), id= idname, start="tdatestart", end="tdatestop", event="event")
	
	splitdata$tdatestart2 <- as.Date(splitdata$tdatestart)
	splitdata$tdatestop2 <- as.Date(splitdata$tdatestop)

	splitdata$dtdatestart <- as.Date(splitdata$tdatestart) - splitdata$mindate
	splitdata$dtdatestop <- as.Date(splitdata$tdatestop) - splitdata$mindate
	
	splitdata$tdatestop1 <- splitdata$mindate + splitdata$dtdatestop -1
	splitdata$tagestop1 <- splitdata$minage + splitdata$dtdatestop -1
	splitdata$oneday <- 1
	

callsurvexpenter <- as.call(list(survexp, formula = dtdatestart  ~ 1 , data=splitdata, rmap=rcall2,
				method="individual.h",
				ratetable=ratetable))

cumratestart <- eval(callsurvexpenter)

callsurvexpend <- as.call(list(survexp, formula = dtdatestop  ~ 1 , data=splitdata, rmap=rcall2,
				method="individual.h",
				ratetable=ratetable))

cumrateend <- eval(callsurvexpend)


callsurvexprate <- as.call(list(survexp, formula = oneday  ~ 1 , data=splitdata, rmap=rcall3,
				method="individual.h",
				ratetable=ratetable))

rateend <- eval(callsurvexprate)



	splitdata$istartperiod <- findInterval(splitdata$tdatestart, cutdate, left.open = FALSE)
	splitdata$iendperiod <- findInterval(splitdata$tdatestop, cutdate, left.open = TRUE)
	
	if(any(splitdata$iendperiod != splitdata$istartperiod)){
		warning("rows lines with splitdata$iendperiod != splitdata$istartperiod; splitdata$istartperiod used")
	}
	
# time since minage ie mindate (unused)
	splitdata$starttime <- splitdata$tdatestart2- splitdata$mindate
	splitdata$endtime <- splitdata$tdatestop2- splitdata$mindate
	
# age at starttime and endtime = date@start - birthdate
	splitdata$tstartage <- splitdata$tdatestart2 - splitdata$startdate + splitdata$startage
	splitdata$tendage <- splitdata$tdatestop2 - splitdata$startdate + splitdata$startage
	
	iminage <- findInterval(agemin * scale, cutage, left.open = FALSE)
	
	splitdata$istartage <- findInterval(splitdata$tstartage, cutage, left.open = FALSE)
	splitdata$iendage <- findInterval(splitdata$tendage, cutage, left.open = TRUE)
	#  print("+3")
	
# match data
	md <- data.frame(eval(rcall, splitdata), stringsAsFactors = TRUE)
	
	
	ret <- data.frame(splitdata[[idname]], rateend, cumratestart, cumrateend, rep(iminage, length(cumratestart)), 
			splitdata$istartage, splitdata$iendage, splitdata$mindate, splitdata$minage, 
			splitdata$tdatestart2, splitdata$tdatestop2, 
			splitdata$tstartage, splitdata$tendage )
	names(ret) <- c( idname, ratename, cumrateentername , cumrateendname, "iminage",
			"istartage", "iendage", "mindate", "minage", 
			"tdatestart", "tdatestop", 
			"tstartage", "tendage" )
	
# brass transformed values
	
	send <- exp(-ret[[cumrateendname]]) 
	logitend <- log((1- send)/send)
	sstart <- exp(-ret[[cumrateentername]]) 
	logitstart <- log((1- sstart)/sstart)
	ret[[cumrateentername]] <- log(1 + exp(predictSpline(SplineBrass, logitstart)))
	ret[[cumrateendname]] <- log(1 + exp(predictSpline(SplineBrass, logitend)))
	ret[[ratename]] <- ret[[ratename]] * (1 + exp(-logitend))/(1 + exp(-predictSpline(SplineBrass, logitend))) * 
			predictSpline(deriv(SplineBrass), logitend)
	
	
	
# extract values at entry end exit of Y
	uid <- sort(unique(ret[[idname]]))
	fId <- factor(ret[[idname]], levels=as.character(uid))
	nline <- tabulate(fId)
	LastId <- cumsum(nline)
	FirstId <- c(1, LastId[-length(LastId)]+1)
	
	attr(ret, "entryexit") <- data.frame(ratename=ret[[ratename]][LastId],
			cumrateentername=ret[[cumrateentername]][FirstId], 
			cumrateendname=ret[[cumrateendname]][LastId])
	names(attr(ret, "entryexit")) <- c( ratename, cumrateentername , cumrateendname)
	attr(ret, "cutpoint") <-cutpoints
	return(ret)
	
}

