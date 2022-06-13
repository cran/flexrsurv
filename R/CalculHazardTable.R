# retourne le rateend, cumrateend et cumrateenter

calculhazardtable <- function(Y, startdate, startage, matchdata=NULL,  ratetable=survival::survexp.us, type.ratetable=c("ageperiod", "pseudogeneration"),
		left.open = FALSE,
		age=1, year=2,
		rmap,
		ageminbrass=16, scale=365.25,
		timename="time", ratename = "rateend", cumrateendname ="cumrateend", cumrateentername ="cumrateenter", 
		make.data.frame=TRUE, Id=1:n,
		epsilon=1e-6, niter=20,
		...
){
# Y object of class Surv with start, end, status columns
# startage, startdate : age et date du T=0
	
	type.ratetable <- match.arg(type.ratetable) 
	
	if (is.ratetable(ratetable)) {
		dimid <- attr(ratetable, "dimid")
		if (is.null(dimid)) dimid <- names(dimnames(ratetable))
	}
	else stop("Invalid rate table")
	
	iage <- match(age, dimid)
	iyear <- match(year, dimid)
	
	# dans cette version, age et year sont les numéros des variables
#iage <- dimid[age]
#  iyear <- dimid[year]
	
#  print(dimid)
#print(c(iage, iyear))
	
	
	minage <- rep(ageminbrass * scale, length(startage) )
	mindate <- startdate - (startage - minage)
	
	
	tmap2 <- paste(dimid[iage], "= minage, ", dimid[iyear], " = mindate", sep="")
#  print("tmap2")
#  print(tmap2)
#  rmap2 <- list(as.name(startage), as.name(startdate)  )
#  names(rmap2) <- dimid[c(iage, iyear)]
#  rcall2<-substitute(rmap2)
	
	if(!is.null(matchdata)){
		if (!missing(rmap)) {
			rcall <- substitute(rmap)
			tmap <- deparse(rcall)
			rcall2 <- parse(text = paste( strsplit(tmap, ")"), ", ", tmap2, ")", sep=""))[[1]]
			if (!is.call(rcall) || rcall[[1]] != as.name("list")) 
				stop("Invalid rcall argument")
		}
		else {
			stop("rmap must be given")
		}
	}
	else {
		rcall2 <- parse(text = paste("list(", tmap2, ")"))[[1]]
	}
	
	temp <- match(names(rcall2)[-1], dimid)
	
	
#  print(rcall2)
	
	
	tmpdata <- data.frame(startdate=startdate,
			dstartage=startage  + Y[,1] - minage,
			dendage = startage + Y[,2] - minage,
			minage=minage, mindate=mindate)
#  tmpdata <- data.frame(date=startdate+Y[,1], age=startage+Y[,1])
	if(!is.null(matchdata)){
		tmpdata <- cbind(tmpdata, matchdata) 
	}
	
	# compute cumulative hazard at ageminbrass car on calcul les cumul à partir de agemin
	
#  callsurvexpmin <- as.call(list(survexp, formula = startage  ~ 1 , data=tmpdata, rmap=rcall2,
#                                      method="individual.h",
#                                      ratetable=ratetable))
#  cumrateminage <- eval(callsurvexpmin)
	
	
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
	callsurvexpend <- as.call(list(survexp, formula = dendage  ~ 1 , data=tmpdata, rmap=rcall2,
					method="individual.h",
					ratetable=ratetable))
	cumrateend <- eval(callsurvexpend)
	
#  cumrateend <- cumrateend0 - cumrateminage
	
# le taux à la sortie
	n <- length(startdate)
	rateend <- rep(0, n)
	
	lengthmessage <- floor(log10(n+1))+1
	backwd <- paste(rep("\b", lengthmessage), sep="", collapse="")
	cat(paste("compute cumulative hazard:", format(0, width=lengthmessage, trim=FALSE), sep="")) 
	for( i in 1:n) {
#    cat(paste(backwd, format(i, width=lengthmessage, trim=FALSE), sep=""))
		if(is.null(matchdata)){
			md <- NULL
		}
		else {
			md <- matchdata[i,, drop=FALSE]
			nmd <- names(matchdata)
		}
		# get the rate matrix, with dim1=age, dim2 = peroide 
		ratematrix <- extract.ratetable(ratetable, which=md, age=age, year=year)
		cutdate <- c(attr(ratetable, "cutpoints")[[iyear]], +Inf)
		cutage  <- c(attr(ratetable, "cutpoints")[[iage]], +Inf)
		cutdate <- c(attr(ratetable, "cutpoints")[[iyear]])
		cutage  <- c(attr(ratetable, "cutpoints")[[iage]])
		cutpoints <- list(cutage=cutage, cutdate=cutdate)
		# get cell index of endage, enddate
		iendage  <- findInterval(ageend[i],  cutage, left.open=left.open)
		
		if(type.ratetable == "ageperiod"){
			# for ageXperiode rate table
			ienddate <- findInterval(dateend[i], cutdate, left.open=left.open)
		} else if(type.ratetable == "pseudogeneration"){
			# for pseudo ageXgeneration rate table
			# find the date at cutage[iendage]
			ienddate <- findInterval(dateend[i]-(ageend[i]-cutage[iendage]), cutdate, left.open=left.open)
		} else if(type.ratetable == "generation"){
			stop("type.ratetable == 'generation' not yet implemented")
		}
		cell <- c(iendage, ienddate)
		rateend[i] <- ratematrix[cell[1], cell[2]]
		
	}
	
	cat("\n")
	
	
	ret <- as.data.frame(cbind(rateend, cumratestart, cumrateend))
	names(ret) <- c( ratename, cumrateentername , cumrateendname )
#  print(cbind(Y, startdate, startage, matchdata, ret)[1:10,])
	
	return(ret)
	
	
	
	
}
