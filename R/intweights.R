intweights_CAV_SIM<-function(nl,step){
	# wehights for cavalieri simpson rule (nl even)
	# weights are (1 4 2 4 2 4 ... 4 2 4 2 4 1 )
	w<-rep(c(2L/3L, 4L/3L)*step, nl%/% 2)
	w[1]<-w[nl+1]<-1L/3L*step
	w
}

intweights_SIM_3_8<-function(nl,step){
	# wehights for 3/8 simpson rule (nl = 3 * i)
	# weights are (1 3 3 2 3 3 2 3 3  ... 3 3 2 3 3 2 3 3 1 )
	w<-rep(c( 3L/4L, 9L/8L, 9L/8L)*step, nl%/%3)
	w[1]<-w[nl+1]<-3L/8L*step
	w 
}

intweights_SIM_3_8b<-function(nl,step){
	# wehights for 3/8 simpson rule (nl = 3 * i)
	# weights are (1 3 3 2 3 3 2 3 3  ... 3 3 2 3 3 2 3 3 1 )
	w<-rep(9L/8L*step, nl+1)
	w[1+3*(1:(nl/3-1))]<-3L/4L*step
	w[1]<-w[nl+1]<-3L/8L*step
	w 
}


intweights_BOOLE<-function(nl,step){
	# wehights for BOOLE rule (nl = 4 * i)
	# weights are (7 32 12 32 14 32 12 32 14 ... 14 32 12 32 14 32 12 32 7 )
	w<-rep(c(28L/45L, 64L/45L, 8L/15L, 64L/45L)*step,  nl %/% 4)
	w[1]<-w[nl+1]<-14L/45L*step
	w
}
intweights_BOOLEb<-function(nl,step){
	# wehights for BOOLE rule (nl = 4 * i)
	# weights are (7 32 12 32 14 32 12 32 14 ... 14 32 12 32 14 32 12 32 7 )
	w<-rep(64L/45L*step, nl+1)
	w[4*1:(nl/4)-1]<-8L/15L*step
	w[4*1:(nl/4)+1]<-28L/45L*step
	w[1]<-w[nl+1]<-14L/45L*step
	w
}
