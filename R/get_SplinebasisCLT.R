# functions get_Splinebasis to get spline basis of the relational model fot the life table correction model

get_SplinebasisCLT <- function(objfit, 
		data=parent.frame()){
	# get spline parameters of the relational model fot the life table correction model
# input
#      objfit : a fit object of class flexrsurvclt 
# output  : "SplineBasis" objects
	
	mfclt <- match.call(flexrsurvclt, call= objfit$call, expand.dots = FALSE)
	mclt <- match(c("formula.table", "logit_start", "logit_end", "knots.table", "degree.table", "Spline.table", "Spline.CLT", "model.correction"),
			names(mfclt), 0L)
	
	if(sum(mclt[-c(2,3)]) > 0) {
		is_correction_model <- TRUE
		Intercept_B <- TRUE
		
		# analysis with correction model both model_correction =="period" & "cohorte" 	  
		

		if(sum(mclt[4:7])>0){
			if (mclt[7]!=0){
				Spline_B <- eval(mfclt$Spline.CLT)
			} else {
				if(mclt[4]!=0 & mclt[5]!=0) {
					Spline_B <- R2bBSplineBasis(knots=eval(mfclt$knots.table), degree=eval(mfclt$degree.table))
				}
				else {
					stop ("With correction of life table models, both 'knots.table' and 'degree.table' are  required.")
				}
			}
			# df for the brass model (fisrt basis has coef equal to one
			nbrass <- getNBases(Spline_B) - 1 - (1 - Intercept_B)
			if(mclt[8]!=0){
				model_correction <- eval(mfclt$model_correction) 			
			} else {
				model_correction <- "cohort" 		
			}
			
			
		}
		else {
			Spline_B <- NULL
			nbrass <- 0
		}
		
		
		
	}
	
	return(list(Spline = Spline_B, Intercept_B = Intercept_B, nbrass = nbrass, model_correction=model_correction))
}

