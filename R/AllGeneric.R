#### Define those generics that we need, if they don't exist;
#### not all will be exported

setGeneric("evaluate",function(object, x,...)standardGeneric("evaluate"))
setGeneric("fevaluate",function(object, x,...)standardGeneric("fevaluate"))

setGeneric("evaluatelc",function(object, x, beta,...)standardGeneric("evaluatelc"))

setGeneric("predictSpline",function(object, x, beta,...)standardGeneric("predictSpline"))
setGeneric("slowpredictSpline",function(object, x, beta,...)standardGeneric("slowpredictSpline"))

setGeneric("integrate",function(object, x,...)standardGeneric("integrate"))
setGeneric("initcoef",function(object, ncol,...)standardGeneric("initcoef"))
setGeneric("initcoefC",function(object, ncol, ...)standardGeneric("initcoefC"))

setGeneric("predictwce",function(object, t, Increment, fromT, ...)standardGeneric("predictwce"))
setGeneric("predictcumwce",function(object, t, Increment, fromT, ...)standardGeneric("predictcumwce"))
setGeneric("gradientwce",function(object, t, Increment, fromT, ...)standardGeneric("gradientwce"))

setGeneric("predictCLT",function(object,...)standardGeneric("predictCLT"))
