## S3 method for class 'flexrsurv'
predict.flexrsurv <- function(object, newdata = NULL,
                              type = c("lp", "link", "risk", "hazard", "hazardrate", "rate", "loghazard", "log",
                                "lograte", "cumulative.rate", "cumulative.hazard", "cumulative", "cum", "survival", "surv", "netsurv"),
                              se.fit = FALSE,
                              na.action = na.pass, ...){

  type <- match.arg(type)

  type <- switch(type,
                 risk = "risk",
                 rate = "risk",
                 hazard = "risk",
                 hazardrate = "risk",
                 lp = "link",
                 log = "link",
                 loghazard = "link",
                 lograte = "link",
                 link = "link",
                 terms = "terms",
                 cumulative.rate="cum",
                 cumulative.hazard="cum",
                 cumulative="cum",
                 cum="cum",
                 surv="surv",
                 survival="surv",
                 netsurv="surv")  

  if (type == "link" | type == "risk" ){
#    if (inherits(object, "glm") & type != "cumulative.rate"){
#      type <- switch(type,
#                     risk = "response",
#                     rate = "response",
#                     hazard = "response",
#                     hazardrate = "response",
#                     lp = "link",
#                     log = "link",
#                     loghazard = "link",
#                     lograte = "link",
#                     link = "link")
#      pred <- predict.glm(object, newdata, type, se.fit, dispersion = NULL, terms, na.action, ...)
#    }
#    else {
      pred <- predictHazard.flexrsurv(object, newdata, type, se.fit,  na.action, ...)
#    }
  } else {
# cummulative rate
    pred <- predictCumulativeHazard.flexrsurv(object, newdata, type, se.fit, na.action, ...)
  }
 pred
}


