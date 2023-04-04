#  modified from File src/library/stats/R/constrOptim.R
#  modification allow for grad to be NULL
# modification to alloaw for non linear constraint
# implementation of the Adaptive Logarithmic  Barrier Methods 
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

# optimize f(par) subject to ui_k(par) >= ci_k
#     b(par) = Sum( ui_k(par.old) ln(ui_k(par) - grad_ui_k(par.old)%*%par)
# barrier algorithm : optimize f(par) - mu b(par)  

nonLinconstrOptim <-
		function(par, f, grad, f_ui, ci, grad_ui, mu = 0.0001, control = list(),
				method = if(is.null(grad)) "Nelder-Mead" else "BFGS",
				outer.iterations = 100, outer.eps = 0.00001, 
				link.barrier = c("log", "power"), alpha.barrier = 1/2, beta.barrier = 1, ..., hessian = FALSE)
{
	
	link.barrier<- match.arg(link.barrier)
	
	if( mu < 0) stop("mu must be positive")
	if (!is.null(control$fnscale) && control$fnscale < 0){
		mu <- -mu ##maximizing
		fnscale <- -1
	}
	else {
		fnscale <- 1
	}
	
	if(link.barrier == "log"){
		g_bar <- function(par, par.old, ...) {
			bari <-  f_ui(par, ...) - ci
			bari.old <- f_ui(par.old, ...) - ci			
			if (any(bari<0)) return(NaN)
			bar <- sum(bari.old*log(bari) - grad_ui(par.old, ...) %*% (par - par.old) )
			if (!is.finite(bar)) bar <-  -Inf
			f(par, ...) - mu*bar
		}
		
		grad_g_bar <- function(par, par.old, ...) {
			bari <- drop(f_ui(par, ...) - ci)
			bari.old <- drop(f_ui(par.old, ...) - ci)
			dbar <- colSums(grad_ui(par, ...)*bari.old/bari - grad_ui(par.old, ...))
			grad(par, ...) - mu*dbar
		}
	} else {
		g_bar <- function(par, par.old, ...) {
			bari <-  f_ui(par, ...) - ci
			if (any(bari<0)) return(NaN)
			bari.old <- f_ui(par.old, ...) - ci
			
			bar <- sum(bari.old^(alpha.barrier+beta.barrier)/bari^alpha.barrier + alpha.barrier * bari.old^(beta.barrier-1) * (grad_ui(par.old, ...) %*% (par -par.old))  )
			if (!is.finite(bar)) bar <-  -Inf
			f(par, ...) - mu*bar
		}
		
		grad_g_bar <- function(par, par.old, ...) {
			bari <- drop(f_ui(par, ...) - ci)
			bari.old <- drop(f_ui(par.old, ...) - ci)
			dbar <- alpha.barrier * colSums(grad_ui(par, ...)*bari.old^(alpha.barrier+beta.barrier)/bari^(alpha.barrier+1) - grad_ui(par.old, ...)* bari.old^(beta.barrier-1) )
			grad(par, ...) - mu*dbar
		}
		
	}
	if (any(f_ui(par, ...) - ci <= 0))
		stop("the initial value is not in the feasible region")
	
	f.par <- f(par, ...) * fnscale
	fbar.par <- g_bar(par, par, ...)
	ffbar.par <- fbar.par
	totCounts <- 0
	optobj <- NULL
	
	for(i in seq_len(outer.iterations)) {
		par.old <- par
		f.par.old <- f.par
		fbar.par.old <- fbar.par
		ffbar.par.old <- ffbar.par
		optobj.old <- optobj
		
		# surrogate objective and gradient functions
		sur.fun <- function(par, ...) g_bar(par, par.old, ...)
		sur.gradient <- if (missing(grad) || is.null(grad)){
					NULL
				} 
				else if (method == "SANN") {
					grad
				}
				else function(par, ...) grad_g_bar(par, par.old, ...)
		
		optobj <- optim(par = par.old, fn = sur.fun, gr = sur.gradient, control = control,
				method = method, hessian = FALSE, ...)
		
		fbar.par <- optobj$value	
		par <- optobj$par
		f.par <- f(par, ...) * fnscale
		
		totCounts <- totCounts + optobj$counts
		
		ffbar.par <- sur.fun(par, par.old, ...)

		if (is.nan(optobj$value) ){
			break
		}
		
		
		if (is.finite(fbar.par) && is.finite(fbar.par.old) &&
				abs(fbar.par - fbar.par.old) < (1e-3 + abs(fbar.par)) * outer.eps){
			break
		}
		if (f.par > f.par.old){
			break
		} 
	}
	
	# the behaviour within optim() is strange when the returned objective value optobj$value is NaN 
	if (is.nan(optobj$value) ){
	# optobj$value = NaN when the all the tested candidates in the search direction are out of the feaseable region
	# In this case several constraints maight be not fulfilled
#		if(any(f_ui(param=optobj$par, Spline=SS, ibrass0=1:13) < ci)){
#			# unable to update in last inner iteration, retrieving results from the start of the outer iteration 
#			if(i > 1){
#				# it was not the first outer iteration, retrieving resulsts from start of the outer iteration (i.e. previous outer iteration)
#				message <- gettextf("Unable to update at last inner iteration of outer iteration %d, retrieving results from outer iteration %d", i, i-1)
#				totCounts <- totCounts - optobj$counts
#				optobj <- optobj.old
#				i <- i-1
#			}
#			else {
#				# it was the fist outer iteration, init value not updated
#				optobj$par <- par.old
#				f.par <- f.par.old  
#				fbar.par <- fbar.par.old 
#				fbar.par <- ffbar.par.old 
#				optobj$value <- ffbar.par
#				message <- gettextf("Unable to update at first outer iteration, init value was not updated")
#				i <- 0
#			}
#			optobj$convergence <- 8
#			optobj$message <- message
#		} else {
#			# optobj$par is a valid parameter but the objective is not computed for this parameter, simply updade optobj$value
#			fbar.par <- optobj$value	
#			par <- optobj$par
#			f.par <- f(par, ...) * fnscale
#			ffbar.par <- sur.fun(par, par.old, ...)
#			optobj$convergence <- 9
#			optobj$message <- gettextf("Unable to update at outer iteration %d", i)
#		}
	}
	
	optobj$outer.iterations <- i
	optobj$counts <- totCounts
	optobj$value <- f.par * fnscale
	optobj$barrier.value <- ffbar.par - optobj$value	
	
	if (i == outer.iterations) {
		optobj$convergence <- 7
		optobj$message <- gettext("Barrier algorithm ran out of iterations and did not converge")
	}
	if (f.par > f.par.old) {
		optobj$convergence <- 11
		if(fnscale > 0){
			optobj$message <- gettextf("Objective function increased at outer iteration %d", i)
		} else {
			optobj$message <- gettextf("Objective function decreased at outer iteration %d", i)
		}
	}
	
	if(hessian){
		optobj$hessian <- optimHess(par=optobj$par, fn = f, gr = grad, control = control,
				method = method, ...)
#		surhessian <- optimHess(par=optobj$par, fn = sur.fun, gr = sur.gradient, control = control,
#				method = method, ...)
		
	}
	
	
	optobj
}


