/*  Routines for computing the final wheighted cummulative exposure  
 *
 *     Copyright (C) 2015 Michel Grzebyk.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * These functions are distributed in the hope that they will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
 * GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 */

#include <R.h>
#include <Rinternals.h>

#include "flexrsurv.h"


#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("splines", String)
#else
#define _(String) (String)
#endif


SEXP eval_wce_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, 
						   SEXP xvals, SEXP beta, SEXP w, SEXP fromT, SEXP outerok)
{
/* evaluate Weighted Cummulative Exposure with weeights defined by non-zero B-spline basis functions at xvals, 
*/
/* 
 knots : vector of ordered unreplicated INTERIOR knots 
Matrices : a vectorized array of dim order X nbases X number_of_intervales(knots) 
  where nbases is the number of bases of the non integrated, non derived splines 
order : order of the splines (see package orthogonalsplinbasis
intercept : wehtehr first basis is included
xvals : vector values at which bases are computed
beta : vector of the linear combination
*  w,fromT : vectors defining the exposure profile: same menght, fromt increasing, w are the increments
*  it is assumed that min(xvals) > max(fromT) 
 */
	R_len_t i, j, nw, nt, nx /*, oo*/;
	double *rxvals, *rdxvals, *rcl, *rwce, *rw, *rfromT;
	SEXP dxvals, cl, wce;
/*	double  outer_val; */
	
	
	PROTECT(knots = coerceVector(knots, REALSXP));
	PROTECT(order = coerceVector(order, INTSXP));
	PROTECT(intercept = coerceVector(intercept, INTSXP));
	PROTECT(Matrices = coerceVector(Matrices, REALSXP));
	PROTECT(xvals = coerceVector(xvals, REALSXP));
	PROTECT(w	= coerceVector(w, REALSXP));
	PROTECT(fromT = coerceVector(fromT, REALSXP));
	PROTECT(outerok = coerceVector(outerok, LGLSXP));
	
	rxvals = REAL(xvals); 
	nx = length(xvals);

	rw = REAL(w); 
	nw = length(w);

	rfromT = REAL(fromT); 
	nt = length(fromT);
	
	if(nw != nt) {
		error("length of 'W' and 'fromT' differ");    
	} 	
	
	PROTECT(wce = allocVector(REALSXP, nx));
	rwce = REAL(wce);
	
	PROTECT(dxvals = allocVector(REALSXP, nx));
	rdxvals = REAL(dxvals);
	
	PROTECT(cl = allocVector(REALSXP, nx));

/* check if necessary to set the value of wce at outer xvals	
	oo = asLogical(outerok);
	
	if(oo == NA_LOGICAL) {
		error("'outer.ok' must be TRUE or FALSE");    
	} else  if (oo) {
		outer_val = 0.0;
	} else {
		outer_val = R_NaN;
	}
*/

	for(i = 0; i < nx; i++) {
		rwce[i] = 0.0;
	}
	for(j = 0; j < nt; j++) {
		Rprintf("j %d w %4.6f fromT %4.6f \n", j, rw[j], rfromT[j]);
		for(i = 0; i < nx; i++) {
			rdxvals[i] = rxvals[i]-rfromT[j];
		Rprintf("i %d t %4.6f dt %4.6f \n", i, rxvals[i], REAL(dxvals)[i]);
		}
		cl = eval_lc_spline_basis(knots, order, Matrices, intercept, 
											 dxvals, beta, outerok);
	 	rcl = REAL(cl);

		for(i = 0; i < nx; i++) {
			rwce[i] = rwce[i] + rcl[i] * rw[j];
		Rprintf("i %d t %4.6f cl %4.6f  wce %4.6f \n", i, rxvals[i], rcl[i], rwce[i]);
		}
	}
/*	
	for(i = 0; i < nx; i++) {
		Rprintf("i %d t %4.6f wce %4.6f \n", i, rxvals[i], rwce[i]);
	}
*/
	unprotect(11);
	return(wce);
}

SEXP eval_wce_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
								SEXP coefs, SEXP degrees, SEXP intercept, 
								SEXP xvals, SEXP beta, SEXP w, SEXP fromT, SEXP outerok)
{
/* evaluate Weighted Cummulative Exposure with weeights defined by truncated power spline basis functions at xvals, 
*knots : vector of ordered unreplicated INTERIOR knots 
*replicates : vector of the number of replicate of knots 
*order : order of the splines
*min, max : working range, outside of [min, max] the value of the bases is 0 if outerok==TRUE, NA if outerok==FALSE
coefs : vector of coef by which each basis is multiplied : b_i(t) = coef[i] * monomial_i(x)
degrees : vector of the degrees of each monimial : monoial_i(x) = (x  ...)^degrees[i]
*intercept : wehtehr first basis is included
*  xvals : vectors of the time at which to compute WCE
*  beta : vector of the coefficient defining the weights
*  w,fromT : vectors defining the exposure profile: same menght, fromt increasing, w are the increments
*  it is assumed that min(xvals) > max(fromT) 
*/






	R_len_t i, j, nw, nt, nx /*, oo*/;
	double *rxvals, *rdxvals, *rcl, *rwce, *rw, *rfromT;
	SEXP dxvals, wce;
/*	double  outer_val; */

	PROTECT(knots = coerceVector(knots, REALSXP));
	PROTECT(replicates = coerceVector(replicates, REALSXP));
	PROTECT(min = coerceVector(min, REALSXP));
	PROTECT(max = coerceVector(max, REALSXP));
	PROTECT(order = coerceVector(order, INTSXP));
	PROTECT(coefs = coerceVector(coefs, REALSXP));
	PROTECT(degrees = coerceVector(degrees, REALSXP));
	PROTECT(intercept = coerceVector(intercept, INTSXP));
	PROTECT(xvals = coerceVector(xvals, REALSXP));
	PROTECT(beta = coerceVector(beta, REALSXP));
	PROTECT(w	= coerceVector(w, REALSXP));
	PROTECT(fromT = coerceVector(fromT, REALSXP));
	PROTECT(outerok = coerceVector(outerok, LGLSXP));
	
	rxvals = REAL(xvals); 
	nx = length(xvals);
	
	rw = REAL(w); 
	nw = length(w);

	rfromT = REAL(fromT); 
	nt = length(fromT);
	
	if(nw != nt) {
		error("length of 'W' and 'fromT' differ");    
	} 	
	
	PROTECT(wce = allocVector(REALSXP, nx));
	rwce = REAL(wce);
	
	PROTECT(dxvals = allocVector(REALSXP, nx));
	rdxvals = REAL(dxvals);
	

/* check if necessary to set the value of wce at outer xvals	
	oo = asLogical(outerok);
	
	if(oo == NA_LOGICAL) {
		error("'outer.ok' must be TRUE or FALSE");    
	} else  if (oo) {
		outer_val = 0.0;
	} else {
		outer_val = R_NaN;
	}
*/	
	for(i = 0; i < nx; i++) {
		rwce[i] = 0.0;
	}
	for(j = 0; j < nt; j++) {
		for(i = 0; i < nx; i++) {
			rdxvals[i] = rxvals[i]-rfromT[j];
		}
		rcl = REAL(eval_lc_trunc_power_basis(knots, replicates, min, max, order, 
											coefs, degrees, intercept, 
											dxvals, beta, outerok));
		for(i = 0; i < nx; i++) {
			rwce[i] = rwce[i] + rcl[i] * rw[j];
		}
	}
	
		
	UNPROTECT(15);
	return(wce);
}




SEXP eval_wce_trunc_power_increasing_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
											SEXP coefs, SEXP degrees, SEXP intercept, 
											SEXP xvals, SEXP beta, SEXP w, SEXP fromT, SEXP outerok)
{
/* evaluate Weighted Cummulative Exposure with weeights defined by truncated power spline basis functions at xvals, 
* knots : vector of ordered unreplicated INTERIOR knots 
* replicates : vector of the number of replicate of knots 
* order : order of the splines
* min, max : working range, outside of [min, max] the value of the bases is 0 if outerok==TRUE, NA if outerok==FALSE
coefs : vector of coef by which each basis is multiplied : b_i(t) = coef[i] * monomial_i(x)
        coefs include scaling due to integration/derivation
degrees : vector of the degrees of each monimial : monoial_i(x) = (x  ...)^degrees[i]
* intercept : wehtehr first basis is included
*  xvals : vectors of the time at which to compute WCE
*  beta : vector of the coefficient defining the weights
*  w,fromT : vectors defining the exposure profile: same menght, fromt increasing, w are the increments
*  it is assumed that min(xvals) > max(fromT) 
*/
/* all the bases are increasing functions.
for truncated bases, 
	if knots[i] => 0, bi = (x - knotq[i])^i if x>= knots[i], bi = 0 if x<knots[i]
	if knots[i] => 0, bi = -(knotq[i]-x)^i if x<= knots[i], bi = 0 if x>knots[i]
*/

	R_len_t i, j, nw, nt, nx /*, oo*/;
	double *rxvals, *rdxvals, *rcl, *rwce, *rw, *rfromT;
	SEXP dxvals, wce;
/*	double  outer_val; */

	PROTECT(knots = coerceVector(knots, REALSXP));
	PROTECT(replicates = coerceVector(replicates, REALSXP));
	PROTECT(min = coerceVector(min, REALSXP));
	PROTECT(max = coerceVector(max, REALSXP));
	PROTECT(order = coerceVector(order, INTSXP));
	PROTECT(intercept = coerceVector(intercept, INTSXP));
	PROTECT(xvals = coerceVector(xvals, REALSXP));
	PROTECT(beta = coerceVector(beta, REALSXP));
	PROTECT(w	= coerceVector(w, REALSXP));
	PROTECT(fromT = coerceVector(fromT, REALSXP));
	PROTECT(outerok = coerceVector(outerok, LGLSXP));
	
	rxvals = REAL(xvals); 
	nx = length(xvals);
	
	rw = REAL(w); 
	nw = length(w);

	rfromT = REAL(fromT); 
	nt = length(fromT);
	
	if(nw != nt) {
		error("length of 'W' and 'fromT' differ");    
	} 	
	
	PROTECT(wce = allocVector(REALSXP, nx));
	rwce = REAL(wce);
	
	PROTECT(dxvals = allocVector(REALSXP, nx));
	rdxvals = REAL(dxvals);
	

/* check if necessary to set the value of wce at outer xvals	
	oo = asLogical(outerok);
	
	if(oo == NA_LOGICAL) {
		error("'outer.ok' must be TRUE or FALSE");    
	} else  if (oo) {
		outer_val = 0.0;
	} else {
		outer_val = R_NaN;
	}
*/	
	for(i = 0; i < nx; i++) {
		rwce[i] = 0.0;
	}
	for(j = 0; j < nt; j++) {
		for(i = 0; i < nx; i++) {
			rdxvals[i] = rxvals[i]-rfromT[j];
		}
		rcl = REAL(eval_lc_trunc_power_increasing_basis(knots, replicates, min, max, order,  
											coefs, degrees, intercept, 
											dxvals, beta, outerok));
		for(i = 0; i < nx; i++) {
			rwce[i] = rwce[i] + rcl[i] * rw[j];
		}
	}
	
		
	UNPROTECT(15);
	return(wce);
}


