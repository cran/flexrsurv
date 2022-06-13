/*  Routines for manipulating B-splines.  These are intended for use with
 *  S or S-PLUS or R.
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
 * The routines are loosely based on the pseudo-code in Schumacher (Wiley,
 * 1981) and the CMLIB library DBSPLINES.
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

#ifndef EVAL_SPLINEPARAM
#include "SplineParam.h"
#endif


SEXP eval_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
		SEXP coefs, SEXP degrees, SEXP intercept, SEXP xvals, SEXP outerok)
{
	/* evaluate the non-zero truncated power spline basis functions at xvals, */
	/*
knots : vector of ordered unreplicated INTERIOR knots 
replicates : vector of the number of replicate of knots 
order : order of the splines
min, max : working range, outside of [min, max] the value of the bases is 0 if outerok==TRUE, NA if outerok==FALSE
coefs : vector of coef by which each basis is multiplied : b_i(t) = coef[i] * monomial_i(x)
degrees : vector of the degrees of each monimial : monoial_i(x) = (x  ...)^degrees[i]
xvals : vector values at which bases are computed
	 */

	R_len_t i, j, k, ibase, icoef, nknots, theorder, nbases, nx, oo;
	R_len_t theinterval, firstbasis, mfl;
	double *rknots, *rcoefs, *rxvals, *rbases, *rreplicates;
	double rmin, rmax;
	SEXP bases;
	double  outer_val;

	PROTECT(knots = coerceVector(knots, REALSXP));
	PROTECT(replicates = coerceVector(replicates, REALSXP));
	PROTECT(min = coerceVector(min, REALSXP));
	PROTECT(max = coerceVector(max, REALSXP));
	PROTECT(coefs = coerceVector(coefs, REALSXP));
	PROTECT(degrees = coerceVector(degrees, REALSXP));
	PROTECT(order = coerceVector(order, INTSXP));
	PROTECT(intercept = coerceVector(intercept, INTSXP));
	PROTECT(xvals = coerceVector(xvals, REALSXP));
	PROTECT(outerok = coerceVector(outerok, LGLSXP));

	rknots = REAL(knots); 
	nknots = length(knots);
	rreplicates = REAL(replicates);
	theorder = INTEGER(order)[0];
	rmin = REAL(min)[0];
	rmax = REAL(max)[0];

	rcoefs = REAL(coefs); 

	/* number of bases */
	nbases = theorder ;	
	for( i=0; i<nknots;  i++){
		nbases = nbases + rreplicates[i];
	}
	rxvals = REAL(xvals); 
	nx = length(xvals);

	firstbasis = (INTEGER(intercept)[0]==0);	

	PROTECT(bases = allocMatrix(REALSXP, nx, nbases-firstbasis));
	rbases = REAL(bases);

	oo = asLogical(outerok);

	if(oo == NA_LOGICAL) {
		error("'outer.ok' must be TRUE or FALSE");    
	} else  if (oo) {
		outer_val = 0.0;
	} else {
		outer_val = R_NaN;
	}

	for(i = 0; i < nx; i++) {
		EVALUATE_one_trunc_power_basis (rxvals[i], rbases, i + nx *)   ;

	}
	unprotect(11);
	return(bases);
}




SEXP eval_trunc_power_increasing_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
		SEXP coefs, SEXP degrees, SEXP intercept, SEXP xvals, SEXP outerok)
{
	/* evaluate the non-zero truncated power spline basis functions at xvals, */
	/* all the bases are increasing functions.
for truncated bases, 
	if knots[i] => 0, bi = (x - knotq[i])^i if x>= knots[i], bi = 0 if x<knots[i]
	if knots[i] => 0, bi = -(knotq[i]-x)^i if x<= knots[i], bi = 0 if x>knots[i]

	 */
	/*
knots : vector of ordered unreplicated INTERIOR knots 
replicates : vector of the number of replicate of knots 
order : order of the splines
min, max : working range, outside of [min, max] the value of the bases is 0 if outerok==TRUE, NA if outerok==FALSE
coefs : vector of coef by which each basis is multiplied : b_i(t) = coef[i] * monomial_i(x)
        coefs include scaling due to integration/derivation
degrees : vector of the degrees of each monimial : monoial_i(x) = (x  ...)^degrees[i]
xvals : vector values at which bases are computed
	 */

	R_len_t i, j, k, ibase, icoef, nknots, theorder, nbases, nx, oo;
	R_len_t theinterval, firstbasis, mfl;
	double *rknots, *rcoefs, *rxvals, *rbases, *rreplicates;
	double rmin, rmax;
	SEXP bases;
	double  outer_val;

	PROTECT(knots = coerceVector(knots, REALSXP));
	PROTECT(replicates = coerceVector(replicates, REALSXP));
	PROTECT(min = coerceVector(min, REALSXP));
	PROTECT(max = coerceVector(max, REALSXP));
	PROTECT(order = coerceVector(order, INTSXP));
	PROTECT(coefs = coerceVector(coefs, REALSXP));
	PROTECT(degrees = coerceVector(degrees, REALSXP));
	PROTECT(intercept = coerceVector(intercept, INTSXP));
	PROTECT(xvals = coerceVector(xvals, REALSXP));
	PROTECT(outerok = coerceVector(outerok, LGLSXP));

	rknots = REAL(knots); 
	nknots = length(knots);
	rreplicates = REAL(replicates);
	theorder = INTEGER(order)[0];
	rmin = REAL(min)[0];
	rmax = REAL(max)[0];

	rcoefs = REAL(coefs); 

	/* number of bases */
	nbases = theorder;	
	for( i=0; i<nknots;  i++){
		nbases = nbases + rreplicates[i];
	}
	rxvals = REAL(xvals); 
	nx = length(xvals);

	firstbasis = (INTEGER(intercept)[0]==0);	

	PROTECT(bases = allocMatrix(REALSXP, nx, nbases-firstbasis));
	rbases = REAL(bases);

	oo = asLogical(outerok);

	if(oo == NA_LOGICAL) {
		error("'outer.ok' must be TRUE or FALSE");    
	} else  if (oo) {
		outer_val = 0.0;
	} else {
		outer_val = R_NaN;
	}

	for(i = 0; i < nx; i++) {
		EVALUATE_one_trunc_power_increasing_basis (rxvals[i], rbases, i + nx *)   ;
	}

	unprotect(11);
	return(bases);
}
