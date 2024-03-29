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

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("splines", String)
#else
#define _(String) (String)
#endif




SEXP eval_lc_linex_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP linexinf, SEXP linexsup, SEXP intercept, SEXP xvals, SEXP beta, SEXP outerok)
{
	/* evaluate linear combination of the non-zero B-spline basis functions at xvals, using matrices generated by packag::orthogonalsplinbasis */
	/*
 knots : vector of ordered unreplicated INTERIOR knots 
Matrices : a vectorized array of dim order X nbases X number_of_intervales(knots) 
  where nbases is the number of bases of the non integrated, non derived splines 
order : order of the splines (see package orthogonalsplinbasis
intercept : wehtehr first basis is included
xvals : vector values at which bases are computed
beta : vector of the linear combination
	 */
	R_len_t i, j, k, nknots, theorder, nbases, nbeta, nx;
	R_len_t theinterval, firstbasis, mfl;
	double *rknots, *rMatrices, *rlinexinf, *rlinexsup, *rxvals, *rbeta, *rcl;
	SEXP cl;
	SEXP dims;
	double temp, tempcl, *U, u;


	PROTECT(knots = coerceVector(knots, REALSXP));
	PROTECT(order = coerceVector(order, INTSXP));
	PROTECT(intercept = coerceVector(intercept, INTSXP));
	PROTECT(Matrices = coerceVector(Matrices, REALSXP));
	PROTECT(linexinf = coerceVector(linexinf, REALSXP));
	PROTECT(linexsup = coerceVector(linexsup, REALSXP));
	PROTECT(xvals = coerceVector(xvals, REALSXP));
	PROTECT(beta = coerceVector(beta, REALSXP));
	PROTECT(outerok = coerceVector(outerok, LGLSXP));


	rknots = REAL(knots); 
	nknots = length(knots);
	theorder = INTEGER(order)[0];

	dims = getAttrib(Matrices, R_DimSymbol);
	if( LENGTH(dims) < 3 ){
		error("'Matrices' must be an array with 3 dim");   
	}
	nbases = INTEGER(dims)[1];

	rxvals = REAL(xvals); 
	nx = length(xvals);

	rbeta = REAL(beta); 
	nbeta = length(beta);


	firstbasis = (INTEGER(intercept)[0]==0);
	rMatrices = REAL(Matrices);
	rlinexinf = REAL(linexinf);
	rlinexsup = REAL(linexsup);

	PROTECT(cl = allocVector(REALSXP, nx));
	rcl = REAL(cl);

	if(nbeta < nbases-firstbasis) {
		error("length of 'beta' must be at least %d", nbases-firstbasis);    
	} 	


	U = (double *) R_alloc( theorder, sizeof(double));
	/* outer is linear
	oo = asLogical(outerok);
	if(oo == NA_LOGICAL) {
		error("'outer.ok' must be TRUE or FALSE");    
	} else  if (oo) {
		outer_val = 0.0;
	} else {
		outer_val = R_NaN;
	}
	 */
	U[0]=1.0;
	for(i = 0; i < nx; i++) {
		if (ISNAN(rxvals[i])) {
			rcl[i] = R_NaN;
		}
		else {
			theinterval= 1;
			mfl = 0;
			/* find the interval within the range of all the knots (which include boundaries)
   of rxvals[i], rightmost_close=TRUE, all_inside = FALSE */ 
			theinterval = findInterval(rknots, nknots, rxvals[i], 1, 0 , theinterval, &mfl );
			if (theinterval == 0) {
				/* before first boundary knot, linear */
				u = (rxvals[i]-rknots[0]);
				for ( j = 1; j < theorder ; j++) {
					U[j] = pow(u, (double)j);
				}
				tempcl = 0;
				for (k = firstbasis; k < nbases; k++) {
					temp = 0;
					for (int j = 0; j < theorder ; j++) {
						temp += U[j] * rlinexinf[theorder*k + j];
					}
					tempcl = tempcl + temp * rbeta[k-firstbasis];
				}
				rcl[i] = tempcl;

				/*
			    if( firstbasis == 0){
				 */
				/* first basis rcl[i + nx * 0] */
				/*
				    rcl[i] = (rlinexinf[0] + rlinexinf[2]*(rxvals[i]-rknots[0])) * rbeta[0];
				 */
				/* second basis rcl[i + nx * i] */
				/*
				    rcl[i] = rcl[i] + ( rlinexinf[1] + rlinexinf[3]*(rxvals[i]-rknots[0])) * rbeta[1];
			    } else {
				 */
				/* first basis discard */
				/* second basis in rcl[i + nx * (1-firstbasis)] */
				/*
				    rcl[i] = ( rlinexinf[1] + rlinexinf[3]*(rxvals[i]-rknots[0])) * rbeta[0];
			    }
				 */
			} else if (theinterval == nknots  ) {
				/* after last boundary knot, linear */
				u = (rxvals[i]-rknots[nknots-1]);
				for ( j = 1; j < theorder ; j++) {
					U[j] = pow(u, (double)j);
				}
				tempcl = 0;
				for (k = firstbasis; k < nbases; k++) {
					temp = 0;
					for (int j = 0; j < theorder ; j++) {
						temp += U[j] * rlinexsup[theorder*k + j];
					}
					tempcl = tempcl + temp * rbeta[k-firstbasis];
				}
				rcl[i] = tempcl;

				/* last basis rcl[i + nx * nbases] */
				/*
			    rcl[i] = (rlinexsup[1] + rlinexsup[3]*(rxvals[i]-rknots[nknots-1])) * rbeta[(nbases-firstbasis-1)];
				 */
				/* last basis rcl[i + nx * nbases] */
				/*
			    rcl[i] = rcl[i] + (rlinexsup[0] + rlinexsup[2]*(rxvals[i]-rknots[nknots-1])) * rbeta[(nbases-firstbasis-2)] ;
				 */
			}
			else {
				if( theinterval == nknots -1) {
					/* xx[i] is the rightmost boundary knot */
					theinterval = nknots - theorder;
				}
				u = (rxvals[i] - rknots[theinterval-1])/(rknots[theinterval]-rknots[theinterval-1]);
				for ( j = 1; j < theorder ; j++) {
					U[j] = pow(u, (double)j);
				}
				/* the usefull matrix is the (theinterval - theorder +1)th matrix of Matrices */
				theinterval = theinterval - theorder;
				tempcl = 0;
				for (k = firstbasis; k < nbases; k++) {
					temp = 0;
					for (int j = 0; j < theorder ; j++) {
						temp += U[j] * rMatrices[theorder*nbases*theinterval+ theorder*k + j];
					}
					tempcl = tempcl + temp * rbeta[k-firstbasis];
				}
				rcl[i] = tempcl;
			}
		}

	}

	UNPROTECT(10);
	return(cl);
}
