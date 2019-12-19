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

#ifndef EVAL_SPLINEPARAM
#include "SplineParam.h"
#endif

SEXP predict_wce_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, 
			   SEXP w, SEXP fromT, SEXP FirstId, SEXP LastId, 
			   SEXP xvals, SEXP xId,  
			   SEXP outerok)
{
/* evaluate Weighted Cummulative Exposure with weights defined by non-zero B-spline basis functions at xvals, 
*/
/* 
 knots : vector of ordered unreplicated INTERIOR knots 
Matrices : a vectorized array of dim order X nbases X number_of_intervales(knots) 
  where nbases is the number of bases of the non integrated, non derived splines 
!!!! The Matrix is scaled by the parameters of the spline function  !!!
order : order of the splines (see package orthogonalsplinbasis
intercept : wehtehr first basis is included
*  w,fromT, FirstId, LastId: vectors defining the exposure profiles if several individual: 
*                     same lenght, fromT,  increasing by Id, w are the increments
*                     fromT[FirstId] is the first time of each individual
*                     fromT[LastId] is the last enter time of each individual
*  xvals : vector values at which wce is evaluated
*  xId : index such that each xvals is in ] fromT[xId], toT[xId]]
*        the Id of each xvals is Id[xId]
*  it is assumed that for each Id, the profile is such that the last interval is right-open [fromT[last] + infty[
*             (ie that xvals <= toT[last])  
 */
	R_len_t nknots, theorder, nbases, nintervals, firstbasis;
	R_len_t i, j, k, nw, nt, nx, istart, ilast, ixId, nxId, nFId, oo, ii;
	int *rxId, *rFirstId, *rLastId;
	double *rxvals, dxval, dwce, *rwce, *rw, *rfromT;
	double *rknots, *rMatrices, *rAddMatrices;
	SEXP dims;
	SEXP wce;
	double  outer_val; 
	R_len_t theinterval, mfl;
	double temppredict, U, u;
	
	
	PROTECT(knots = coerceVector(knots, REALSXP));
/*	PROTECT(order = coerceVector(order, INTSXP)); */
	PROTECT(Matrices = coerceVector(Matrices, REALSXP));
/*	PROTECT(intercept = coerceVector(intercept, INTSXP)); */
	PROTECT(xvals = coerceVector(xvals, REALSXP));
	PROTECT(xId = coerceVector(xId, INTSXP));
	PROTECT(w	= coerceVector(w, REALSXP));
	PROTECT(fromT = coerceVector(fromT, REALSXP));
	PROTECT(FirstId = coerceVector(FirstId, INTSXP));
	PROTECT(LastId = coerceVector(LastId, INTSXP));
/*	PROTECT(outerok = coerceVector(outerok, LGLSXP)); */
	
	rknots = REAL(knots); 
	nknots = length(knots);
	theorder = asInteger(order);
	
	dims = getAttrib(Matrices, R_DimSymbol);
	if( LENGTH(dims) < 3 ){
		error("'Matrices' must be an array with 3 dim");   
	}
	nbases = INTEGER(dims)[1];
	nintervals = INTEGER(dims)[2];
	
/* first basis to start with : 
 *     if intercept =TRUE, all bases                -> firstbasis = 0, 
 *     if intercept =FALSE, remove base number one  -> firstbasis = 1   */
	ii = asLogical(intercept);
	if(ii == NA_LOGICAL) {
		error("'intercept' must be TRUE or FALSE");    
	} else   if (ii) {
		firstbasis = 0;
	} else {
		firstbasis = 1;
	}
	rMatrices = REAL(Matrices);
        
	rxvals = REAL(xvals); 
	nx = length(xvals);

	rxId = INTEGER(xId);
	nxId = length(xId);
	
	rw = REAL(w); 
	nw = length(w);

	rfromT = REAL(fromT); 
	nt = length(fromT);
	
	rFirstId = INTEGER(FirstId);
	nFId = length(FirstId);

	rLastId = INTEGER(LastId);
	
	if(nw != nt) {
		error("length of 'W' and 'fromT' differ");    
	} 	
	
	if(nFId != nt) {
		error("length of 'FirstId' and 'fromT' differ");    
	} 	
	
	if(nxId != nx) {
		error("length of 'xId' and 'xvals' differ");    
	} 	

	PROTECT(wce = allocVector(REALSXP, nx));
	rwce = REAL(wce);

/* check if necessary to set the value of wce at outer xvals */	
	oo = asLogical(outerok);
	
	if(oo == NA_LOGICAL) {
		error("'outer.ok' must be TRUE or FALSE");    
	} else  if (oo) {
		outer_val = 0.0;
	} else {
		outer_val = R_NaN;
	}

/* first add cols of each matrix */
    rAddMatrices = (double *) R_alloc( theorder * nintervals, sizeof(double));
    for (i = 0; i < nintervals; i++) {
        for (int j = 0; j < theorder ; j++) {
            rAddMatrices[theorder*i + j] = 0.0;
            for (k = firstbasis; k < nbases; k++) {
                rAddMatrices[theorder*i + j] += rMatrices[theorder*nbases*i+ theorder*k + j];
            }
        }
    }

	for(i = 0; i < nx; i++) {
		rwce[i] = 0.0;

		/* starting index */
		ixId = rxId[i]-1;
		istart = rFirstId[ixId] -1;
		ilast  = rLastId[ixId];
		for( j = istart ; j<ilast && rfromT[j] < rxvals[i] ; j++){
			dxval = rxvals[i]-rfromT[j];
			PREDIC_one_espline_basis(dxval, dwce) ;
			rwce[i] += rw[j] * dwce;
		}
	}
	UNPROTECT(9);
	return(wce);
}

SEXP predict_wce_espline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, 
			   SEXP w, SEXP fromT, SEXP FirstId, SEXP LastId,  
			   SEXP xvals, SEXP xId)
{
/* evaluate Weighted Cummulative Exposure with weights defined by non-zero B-spline basis functions at xvals, 
*/
/* 
 knots : vector of ordered unreplicated INTERIOR knots 
Matrices : a vectorized array of dim order X nbases X number_of_intervales(knots) 
  where nbases is the number of bases of the non integrated, non derived splines 
!!!! The Matrix is scaled by the parameters of the spline function  !!!
order : order of the splines (see package orthogonalsplinbasis
intercept : wehtehr first basis is included
*  w,fromT, FirstId, LastId, Id : vectors defining the exposure profiles if several individual: 
*                     same lenght, fromt, toT  increasing by Id, w are the increments
*                     fromT[FirstId] is the first time of each individual
*  xvals : vector values at which wce is evaluated
*  xId : index such that each xvals is in ] fromT[xId], toT[xId]]
*        the Id of each xvals is Id[xId]
*  it is assumed that each xvals <= toT[xId] 
 */
	R_len_t nknots, theorder, nbases, nintervals, firstbasis;
	R_len_t i, j, k, nw, nt, nx, istart, ilast, ixId, nxId, nFId, ii;
	int *rxId, *rFirstId, *rLastId;
	double *rxvals, dxval, dwce, *rwce, *rw, *rfromT;
	double *rknots, *rMatrices, *rAddMatrices;
	SEXP dims;
	SEXP wce;
	R_len_t theinterval, mfl;
	double temppredict, U, u;
	
	double outer_val;
	outer_val = 0.0;

	PROTECT(knots = coerceVector(knots, REALSXP));
/*	PROTECT(order = coerceVector(order, INTSXP)); */
	PROTECT(Matrices = coerceVector(Matrices, REALSXP));
/*	PROTECT(intercept = coerceVector(intercept, INTSXP)); */
	PROTECT(xvals = coerceVector(xvals, REALSXP));
	PROTECT(xId = coerceVector(xId, INTSXP));
	PROTECT(w	= coerceVector(w, REALSXP));
	PROTECT(fromT = coerceVector(fromT, REALSXP));
	PROTECT(FirstId = coerceVector(FirstId, INTSXP));
	PROTECT(LastId = coerceVector(LastId, INTSXP));
	
	rknots = REAL(knots); 
	nknots = length(knots);
	theorder = INTEGER(order)[0];
	
	dims = getAttrib(Matrices, R_DimSymbol);
	if( LENGTH(dims) < 3 ){
		error("'Matrices' must be an array with 3 dim");   
	}
	nbases = INTEGER(dims)[1];
	nintervals = INTEGER(dims)[2];
	
/* first basis to start with : 
 *     if intercept =TRUE, all bases                -> firstbasis = 0, 
 *     if intercept =FALSE, remove base number one  -> firstbasis = 1   */
	ii = asLogical(intercept);
	if(ii == NA_LOGICAL) {
		error("'intercept' must be TRUE or FALSE");    
	} else   if (ii) {
		firstbasis = 0;
	} else {
		firstbasis = 1;
	}
	rMatrices = REAL(Matrices);

	rxvals = REAL(xvals); 
	nx = length(xvals);

	rxId = INTEGER(xId);
	nxId = length(xId);
	
	rw = REAL(w); 
	nw = length(w);

	rfromT = REAL(fromT); 
	nt = length(fromT);
	
	rFirstId = INTEGER(FirstId);
	nFId = length(FirstId);
	
	rLastId = INTEGER(LastId);

	if(nw != nt) {
		error("length of 'W' and 'fromT' differ");    
	} 	
	
	if(nFId != nt) {
		error("length of 'FirstId' and 'fromT' differ");    
	} 	
	
	if(nxId != nx) {
		error("length of 'xId' and 'xvals' differ");    
	} 	

	PROTECT(wce = allocVector(REALSXP, nx));
	rwce = REAL(wce);
	
	
/* first add cols of each matrix */
    rAddMatrices = (double *) R_alloc( theorder * nintervals, sizeof(double));
    for (i = 0; i < nintervals; i++) {
        for (int j = 0; j < theorder ; j++) {
            rAddMatrices[theorder*i + j] = 0.0;
            for (k = firstbasis; k < nbases; k++) {
                rAddMatrices[theorder*i + j] += rMatrices[theorder*nbases*i+ theorder*k + j];
            }
        }
    }



	for(i = 0; i < nx; i++) {
		rwce[i] = 0.0;

		/* starting index */
		ixId = rxId[i];
		istart = rFirstId[ixId]-1;
		ilast  = rLastId[ixId];

		for( j = istart ; j<ilast && rfromT[j] < rxvals[i] ; j++){
			dxval = rxvals[i]-rfromT[j];
			PREDIC_one_espline_basis(dxval, dwce);
			rwce[i] += rw[j] * dwce;
		}
	}
	unprotect(9);
	return(wce);
}

SEXP predict_wce_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
				SEXP coefs, SEXP degrees, SEXP intercept, 
				SEXP w, SEXP fromT, SEXP FirstId, SEXP LastId, 
				SEXP xvals, SEXP xId,  
				SEXP outerok)
{
/* evaluate Weighted Cummulative Exposure with weeights defined by truncated power spline basis functions at xvals, 
*knots : vector of ordered unreplicated INTERIOR knots 
*replicates : vector of the number of replicate of knots 
*order : order of the splines
*min, max : working range, outside of [min, max] the value of the bases is 0 if outerok==TRUE, NA if outerok==FALSE
coefs : vector of coef by which each basis is multiplied : b_i(t) = coef[i] * monomial_i(x)
!!!! The coefs is scaled by the parameters of the spline function  !!! Spline(x) = sum(coef[i] * monomial_i(x))
degrees : vector of the degrees of each monimial : monoial_i(x) = (x  ...)^degrees[i]
*intercept : wehtehr first basis is included
*  w,fromT, FirstId, LastId, Id : vectors defining the exposure profiles if several individual: 
*                     same lenght, fromT,  increasing by Id, w are the increments
*                     fromT[FirstId] is the first time of each individual
*                     fromT[LastId] is the last enter time of each individual
*  xvals : vector values at which wce is evaluated
*  xId : index such that each xvals is in ] fromT[xId], toT[xId]]
*        the Id of each xvals is Id[xId]
*  it is assumed that for each Id, the profile is such that the last interval is right-open [fromT[last] + infty[
*             (ie that x_splinevals <= toT[last])  
*/



	R_len_t i, j, k, ibase, icoef, nknots, theorder, nbases;
	R_len_t nw, nt, nx, istart, ilast, ixId, nxId, nFId, oo, ii;
	int *rxId, *rFirstId, *rLastId;
	double *rxvals, dxval, dwce, *rwce, *rw, *rfromT;
	SEXP wce;
	double  outer_val; 

	R_len_t theinterval, firstbasis, mfl;
	double *rknots, *rcoefs, *rdegrees, *rreplicates;
	double rmin, rmax, temppredict;

	PROTECT(knots = coerceVector(knots, REALSXP));
	PROTECT(replicates = coerceVector(replicates, REALSXP));
/*	PROTECT(min = coerceVector(min, REALSXP));  */
/*	PROTECT(max = coerceVector(max, REALSXP)); */
/*	PROTECT(order = coerceVector(order, INTSXP)); */
	PROTECT(coefs = coerceVector(coefs, REALSXP));
	PROTECT(degrees = coerceVector(degrees, REALSXP));
/*	PROTECT(intercept = coerceVector(intercept, INTSXP)); */
	PROTECT(xvals = coerceVector(xvals, REALSXP));
	PROTECT(xId = coerceVector(xId, INTSXP));
	PROTECT(w	= coerceVector(w, REALSXP));
	PROTECT(fromT = coerceVector(fromT, REALSXP));
	PROTECT(FirstId = coerceVector(FirstId, INTSXP));
	PROTECT(LastId = coerceVector(LastId, INTSXP));
/*	PROTECT(outerok = coerceVector(outerok, LGLSXP)); */
	
	rknots = REAL(knots); 
	nknots = length(knots);
	rreplicates = REAL(replicates);
	theorder = asInteger(order);
	rmin = asReal(min);
	rmax = asReal(max);
	
	rcoefs = REAL(coefs); 
	rdegrees = REAL(degrees); 

	/* number of bases */
	nbases = nknots;
	for( i=0; i<nknots;  i++){
		nbases = nbases + rreplicates[i];
	}
/* first basis to start with : 
 *     if intercept =TRUE, all bases                -> firstbasis = 0, 
 *     if intercept =FALSE, remove base number one  -> firstbasis = 1   */
	ii = asLogical(intercept);
	if(ii == NA_LOGICAL) {
		error("'intercept' must be TRUE or FALSE");    
	} else   if (ii) {
		firstbasis = 0;
	} else {
		firstbasis = 1;
	}
	
	rxvals = REAL(xvals); 
	nx = length(xvals);
	
	rxId = INTEGER(xId);
	nxId = length(xId);
	
	rw = REAL(w); 
	nw = length(w);

	rfromT = REAL(fromT); 
	nt = length(fromT);
	
	rFirstId = INTEGER(FirstId);
	nFId = length(FirstId);

	rLastId = INTEGER(LastId);

	if(nw != nt) {
		error("length of 'W' and 'fromT' differ");    
	} 	
	
	if(nFId != nt) {
		error("length of 'FirstId' and 'fromT' differ");    
	} 	
	
	if(nxId != nx) {
		error("length of 'xId' and 'xvals' differ");    
	} 	
	
	PROTECT(wce = allocVector(REALSXP, nx));
	rwce = REAL(wce);
	
/* check if necessary to set the value of wce at outer xvals	*/
	oo = asLogical(outerok);
	
	if(oo == NA_LOGICAL) {
		error("'outer.ok' must be TRUE or FALSE");    
	} else  if (oo) {
		outer_val = 0.0;
	} else {
		outer_val = R_NaN;
	}

	for(i = 0; i < nx; i++) {
		rwce[i] = 0.0;

		/* starting index */
		ixId = rxId[i];
		istart = rFirstId[ixId] -1;
		ilast  = rLastId[ixId];

		for( j = istart ; j<ilast && rfromT[j] < rxvals[i] ; j++){
			dxval = rxvals[i]-rfromT[j];
			PREDIC_one_trunc_power_basis(dxval, dwce);
			rwce[i] += rw[j] * dwce;
		}
	}

	UNPROTECT(11);
	return(wce);
}



