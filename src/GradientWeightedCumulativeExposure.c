/*  Routines for computing the gradient of wheighted cummulative exposure with respect to the spline parameters   
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

SEXP grad_wce_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, 
		SEXP w, SEXP fromT, SEXP FirstId, SEXP LastId,
		SEXP xvals, SEXP xId,
		SEXP outerok)
{
	/* evaluate the gradient of the Weighted Cummulative Exposure with weights defined by non-zero B-spline basis functions at xvals,
	 */
	/*
 knots : vector of ordered unreplicated INTERIOR knots 
Matrices : a vectorized array of dim order X nbases X number_of_intervales(knots) 
  where nbases is the number of bases of the non integrated, non derived splines 
!!!! The Matrix is UNSCALES by the parameters of the spline function  !!!
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
	R_len_t nknots, theorder, nbases,  firstbasis;
	R_len_t i, j, k, l, nw, nt, nx, istart, ilast, ixId, nxId, nFId, oo;
	int *rxId, *rFirstId, *rLastId;
	double *rxvals, dxval, *dgwce, *rw, *rfromT, *rgradwce;
	double *rknots, *rMatrices;
	SEXP dims;
	SEXP gradwce;
	double  outer_val; 
	R_len_t theinterval, mfl;
	double temp, *U, u;

	PROTECT(knots = coerceVector(knots, REALSXP));
	PROTECT(order = coerceVector(order, INTSXP));
	PROTECT(Matrices = coerceVector(Matrices, REALSXP));
	PROTECT(intercept = coerceVector(intercept, INTSXP));
	PROTECT(xvals = coerceVector(xvals, REALSXP));
	PROTECT(xId = coerceVector(xId, INTSXP));
	PROTECT(w	= coerceVector(w, REALSXP));
	PROTECT(fromT = coerceVector(fromT, REALSXP));
	PROTECT(FirstId = coerceVector(FirstId, INTSXP));
	PROTECT(LastId = coerceVector(LastId, INTSXP));
	PROTECT(outerok = coerceVector(outerok, LGLSXP));

	rknots = REAL(knots); 
	nknots = length(knots);
	theorder = INTEGER(order)[0];

	dims = getAttrib(Matrices, R_DimSymbol);
	if( LENGTH(dims) < 3 ){
		error("'Matrices' must be an array with 3 dim");   
	}
	nbases = INTEGER(dims)[1];

	firstbasis = (INTEGER(intercept)[0]==0);
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

	/* check if necessary to set the value of wce at outer xvals */
	oo = asLogical(outerok);

	if(oo == NA_LOGICAL) {
		error("'outer.ok' must be TRUE or FALSE");    
	} else  if (oo) {
		outer_val = 0.0;
	} else {
		outer_val = R_NaN;
	}

	PROTECT(gradwce = allocMatrix(REALSXP, nx, nbases-firstbasis));
	rgradwce = REAL(gradwce);

	dgwce = (double *) R_alloc(nbases - firstbasis , sizeof(double));

	U = (double *) R_alloc( theorder, sizeof(double));

	U[0]=1.0;
	for(i = 0; i < nx; i++) {
		for (k = 0 ; k < nbases - firstbasis; k++) {	
			rgradwce[i + nx *k] = 0.0;
		}

		/* starting index */
		ixId = rxId[i]-1;
		istart = rFirstId[ixId] -1;
		ilast  = rLastId[ixId];
		for( l = istart ; l<ilast && rfromT[l] < rxvals[i] ; l++){
			dxval = rxvals[i]-rfromT[l];
			EVALUATE_one_spline_basis(dxval, dgwce,  ) ;                                  
			/*#define EVALUATE_one_spline_basis(_X_, _P_, _I_)			\                               */
			/*        if (ISNAN(dxval)) {                                                                         */
			/*            for (k = firstbasis; k < nbases; k++) {                                                 */
			/*                dgwce[  (k-firstbasis)] = R_NaN;                                                    */
			/*            }                                                                                       */
			/*        } else {                                                                                    */
			/* find the interval within the range of all the knots (which include boundaries)                   */
			/*   of dxval, rightmost_close=TRUE, all_inside = FALSE                                             */
			/*            if (dxval < rknots[theorder-1] || dxval > rknots[nknots - theorder]  ) {                */
			/* out of the boundary knots interval                                                */
			/*                for (k = firstbasis; k < nbases; k++) {                                             */
			/*                    dgwce[  (k-firstbasis)] = outer_val;                                            */
			/*                }                                                                                   */
			/*            } else {                                                                                */
			/*                mfl = 0;                                                                            */
			/*                theinterval = findInterval2(rknots, nknots, dxval, 1, 0 , FALSE, theorder, &mfl );  */
			/*                if( theinterval > nknots - theorder) {                                              */
			/* xx[i] is the rightmost boundary knot */
			/*                    theinterval = nknots - theorder;                                                */
			/*                }                                                                                   */
			/*                u = (dxval - rknots[theinterval-1])/(rknots[theinterval]-rknots[theinterval-1]);    */
			/*                for ( j = 1; j < theorder ; j++) {                                                  */
			/*                    U[j] = pow(u, (double)j);                                                       */
			/*                }                                                                                   */
			/* the usefull matrix is the (theinterval - theorder +1)th matrix of Matrices */
			/*                theinterval = theinterval - theorder;                                               */
			/*                for (k = firstbasis; k < nbases; k++) {                                             */
			/*                    temp = 0;                                                                       */
			/*                    for (int j = 0; j < theorder ; j++) {                                           */
			/*                        temp += U[j] * rMatrices[theorder*nbases*theinterval+ theorder*k + j];      */
			/*                    }                                                                               */
			/*                    dgwce[  (k-firstbasis)] = temp;                                                 */
			/*                }                                                                                   */
			/*            }                                                                                       */
			/*        }                                                                                           */

			for (k = 0; k < nbases - firstbasis; k++) {	
				rgradwce[i + nx *k ] += rw[l] * dgwce[k];
			}
		}
	}
	UNPROTECT(12);
	return(gradwce);
}

SEXP grad_wce_espline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, 
		SEXP w, SEXP fromT, SEXP FirstId, SEXP LastId,
		SEXP xvals, SEXP xId)
{
	/* evaluate Weighted Cummulative Exposure with weights defined by non-zero B-spline basis functions at xvals,
	 */
	/*
 knots : vector of ordered unreplicated INTERIOR knots 
Matrices : a vectorized array of dim order X nbases X number_of_intervales(knots) 
  where nbases is the number of bases of the non integrated, non derived splines 
!!!! The Matrix is UNSCALED by the parameters of the spline function  !!!
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
	R_len_t nknots, theorder, nbases,  firstbasis;
	R_len_t i, j, k, nw, nt, nx, istart, ilast, ixId, nxId, nFId;
	int *rxId, *rFirstId, *rLastId;
	double *rxvals, dxval, *dgwce, *rw, *rfromT, *rgradwce;
	double *rknots, *rMatrices;
	SEXP dims;
	SEXP gradwce;
	R_len_t theinterval, mfl;
	double temp, *U, u;

	double outer_val;
	outer_val = 0.0;

	PROTECT(knots = coerceVector(knots, REALSXP));
	PROTECT(order = coerceVector(order, INTSXP));
	PROTECT(Matrices = coerceVector(Matrices, REALSXP));
	PROTECT(intercept = coerceVector(intercept, INTSXP));
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

	firstbasis = (INTEGER(intercept)[0]==0);
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

	PROTECT(gradwce = allocMatrix(REALSXP, nx, nbases-firstbasis));
	rgradwce = REAL(gradwce);

	dgwce = (double *) R_alloc(nbases - firstbasis , sizeof(double));

	U = (double *) R_alloc( theorder, sizeof(double));
	U[0]=1.0;
	for(i = 0; i < nx; i++) {
		for (k = 0 ; k < nbases - firstbasis; k++) {	
			rgradwce[i + nx * k] = 0.0;
			/*			dgwce[k] = 0.0; */
		}

		/* starting index */
		ixId = rxId[i]-1;
		istart = rFirstId[ixId]-1;
		ilast  = rLastId[ixId];
		for( j = istart ; j<ilast && rfromT[j] < rxvals[i] ; j++){
			dxval = rxvals[i]-rfromT[j];
			EVALUATE_one_espline_basis(dxval, dgwce, );
			for (k = 0; k < nbases - firstbasis; k++) {	
				rgradwce[i + nx *k ] += rw[j] * dgwce[k];
			}
		}
	}
	unprotect(11);
	return(gradwce);
}

SEXP grad_wce_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
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
!!!! The coefs is UNSCALED by the parameters of the spline function  !!! Spline(x) = sum(coef[i] * monomial_i(x))
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
	R_len_t nw, nt, nx, istart, ilast, ixId, nxId, nFId, oo;
	int *rxId, *rFirstId, *rLastId;
	double *rxvals, dxval, *dgwce, *rw, *rfromT, *rgradwce;
	SEXP gradwce;
	double  outer_val; 

	R_len_t theinterval, firstbasis, mfl;
	double *rknots, *rcoefs, *rreplicates;
	double rmin, rmax;

	PROTECT(knots = coerceVector(knots, REALSXP));
	PROTECT(replicates = coerceVector(replicates, REALSXP));
	PROTECT(min = coerceVector(min, REALSXP));
	PROTECT(max = coerceVector(max, REALSXP));
	PROTECT(order = coerceVector(order, INTSXP));
	PROTECT(coefs = coerceVector(coefs, REALSXP));
	PROTECT(degrees = coerceVector(degrees, REALSXP));
	PROTECT(intercept = coerceVector(intercept, INTSXP));
	PROTECT(xvals = coerceVector(xvals, REALSXP));
	PROTECT(xId = coerceVector(xId, INTSXP));
	PROTECT(w	= coerceVector(w, REALSXP));
	PROTECT(fromT = coerceVector(fromT, REALSXP));
	PROTECT(FirstId = coerceVector(FirstId, INTSXP));
	PROTECT(LastId = coerceVector(LastId, INTSXP));
	PROTECT(outerok = coerceVector(outerok, LGLSXP));

	rknots = REAL(knots); 
	nknots = length(knots);
	rreplicates = REAL(replicates);
	theorder = INTEGER(order)[0];
	rmin = REAL(min)[0];
	rmax = REAL(max)[0];

	rcoefs = REAL(coefs); 

	/* number of bases */
	nbases = nknots;
	for( i=0; i<nknots;  i++){
		nbases = nbases + rreplicates[i];
	}
	firstbasis = (INTEGER(intercept)[0]==0);

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

	/* check if necessary to set the value of wce at outer xvals	*/
	oo = asLogical(outerok);

	if(oo == NA_LOGICAL) {
		error("'outer.ok' must be TRUE or FALSE");    
	} else  if (oo) {
		outer_val = 0.0;
	} else {
		outer_val = R_NaN;
	}

	PROTECT(gradwce = allocMatrix(REALSXP, nx, nbases-firstbasis));
	rgradwce = REAL(gradwce);

	dgwce = (double *) R_alloc(nbases - firstbasis , sizeof(double));

	for(i = 0; i < nx; i++) {
		for (k = 0 ; k < nbases - firstbasis; k++) {	
			rgradwce[i + nx * k] = 0.0;
			/*			dgwce[k] = 0.0; */
		}

		/* starting index */
		ixId = rxId[i];
		istart = rFirstId[ixId] -1;
		ilast  = rLastId[ixId];

		for( j = istart ; j<ilast && rfromT[j] < rxvals[i] ; j++){
			dxval = rxvals[i]-rfromT[j];
			EVALUATE_one_trunc_power_basis(dxval, dgwce, );
			for (k = 0; k < nbases - firstbasis; k++) {	
				rgradwce[i + nx *k ] += rw[j] * dgwce[k];
			}
		}
	}

	UNPROTECT(16);
	return(gradwce);
}




