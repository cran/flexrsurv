/*
 * derived from stats:integrate.c
 */

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-2016  the R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <math.h> // for isfinite
#include <Rinternals.h>
#include <R_ext/Applic.h>
/*--- typedef void integr_fn(double *x, int n, void *ex) ---
 * vectorizing function   f(x[1:n], ...) -> x[]  {overwriting x[]}.
 * Vectorization can be used to speed up the integrand
 * instead of calling it  n  times.
*/


typedef struct int_struct
{
    SEXP f;    /* function */
    SEXP env;  /* where to evaluate the calls */
} int_struct, *IntStruct;


SEXP mkans(double x)
{
// no need for PROTECT() here, as REAL(.) does not allocate:
SEXP ans = allocVector(REALSXP, 1);
REAL(ans)[0] = x;
return ans;
}

SEXP mkians(int x)
{
// no need for PROTECT() here, as REAL(.) does not allocate:
SEXP ans = allocVector(INTSXP, 1);
INTEGER(ans)[0] = x;
return ans;
}

SEXP mkvans(double *x, int n)
{
// no need for PROTECT() here, as REAL(.) does not allocate:
	int i;
SEXP ans = allocVector(REALSXP, n);
    for(i = 0; i < n; i++) REAL(ans)[i] = x[i];
return ans;
}


/* the evaluator f(x) */

/*  with install in env */
double feval(double x, SEXP f, SEXP rho)
{
// a version with (too) much PROTECT()ion .. "better safe than sorry"
SEXP symbol, value;
PROTECT(symbol = install("x"));
PROTECT(value = mkans(x));
defineVar(symbol, value, rho);
UNPROTECT(2);
return(REAL(eval(f, rho))[0]);
}	

/*  vectorised version with lang3*/
double *fveval3(SEXP f, double *x, int i, int n, SEXP rho)
{
// a version with (too) much PROTECT()ion .. "better safe than sorry"
	SEXP value, index;
PROTECT(value = mkvans(x, n));
PROTECT(index = mkians(i));

UNPROTECT(2);
return(REAL(eval(lang3(f, value, index), rho)));
}	


/*  vectorised version with lang2*/
double *fveval2(SEXP f, double *x, int n, SEXP rho)
{
// a version with (too) much PROTECT()ion .. "better safe than sorry"
	SEXP value;
PROTECT(value = mkvans(x, n));

UNPROTECT(1);
return(REAL(eval(lang2(f, value), rho)));
}	



/* computes weights for Newton-Cotes quadrature */
void ncweights(int nstep, double step, int rnctype , double *weights)
/* nstep is the number of subdivisions
   length(weights) = nstep+1  
*/ 
{
	int i;
	if(rnctype == 2){
		/* wehights for cavalieri simpson rule (nstep even) */
		double p1=2.0/3.0 * step;
		double p2=4.0/3.0 * step;
		for(i = 0; i< nstep/2; i++){
			weights[2*i]   = p1;
			weights[2*i+1] = p2;
		}
		weights[0]       =  step/3.0;
		weights[nstep]   =  step/3.0;
		return;
	}
	else 	if(rnctype == 3){
		/* wehights for simpson's 3/8 rule (nstep = 3*ii) */
		double p1=3.0/4.0 * step;
		double p2=9.0/8.0 * step;
		for(i = 0; i< nstep/3; i++){
			weights[3*i]   = p1;
			weights[3*i+1] = p2;
			weights[3*i+2] = p2;
		}
		weights[0]        = 3.0/8.0 * step;
		weights[nstep]   = 3.0/8.0 * step;
		return;
	}
	else	if(rnctype == 4){
		/* wehights for boole rule (nstep = 4*i) */
		double p1=28.0/45.0 * step;
		double p2=64.0/45.0 * step;
		double p3=8.0/15.0 * step;
		for(i = 0; i< nstep/4; i++){
			weights[4*i]   = p1;
			weights[4*i+1] = p2;
			weights[4*i+2] = p3;
			weights[4*i+3] = p2;
		}
		weights[0]       = 14.0/45.0 * step;
		weights[nstep]   = 14.0/45.0 * step;
		return;
	}
}	


SEXP intTDft_NC(SEXP f, SEXP from, SEXP to, SEXP step, SEXP nstep, SEXP nstepmax, SEXP nctype, SEXP rho)
{
int i, j;
int rnstepmax;
double *rfrom, *rto, *rstep, *xwork, *wwork, *fwork, *rres, thestep;
int *rnstep, rnctype, nres;
SEXP res;
PROTECT(from  = coerceVector(from , REALSXP));
PROTECT(to    = coerceVector(to   , REALSXP));
PROTECT(step  = coerceVector(step , REALSXP));
PROTECT(nstep = coerceVector(nstep, INTSXP));

rfrom = REAL(from);
rto = REAL(to); 
rstep = REAL(step);
rnstep = INTEGER(nstep);
/* rdnstep = 1.0 * rnstep; */

	rnctype = asInteger(nctype);
	rnstepmax = asInteger(nstepmax);



	nres = length(from);


	res = PROTECT(allocVector(REALSXP, nres));
	rres = REAL(res);

/* work arrays
 * xwork : for points to evaluate
 * wwork : for weights
 * fwork : for evaluted f(x)
 */
    const void *vmax = vmaxget();

/*    double *rdnstep;                                                */
/*    int method = 2;                                                 */
/*    int nc = 1;                                                     */
/*    R_max_col(rdnstep, &nres, &nc, &Zalphabeta[iT,])imax, &method); */
/*   imax = *std::max_element(rdnstep , rdnstep+nres)
     maxx=rnstep[imax]+1; */
    xwork = (double *) R_alloc((size_t) rnstepmax, sizeof(double));
    wwork = (double *) R_alloc((size_t) rnstepmax, sizeof(double));
    fwork = (double *) R_alloc((size_t) rnstepmax, sizeof(double));

	for( i=0; i<nres;  i++){
/* sets points to evaluate */
		if(rnstep[i] % rnctype) error("inconsistency in the number of subdivisions in the Newton-Cotes formula");
		thestep = rstep[i];

		xwork[0] = rfrom[i];
		for( j = 1; j<rnstep[i]; j++){
/*			xwork[j] = xwork[j-1] + *thestep; */
			xwork[j] = rfrom[i] + j * thestep;
		}
		xwork[rnstep[i]] = rto[i];
/* compute f(xwork) with fveval3/lang3*/
		fwork = fveval3(f, xwork, i+1, rnstep[i]+1, rho);
/* compute weights in wwork*/
		ncweights(rnstep[i], rstep[i], rnctype, wwork);
/* computes and assigne the res */
		rres[i] = 0.0;
		for( j = 0; j<rnstep[i]+1; j++){
			rres[i] += fwork[j] * wwork[j];
		}
/*		rres[i] *= rstep[i]; */
	}
    vmaxset(vmax);
    UNPROTECT(5);
/*  returned value */
return(res);
}


SEXP intTDftbase_NC(SEXP f1, SEXP f2, SEXP from, SEXP to, SEXP step, SEXP nstep, SEXP nstepmax, SEXP nctype,SEXP nbase,  SEXP rho)
{

	int i, j, k;
double *rfrom, *rto, *rstep, *xwork, *wwork, *fwork, *bwork, *tmprres, *rres, thestep;
int *rnstep, rnctype, rnbase, nres, rnstepmax;
SEXP res;
PROTECT(from  = coerceVector(from , REALSXP));
PROTECT(to    = coerceVector(to   , REALSXP));
PROTECT(step  = coerceVector(step , REALSXP));
PROTECT(nstep = coerceVector(nstep, INTSXP));

rfrom = REAL(from);
rto = REAL(to); 
rstep = REAL(step);
rnstep = INTEGER(nstep);
/* rnbase = INTEGER(nbase); */
rnbase = asInteger(nbase);
/* rdnstep = 1.0 * rnstep; */

	rnctype = asInteger(nctype);
	rnstepmax = asInteger(nstepmax);



	nres = length(from);


	res = PROTECT(allocMatrix(REALSXP, nres, rnbase));
	rres = REAL(res);

/* work arrays
 * xwork : for points to evaluate
 * wwork : for weights
 * fwork : for evaluted f(x)
 * bwork : for evaluted splinebase(x)
 */
    const void *vmax = vmaxget();

/*    double *rdnstep;                                                */
/*    int method = 2;                                                 */
/*    int nc = 1;                                                     */
/*    R_max_col(rdnstep, &nres, &nc, &Zalphabeta[iT,])imax, &method); */
/*   imax = *std::max_element(rdnstep , rdnstep+nres)
     maxx=rnstep[imax]+1; */
    xwork = (double *) R_alloc((size_t) rnstepmax, sizeof(double));
    wwork = (double *) R_alloc((size_t) rnstepmax, sizeof(double));
    fwork = (double *) R_alloc((size_t) rnstepmax, sizeof(double));
	for( i=0; i<nres;  i++){
/* sets points to evaluate */
		if(rnstep[i] % rnctype) error("inconsistency in the number of subdivisions in the Newton-Cotes formula");
		thestep = rstep[i];

		xwork[0] = rfrom[i];
		for( j = 1; j<rnstep[i]; j++){
/*			xwork[j] = xwork[j-1] + *thestep; */
			xwork[j] = rfrom[i] + j * thestep;
		}
		xwork[rnstep[i]] = rto[i];
/* compute f1(xwork, i) with fveval3 with mlang3 */
		fwork = fveval3(f1, xwork, i+1, rnstep[i]+1, rho);
/* compute base_matrix_spline(x)  with fveval2 with lang2*/
/* bwork is (rnstep[i]+1)%*%nabse matrix */ 
		bwork = fveval2(f2, xwork, rnstep[i]+1, rho);

/* compute weights in wwork*/
		ncweights(rnstep[i], rstep[i], rnctype, wwork);
/* computes and assigne the res */
		for( k=0; k<rnbase;  k++){
/* point to the colomn of th kth base */
			tmprres = rres + i + k*nres;
			*tmprres = 0.0;
			for( j = 0; j<rnstep[i]+1; j++){
				rres[i + k*nres] += fwork[j] * wwork[j] * bwork[j+k*(rnstep[i]+1)];
			}
		}
	}

    vmaxset(vmax);
    UNPROTECT(5);
/*  returned value */
return(res);
}



SEXP intTDftwcebase_NC(SEXP f1, SEXP f2, SEXP from, SEXP to, SEXP step, SEXP nstep, SEXP nstepmax, SEXP nctype,SEXP nbase,  SEXP rho)
{
	int i, j, k;
int rnstepmax;
double *rfrom, *rto, *rstep, *xwork, *wwork, *fwork, *bwork, *tmpbwork, *tmprres, *rres, *thestep;
int *rnstep, rnctype, *rnbase, nres;
SEXP res;

PROTECT(from  = coerceVector(from , REALSXP));
PROTECT(to    = coerceVector(to   , REALSXP));
PROTECT(step  = coerceVector(step , REALSXP));
PROTECT(nstep = coerceVector(nstep, INTSXP));

rfrom = REAL(from);
rto = REAL(to); 
rstep = REAL(step);
rnstep = INTEGER(nstep);
rnbase = INTEGER(nbase);
/* rdnstep = 1.0 * rnstep; */

	rnctype = asInteger(nctype);
	rnstepmax = asInteger(nstepmax);



	nres = length(from);


	res = PROTECT(allocMatrix(REALSXP, nres, *rnbase));
	rres = REAL(res);

/* work arrays
 * xwork : for points to evaluate
 * wwork : for weights
 * fwork : for evaluted f(x)
 * bwork : for evaluted splinebase(x)
 */
    const void *vmax = vmaxget();

/*    double *rdnstep;                                                */
/*    int method = 2;                                                 */
/*    int nc = 1;                                                     */
/*    R_max_col(rdnstep, &nres, &nc, &Zalphabeta[iT,])imax, &method); */
/*   imax = *std::max_element(rdnstep , rdnstep+nres)
     maxx=rnstep[imax]+1; */
    xwork = (double *) R_alloc((size_t) rnstepmax, sizeof(double));
    wwork = (double *) R_alloc((size_t) rnstepmax, sizeof(double));
    fwork = (double *) R_alloc((size_t) rnstepmax, sizeof(double));

	for( i=0; i<nres;  i++){
/* sets points to evaluate */
		if(rnstep[i] % rnctype) error("inconsistency in the number of subdivisions in the Newton-Cotes formula (row %d ; Nstep %d ; degree %d)", i, rnstep[i], rnctype);
		thestep = rstep+i;

		xwork[0] = rfrom[i];
		for( j = 1; j<rnstep[i]; j++){
			xwork[j] = xwork[j-1] + *thestep;
		}
		xwork[rnstep[i]] = rto[i];
/* compute f1(xwork, i) with fveval3 with mlang3 */
		fwork = fveval3(f1, xwork, i+1, rnstep[i]+1, rho);
/* compute base_matrix_spline(x)  with fveval2 with lang2*/
/* bwork is (rnstep[i]+1)%*%nabse matrix */ 
		bwork = fveval3(f2, xwork, i+1, rnstep[i]+1, rho);
/* compute weights in wwork*/
		ncweights(rnstep[i], rstep[i], rnctype, wwork);
/* computes and assigne the res */
		for( k=0; k<*rnbase;  k++){
/* point to the colomn of th kth base */
			tmprres = rres + i + k*nres;
			tmpbwork = bwork + k*(rnstep[i]+1);
			*tmprres = 0.0;
			for( j = 0; j<rnstep[i]+1; j++){
				*tmprres += fwork[j] * wwork[j] * tmpbwork[j];
			}
/*			*tmprres *= rstep[i]; */
		}
	}
    vmaxset(vmax);
    UNPROTECT(5);
/*  returned value */
return(res);
}






