/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-2012   The R Core Team.
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
 *  http://www.r-project.org/Licenses/
 */

#include <R.h>
#include <Rinternals.h>

#include "flexrsurv.h"

#include <R_ext/Rdynload.h>


#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef callMethods[] = {
		CALLDEF(eval_spline_basis, 6),
		CALLDEF(eval_linex_spline_basis, 9),
		CALLDEF(eval_trunc_power_basis, 10),
		CALLDEF(eval_trunc_power_increasing_basis, 10),
		CALLDEF(eval_lc_spline_basis, 7),
		CALLDEF(eval_lc_linex_spline_basis, 9),
		CALLDEF(eval_lc_trunc_power_basis, 11),
		CALLDEF(eval_lc_trunc_power_increasing_basis, 11),
		CALLDEF(slow_predict_spline_basis, 6),
		CALLDEF(predict_spline_basis, 6),
		CALLDEF(slow_predict_linex_spline_basis, 9),
		CALLDEF(predict_linex_spline_basis, 9),
		CALLDEF(predict_trunc_power_increasing_basis, 10),
		CALLDEF(predict_trunc_power_basis, 10),
		CALLDEF(predict_wce_spline_basis, 11),
		CALLDEF(predict_wce_espline_basis, 10),
		CALLDEF(predict_wce_trunc_power_basis, 15),
		CALLDEF(grad_wce_spline_basis, 11),
		CALLDEF(grad_wce_espline_basis, 10),
		CALLDEF(grad_wce_trunc_power_basis, 15),
		CALLDEF(intTDft_NC, 8),
		CALLDEF(intTDftbase_NC, 10),
		CALLDEF(intTDftwcebase_NC, 10),
		{NULL, NULL, 0}

};


void R_init_flexrsurv(DllInfo *info)
{
	/* Register the .Call routines.
No  .C() .Fortran() or .External() routines,
so pass those arrays as NULL.
	 */
	R_registerRoutines(info, 
			NULL,
			callMethods,
			NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
}
