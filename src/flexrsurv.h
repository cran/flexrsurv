SEXP eval_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, SEXP xvals, SEXP outerok);

SEXP eval_espline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, SEXP xvals, SEXP outerok);

SEXP eval_linex_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP linexinf, SEXP linexsup,  SEXP orderextrapol,
			     SEXP intercept, SEXP xvals, SEXP outerok);

SEXP eval_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, SEXP coefs, 
                            SEXP degrees, SEXP intercept, SEXP xvals, SEXP outerok);

SEXP eval_trunc_power_increasing_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
                                       SEXP coefs, SEXP degrees, SEXP intercept, SEXP xvals, SEXP outerok);

SEXP eval_lc_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, SEXP xvals, SEXP beta, SEXP outerok);

SEXP eval_lc_linex_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP linexinf, SEXP linexsup, 
                                SEXP intercept, SEXP xvals, SEXP beta, SEXP outerok);

SEXP eval_lc_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
					SEXP coefs, SEXP degrees, SEXP intercept, 
					SEXP xvals, SEXP beta, SEXP outerok);

SEXP eval_lc_trunc_power_increasing_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
					SEXP coefs, SEXP degrees, SEXP intercept, 
					SEXP xvals, SEXP beta, SEXP outerok);

SEXP slow_predict_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, SEXP xvals, SEXP outerok);

SEXP predict_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, SEXP xvals, SEXP outerok);

SEXP slow_predict_linex_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP linexinf, SEXP linexsup, SEXP intercept, SEXP xvals, SEXP outerok);

SEXP predict_linex_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP linexinf, SEXP linexsup, SEXP intercept, SEXP xvals, SEXP outerok);


SEXP predict_trunc_power_increasing_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
							SEXP coefs, SEXP degrees, SEXP intercept, 
							SEXP xvals, SEXP outerok);

SEXP predict_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
					SEXP coefs, SEXP degrees, SEXP intercept, 
					SEXP xvals, SEXP outerok);

/* WCE functions */
SEXP eval_wce_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, 
			   SEXP w, SEXP fromT, SEXP FirstId, SEXP LastId, 
			   SEXP xvals, SEXP xId,  
			   SEXP outerok);


SEXP eval_wce_espline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, 
			   SEXP w, SEXP fromT, SEXP FirstId, SEXP LastId,  
			    SEXP xvals, SEXP xId);



SEXP eval_wce_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
				SEXP coefs, SEXP degrees, SEXP intercept, 
				SEXP w, SEXP fromT, SEXP FirstId, SEXP LastId, 
				SEXP xvals, SEXP xId,  
				SEXP outerok);

