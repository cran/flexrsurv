SEXP eval_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, SEXP xvals, SEXP outerok);

SEXP eval_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, SEXP coefs, SEXP degrees, SEXP intercept, SEXP xvals, SEXP outerok);

SEXP eval_trunc_power_increasing_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, SEXP coefs, SEXP degrees, SEXP intercept, SEXP xvals, SEXP outerok);

SEXP eval_lc_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, SEXP xvals, SEXP beta, SEXP outerok);

SEXP eval_lc_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, SEXP coefs, SEXP degrees, SEXP intercept, SEXP xvals, SEXP beta, SEXP outerok);

SEXP eval_lc_trunc_power_increasing_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
					SEXP coefs, SEXP degrees, SEXP intercept, 
					SEXP xvals, SEXP beta, SEXP outerok);

SEXP predict_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, SEXP xvals, SEXP outerok);

SEXP predict_trunc_power_increasing_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
							SEXP coefs, SEXP degrees, SEXP intercept, 
							SEXP xvals, SEXP outerok);

SEXP predict_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
					SEXP coefs, SEXP degrees, SEXP intercept, 
					SEXP xvals, SEXP outerok);

SEXP eval_wce_spline_basis(SEXP knots, SEXP order, SEXP Matrices, SEXP intercept, 
						   SEXP xvals, SEXP beta, SEXP w, SEXP fromT, SEXP outerok);

SEXP eval_wce_trunc_power_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
								SEXP coefs, SEXP degrees, SEXP intercept, 
								SEXP xvals, SEXP beta, SEXP w, SEXP fromT, SEXP outerok);

SEXP eval_wce_trunc_power_increasing_basis(SEXP knots, SEXP replicates, SEXP min, SEXP max, SEXP order, 
											SEXP coefs, SEXP degrees, SEXP intercept, 
											SEXP xvals, SEXP beta, SEXP w, SEXP fromT, SEXP outerok);

