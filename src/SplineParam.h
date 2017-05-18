/* macros for spline evaluations and prediction                                                   */


#define EVAL_SPLINEPARAM 1

/******************************************************************************************************/
/* evaluation at value _X_, results in _P_                                                        */
/* _P_ is a matrix nx * nbases, evaluation assignement at the ith line                            */
/* _I_ is the expression of the index in _P_ where evaluations are performed                      */
/******************************************************************************************************/

/******************************************************************************************************/
/* for SplineBasis                                                                                */
#define EVALUATE_one_spline_basis(_X_, _P_, _I_)				                    \
        if (ISNAN(_X_)) {                                                                           \
            for (k = firstbasis; k < nbases; k++) {                                                 \
                _P_[ _I_ (k-firstbasis)] = R_NaN;                                                   \
            }                                                                                       \
        } else {                                                                                    \
/* find the interval within the range of all the knots (which include boundaries)                   \
   of _X_, rightmost_close=TRUE, all_inside = FALSE */                                              \
            if (_X_ < rknots[theorder-1] || _X_ > rknots[nknots - theorder]  ) {                    \
/* out of the boundary knots interval                                                */             \
                for (k = firstbasis; k < nbases; k++) {                                             \
                    _P_[ _I_ (k-firstbasis)] = outer_val;                                           \
                }                                                                                   \
            } else {                                                                                \
                mfl = 0;                                                                            \
                theinterval = findInterval2(rknots, nknots, _X_, 1, 0 , FALSE, theorder, &mfl );    \
               if( theinterval > nknots - theorder) {                                               \
                    /* xx[i] is the rightmost boundary knot */                                      \
                    theinterval = nknots - theorder;                                                \
                }                                                                                   \
                u = (_X_ - rknots[theinterval-1])/(rknots[theinterval]-rknots[theinterval-1]);      \
                for ( j = 1; j < theorder ; j++) {                                                  \
                    U[j] = pow(u, (double)j);                                                       \
                }                                                                                   \
                /* the usefull matrix is the (theinterval - theorder +1)th matrix of Matrices */    \
                theinterval = theinterval - theorder;                                               \
                for (k = firstbasis; k < nbases; k++) {                                             \
                    temp = 0;                                                                       \
                    for (int j = 0; j < theorder ; j++) {                                           \
                        temp += U[j] * rMatrices[theorder*nbases*theinterval+ theorder*k + j];      \
                    }                                                                               \
                    _P_[ _I_ (k-firstbasis)] = temp;                                                \
                }                                                                                   \
            }                                                                                       \
        }                                                                                           



/******************************************************************************************************/
/* for extended SplineBasis                                                                        */
#define EVALUATE_one_espline_basis(_X_, _P_, _I_)				\
        if (ISNAN(_X_)) {                                                                           \
            for (j = 0; j < nbases - firstbasis; j++) {                                             \
                _P_[ _I_ j] = R_NaN;                                                             \
            }                                                                                       \
        }                                                                                           \
        else {                                                                                      \
/* find the interval within the range of all the knots (which include boundaries)                   \
   of _X_, rightmost_close=TRUE, all_inside = FALSE */                                              \
            if (_X_ < rknots[theorder-1] ) {                                                         \
 /* before boundary_inf                                               */                            \
                for (j = 0; j < nbases - firstbasis; j++) {                                         \
                    _P_[_I_ j] = outer_val;                                                 \
                }                                                                                   \
            } else {                                                                                \
            if (_X_ > rknots[nknots - theorder]  ) {                                                 \
 /* after boundary_sup                                                */                            \
                    u = (_X_ - rknots[theinterval-1]);                                              \
                    theinterval = nknots - 2 * theorder + 1;                                        \
                    }                                                                               \
                else {                                                                              \
                mfl = 0;                                                                            \
                theinterval = findInterval2(rknots, nknots, _X_, 1, 0 , FALSE, theorder, &mfl );    \
                if( theinterval > nknots - theorder) {                                              \
                /* _X_ is the rightmost boundary knot , last interval*/                             \
                    theinterval = nknots - 2 * theorder;                                            \
                    u=1.0;                                                                          \
                }                                                                                   \
                else {                                                                              \
                    u = (_X_ - rknots[theinterval-1])/(rknots[theinterval]-rknots[theinterval-1]);  \
            /* the usefull matrix is the (theinterval - theorder +1)th matrix of Matrices */        \
            /*            with index (theinterval - theorder)  in the C indexing notation */        \
                    theinterval = theinterval - theorder;                                           \
                }                                                                                   \
                for ( j = 1; j < theorder ; j++) {                                                  \
                    U[j] = pow(u, (double)j);                                                       \
                }                                                                                   \
                for (k = firstbasis; k < nbases; k++) {                                             \
                    temp = 0;                                                                       \
                    for (int j = 0; j < theorder ; j++) {                                           \
                        temp += U[j] * rMatrices[theorder*nbases*theinterval+ theorder*k + j];      \
                    }                                                                               \
                    _P_[_I_ (k-firstbasis)] = temp;                                            \
                }                                                                                   \
		}							\
            }								\
	} 



/******************************************************************************************************/
/* for TruncatedSplineBasis                                                                        */
#define EVALUATE_one_trunc_power_basis(_X_, _P_, _I_)				                        \
    if (ISNAN(_X_)) {											\
	    for (j = 0; j < nbases -firstbasis; j++) {							\
		    _P_[_I_ j] = R_NaN;									\
	    }												\
    } else {												\
	    if (_X_< rmin || _X_ > rmax) {								\
		    for (j = 0; j < nbases - firstbasis; j++) {						\
			    _P_[_I_ j] = outer_val;							\
		    }											\
	    } else {											\
		    theinterval= 1;									\
		    mfl = 0;										\
/* find the interval within interior knots (which exclude boundaries(min max)) of _X_, 			\
   rightmost_close=TRUE, all_inside = FALSE 								\
   if theinterval == 0, xvals[i]<knots[0] first interior knot						\
   if theinterval == nknots , xvals[i]>knots[nknots] last interior knot	   				\
*/ 													\
		    theinterval = findInterval(rknots, nknots, _X_, 1, 0 , theinterval, &mfl );		\
		    /* the first theorder bases are powers of xvals */					\
		    ibase=0;										\
		    for ( j = firstbasis; j < theorder ; j++) {						\
			    _P_[_I_ ibase] = pow(_X_, j)*rcoefs[j];					\
			    ibase++;									\
		    }  											\
/*				ibaseb=theorder-firstbasis; */						\
		    icoef=theorder;									\
		    for (k = 0; k < theinterval; k++) {							\
			    for (j = rreplicates[k] ; j > 0 ; j--) {                                    \
				    _P_[_I_ ibase] = pow((_X_ - rknots[k]), theorder-j)*rcoefs[icoef];	\
				    ibase++;								\
				    icoef++;								\
			    }										\
		    }											\
		    while(ibase < nbases - firstbasis) {						\
			    _P_[_I_ ibase] = 0.0;							\
			    ibase++;									\
		    }											\
	    } 												\
    }

/******************************************************************************************************/
/* for increasing TruncatedSplineBasis                                                                        */
#define EVALUATE_one_trunc_power_increasing_basis(_X_, _P_, _I_)							    \
    if (ISNAN(_X_)) {                                                                                                       \
	    for (j = 0; j < nbases -firstbasis; j++) {									    \
		    _P_[_I_ j] = R_NaN;											    \
	    }														    \
    } else {														    \
	    if (_X_< rmin || _X_ > rmax) {										    \
		    for (j = 0; j < nbases - firstbasis; j++) {								    \
			    _P_[_I_ j] = outer_val;									    \
		    }													    \
	    } else {													    \
		    theinterval= 1;											    \
		    mfl = 0;												    \
           /* find the interval within interior knots (which exclude boundaries(min max)) of _X_, 			    \
	      rightmost_close=TRUE, all_inside = FALSE 									    \
	      if theinterval == 0, xvals[i]<knots[0] first interior knot						    \
	      if theinterval == nknots , xvals[i]>knots[nknots] last interior knot	   				    \
	   */ 														    \
		    theinterval = findInterval(rknots, nknots, _X_, 1, 0 , theinterval, &mfl );				    \
		    /* the first theorder bases are powers of xvals */							    \
		    ibase=0;												    \
		    for ( j = firstbasis; j < theorder ; j++) {								    \
			    _P_[_I_ ibase] = pow(_X_, j)*rcoefs[j];							    \
			    ibase++;											    \
		    }  													    \
/*				ibaseb=theorder-firstbasis; */								    \
		    icoef=theorder;											    \
		    for (k = 0; k < nknots; k++) {									    \
			    if( rknots[k] < 0 ){									    \
				    if( theinterval <= k ){								    \
					    for (j = rreplicates[k]; j > 0 ; j--) {					    \
						    _P_[_I_ ibase] = - pow((rknots[k] - _X_), theorder-j)*rcoefs[icoef];    \
						    ibase++;								    \
						    icoef++;								    \
					    }										    \
				    } else {										    \
					    for (j = rreplicates[k]; j > 0 ; j--) {					    \
						    _P_[_I_ ibase] = 0.0;						    \
						    ibase++;								    \
						    icoef++;								    \
					    }										    \
				    }											    \
			    } else {											    \
				    if( theinterval > k ){								    \
					    for (j = rreplicates[k]; j > 0 ; j--) {					    \
						    _P_[_I_ ibase] = pow((_X_ - rknots[k]), theorder-j)*rcoefs[icoef];	    \
						    ibase++;								    \
						    icoef++;								    \
					    }										    \
				    } else {										    \
					    for (j = rreplicates[k]; j > 0 ; j--) {					    \
						    _P_[_I_ ibase] = 0.0;						    \
						    ibase++;								    \
						    icoef++;								    \
					    }										    \
				    }											    \
			    }												    \
		    }												            \
	    } 												                    \
    }


/******************************************************************************************************/
/* prediction at value _X_, results in _P_                                                        */
/*  _P_ is a scalar                                                                               */
/******************************************************************************************************/

/******************************************************************************************************/
/* for SplineBasis                                                                                */
#define PREDIC_one_spline_basis(_X_, _P_)                                                               \
    if (ISNAN(_X_)) {                                                                                   \
        _P_ = R_NaN;                                                                                    \
    }                                                                                                   \
    else {                                                                                              \
/* find the interval within the range of all the knots (which include boundaries)                       \
of rxvals[i], rightmost_close=TRUE, all_inside = FALSE */                                               \
	    if (_X_ < rknots[theorder-1] || _X_ > rknots[nknots - theorder]  ) {                        \
/* out of the boundary knots interval                                                */                 \
		    _P_ = outer_val;                                                                    \
	    } else {                                                                                    \
		    mfl = 0;                                                                            \
		    theinterval = findInterval2(rknots, nknots, _X_, 1, 0 , FALSE, theorder, &mfl );    \
		    if( theinterval > nknots - theorder) {                                              \
                /* xx[i] is the rightmost boundary knot */                                              \
			    theinterval = nknots - theorder;                                            \
		    }                                                                                   \
		    u = (_X_ - rknots[theinterval-1])/(rknots[theinterval]-rknots[theinterval-1]);      \
            /* the usefull matrix is the (theinterval - theorder +1)th matrix of Matrices */            \
            /*            with index (theinterval - theorder)  in the C indexing notation */            \
		    theinterval = theinterval - theorder;                                               \
		    temppredict = rAddMatrices[theorder*theinterval];                                   \
		    U = 1.0;                                                                            \
		    for (int j = 1; j < theorder ; j++) {                                               \
			    U *= u;                                                                     \
			    temppredict += U * rAddMatrices[theorder*theinterval + j];                  \
		    }                                                                                   \
		    _P_ = temppredict;                                                                  \
	    }                                                                                           \
    }                                                                                                   



/******************************************************************************************************/
/* for extended SplineBasis                                                                           */
#define PREDIC_one_espline_basis(_X_, _P_)                                                                   \
    if (ISNAN(_X_)) {                                                                               	     \
        _P_ = R_NaN;                                                                                	     \
    }                                                                                               	     \
    else if (_X_ < rknots[theorder-1] ) {                                                                     \
 /* before boundary_inf                                               */                            	     \
	    _P_ = outer_val;                                                                               	     \
    } else {                                                                                          	     \
/* find the interval within the range of all the knots (which include boundaries)                      	     \
   of rxvals[i], rightmost_close=TRUE, all_inside = FALSE */                                        	     \
	    if (_X_ > rknots[nknots - theorder]  ) {                                                   	     \
 /* after boundary_sup                                                */                            	     \
		    u = (_X_ - rknots[nknots - theorder]);                                                   \
		    theinterval = nknots - 2 * theorder + 1 ;                                                \
	    } else {                                                                                         \
		    mfl = 0;                                                                                 \
		    theinterval = findInterval2(rknots, nknots, _X_, 1, 0 , FALSE, theorder, &mfl );         \
		    if( theinterval > nknots - theorder) {                                                   \
                /* _X_ is the rightmost boundary knot , last interval*/                             	     \
			    theinterval = nknots - 2 * theorder;                                             \
			    u=1.0;                                                                           \
		    } else {                                                                                 \
			    u = (_X_ - rknots[theinterval-1])/(rknots[theinterval]-rknots[theinterval-1]);   \
            /* the usefull matrix is the (theinterval - theorder +1)th matrix of Matrices */        	     \
            /*            with index (theinterval - theorder)  in the C indexing notation */        	     \
			    theinterval = theinterval - theorder;                                            \
		    }                                                                                        \
	    }                                                                             		     \
	    temppredict = rAddMatrices[theorder*theinterval];                                   	     \
            U = 1.0;                                                                                	     \
            for (int j = 1; j < theorder ; j++) {                                                   	     \
		    U *= u;                                                                                  \
		    temppredict += U * rAddMatrices[theorder*theinterval + j];                               \
            }                                                                                       	     \
            _P_ = temppredict;                                                                      	     \
    }                                                                                           	     





/******************************************************************************************************/
/* for truncated power basis                                                                          */
#define PREDIC_one_trunc_power_basis(_X_, _P_)                                                                     \
    if (ISNAN(_X_)) {                                                                                              \
	    _P_ = R_NaN;											   \
    } else if (_X_< rmin || _X_ > rmax) {									   \
	    _P_ = outer_val;											   \
    } else {													   \
	    theinterval= 1;											   \
	    mfl = 0;												   \
														   \
/* find the interval within interior knots (which exclude boundaries(min max)) of _X_, 				   \
	rightmost_close=TRUE, all_inside = FALSE 								   \
	if theinterval == 0, xvals[i]<knots[0] first interior knot						   \
	if theinterval == nknots , xvals[i]>knots[nknots] last interior knot					   \
*/ 														   \
	    theinterval = findInterval2(rknots, nknots, _X_, 1, 0 , FALSE, theinterval, &mfl ); 		   \
	    icoef=0;												   \
	    temppredict = 0;											   \
	    /* the first theorder bases are powers of xvals */							   \
	    for ( j = firstbasis; j < theorder ; j++) {								   \
		    temppredict = temppredict + pow(_X_, j) * rcoefs[j];					   \
		    icoef++; 											   \
	    }  													   \
	    ibase=0;												   \
	    for (k = 0; k < theinterval; k++) {									   \
		    for (j = rreplicates[k]; j > 0 ; j--) {							   \
			    temppredict = temppredict + pow((_X_ - rknots[k]), rdegrees[ibase]) * rcoefs[icoef];   \
			    ibase++;										   \
			    icoef++;										   \
		    }												   \
	    }													   \
	    _P_ = temppredict;											   \
    } 
		


/******************************************************************************************************/
/* for truncated extended power basis                                                                          */
#define PREDIC_one_etrunc_power_basis(_X_, _P_)                                                                                \
    if (ISNAN(_X_)) {													       \
	    _P_ = R_NaN;												       \
    } else if (_X_< rmin ) {												       \
	    _P_ = 0.0;													       \
    } else {                                                                                          	     		       \
/* find the interval within the range of all the knots (which include boundaries)                      	     		       \
   of rxvals[i], rightmost_close=TRUE, all_inside = FALSE                                   */                                 \
	    if (_X_ > rmax){												       \
		    													       \
	    } else {													       \
		    theinterval= 1;											       \
		    mfl = 0;												       \
	    														       \
/* find the interval within interior knots (which exclude boundaries(min max)) of _X_, 					       \
	rightmost_close=TRUE, all_inside = FALSE 									       \
	if theinterval == 0, xvals[i]<knots[0] first interior knot							       \
	if theinterval == nknots , xvals[i]>knots[nknots] last interior knot						       \
*/ 															       \
		    theinterval = findInterval2(rknots, nknots, _X_, 1, 0 , FALSE, theinterval, &mfl ); 		       \
		    icoef=0;												       \
		    temppredict = 0;											       \
		    /* the first theorder bases are powers of xvals */							       \
		    for ( j = firstbasis; j < theorder ; j++) {								       \
			    temppredict = temppredict + pow(_X_, j) * rcoefs[j];					       \
			    icoef++; 											       \
		    }  													       \
		    ibase=0;												       \
		    for (k = 0; k < theinterval; k++) {									       \
			    for (j = rreplicates[k]; j > 0 ; j--) {							       \
				    temppredict = temppredict + pow((_X_ - rknots[k]), rdegrees[ibase]) * rcoefs[icoef];       \
				    ibase++;										       \
				    icoef++;										       \
			    }												       \
		    }													       \
		    _P_ = temppredict;											       \
	    } 														       
