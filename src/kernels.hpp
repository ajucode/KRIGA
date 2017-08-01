/*__________________________________________________________________________________________________

  author: Pedro Guarderas
  email: pedro.felipe.guarderas@gmail.com
  date: 02-04-2013
  file: kernels.hpp

  This program is free software; you can redistribute it and/or modify it under the 
  terms of the GNU General Public License as published by the Free Software Foundation; 
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#ifndef KERNELS
#define KERNELS

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_math.h>


namespace Kriga {
/*__________________________________________________________________________________________________
  Isotropic covariance models
*/
// Linear
inline double kf_lin( const double& h, 
												const double& alpha = 1.0 ) {
    return alpha * h;
}

// Square
inline double kf_sqr( const double& h, 
		      const double& alpha = 1.0 ) {
    return alpha * h * h;
}
		       
// Triangular
inline double kf_tri( const double& h,
		      const double& c = 1.0,
		      const double& alpha = 1.0 ) {
    return c * GSL_MAX_DBL( alpha - h, 0 );
}

// Exponential
inline double kf_exp( const double& h,
		      const double& sigma = 1.0,
		      const double& theta = 1.0 ) {
    return sigma * sigma * gsl_sf_exp( -h / theta );
}

// Square exponential ( a.k.a Gaussian )
inline double kf_sqexp( const double& h,
			const double& sigma = 1.0,
			const double& theta = 1.0 ) {
    return sigma * sigma * gsl_sf_exp( -h * h / ( theta * theta ) );
}

// Matérn
inline double kf_matern( const double& h,
			 const double& v = 2.0,
			 const double& sigma = 1.0,
			 const double& theta = 1.0 ) {
    double H = sqrt( v ) * h / theta;
    if ( H > 0 ) {
	return sigma * sigma * 2 * pow( H, v ) * gsl_sf_bessel_Knu( v, 2 * H ) / 
	    gsl_sf_gamma( v );
    }
    else {
	return 1.0;
    }
}

// Multilog
inline double kf_multilog( const double& h,
			   const double& R = 1.0 ) {
    return gsl_sf_log( h * h + R * R );
}

// Natural cubic slpine
inline double kf_cubspline( const double& h,
			    const double& R = 1.0 ) {
    return pow( h * h + R * R, 2.0 / 3.0 );
}

// Thin plate spline
inline double kf_tpspline( const double& h,
			   const double& R = 1.0 ) {
    return ( h * h + R * R ) * gsl_sf_log( h * h + R * R );
}

// 
inline double kf_mix( const double& h,
		      const double& sigma = 1.0, 
		      const double& theta = 1.0 ) {
    double ht = h / theta;
    return sigma * sigma * ( 1 + ( sqrt(5.0) + 5.0 / 3.0 * ht ) * ht ) * 
	gsl_sf_exp( -sqrt(5.0) * ht );
}

} // namespace Kriga

#endif // KERNELS
