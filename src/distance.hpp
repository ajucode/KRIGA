/*__________________________________________________________________________________________________

  author: Pedro Guarderas
  email: pedro.felipe.guarderas@gmail.com
  date: 02-04-2013
  file: distance.hpp

  This program is free software; you can redistribute it and/or modify it under the 
  terms of the GNU General Public License as published by the Free Software Foundation; 
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#ifndef DISTANCE
#define DISTANCE

#include <gsl/gsl_math.h>

namespace Kriga {

/*__________________________________________________________________________________________________
  Distance functions
*/
/*
  template <typename T, size_t N>
  inline size_t countof(const T (&arr)[N]) {
  return N;
  }
*/
template <typename T, size_t N >
char ( &_ArraySizeHelper( T (&arr)[N] ))[N];
#define countof( arr ) (sizeof( _ArraySizeHelper( arr ) ))


template< size_t N >
double dst_weighted( const double (&x)[N], const double (&y)[N],
		     const double (&w)[N], const double (&p)[N] ) {

    double d = 0.0;
    int size_x = countof( x );
    int size_y = countof( y );
    int size_w = countof( w );
    int size_p = countof( p );
  
    if ( size_x > 0 && size_x == size_y && size_y == size_w && size_w == size_p ) {
	int i;
	for( i = 0; i < size_x ; i++ ) {
	    d += w[i] * pow( abs( x[i] - y[i] ), p[i] );
	}
    }
    return d;
}

double dst_weighted( const int& size, double* x, double* y,
		     double* w, double* p ) {

    double d = 0.0;
    int i;

    for( i = 0; i < size ; i++ ) {
	if ( (x + i) == nullptr || (y + i) == nullptr ||
	     (w + i) == nullptr || (p + i) == nullptr ) {
	    d = 0.0;
	    break;
	}
	else {
	    d += *(w + i) * pow( abs( *(x + i) - *(y + i) ), *(p + i) );
	}
    }

    return d;
}

template< size_t N >
double dst_weighted( double (&x)[N], double (&y)[N],
		     const double& w, const double& p ) {

    double d = 0.0;
    int size_x = countof( x );
    int size_y = countof( y );
  
    if ( size_x > 0 && size_x == size_y ) {
	int i;
	for( i = 0; i < size_x; i++ ) {
	    d += pow( abs( x[i] - y[i] ), p );
	}
    }
    return pow( abs( w * d ), 1.0 / p );
}

double dst_weighted( const int& size, double* x, double* y,
		     const double& w, const double& p ) {

    double d = 0.0;
    int i;
    for( i = 0; i < size; i++ ) {
	if ( (x + i) == nullptr || (y + i) == nullptr ) {
	    d = 0.0;
	    break;
	}
	else {
	    d += pow( abs( *(x + i) - *(y + i) ), p );
	}
    }

    return pow( abs( w * d ), 1.0 / p );
}

} // namespace Kriga

#endif // DISTANCE
