/*__________________________________________________________________________________________________

  author: Pedro Guarderas
  email: pedro.felipe.guarderas@gmail.com
  date: 02-04-2013
  file: kriging.hpp

  This program is free software; you can redistribute it and/or modify it under the 
  terms of the GNU General Public License as published by the Free Software Foundation; 
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#ifndef KRIGA
#define KRIGA

#include <exception>
#include <ga++.h>
#include <iostream>
#include <list>

namespace Kriga {

  typedef std::list< GA::GlobalArray* > GA_array;

/*__________________________________________________________________________________________________
  Kriging exception class
*/
  class KrigaError {
  public:
    typedef enum Exeptions {
      type,
      dimension,
      name
    } Exep;
      
    KrigaError( const Exep e, const GA_array& array );
      
    const char* what() const throw ();

    const Exep _error;
    const GA_array _garray;
  };

  std::ostream& operator << ( std::ostream& out, const KrigaError& Error  );


/*__________________________________________________________________________________________________
   Covariance matrix constructor based in Global Arrays
*/
  GA::GlobalArray* covariance( GA::GlobalArray* X, 
                               GA::GlobalArray* w,
                               double (* kernel)( const double& ) );

  /* Convariance
     w is used to calculate the weigthed distance function
     distance( x, y ) = sum_i w[i] * ( x[i] - y[i] )^2
     return:
     k0 = [ kernel( x0[.][ j ], X[.][ i ] ) ]
     where i = rows, j = columns
  */
  GA::GlobalArray* covariance_prediction( GA::GlobalArray* x0, 
                                          GA::GlobalArray* X, 
                                          GA::GlobalArray* w,
                                          double (* kernel)( const double& ) );
          

/*__________________________________________________________________________________________________
  Kringing result structure
*/
  typedef enum KrigingType {
    simple,
    ordinary,
    universal
  } KrigType;
    
  typedef struct KrigingResult {
    GA::GlobalArray* z;
    GA::GlobalArray* k;
    GA::GlobalArray* K;
    GA::GlobalArray* MSE;
    GA::GlobalArray* F;
    KrigType type;
  } KrigResult;

/*__________________________________________________________________________________________________
  Simple Kriging
*/
  KrigResult simple_kriging( GA::GlobalArray* x0, 
                             GA::GlobalArray* X,
                             GA::GlobalArray* Z,
                             GA::GlobalArray* w,
                             double (* kernel)( const double& ) );

/*__________________________________________________________________________________________________
  Ordinary Kriging
*/
  KrigResult ordinary_kriging( GA::GlobalArray* x0,
                               GA::GlobalArray* X,
                               GA::GlobalArray* Z,
                               GA::GlobalArray* w,
                               double (* kernel)( const double& ) );

/*__________________________________________________________________________________________________
  Universal Kriging
*/
  KrigResult universal_kriging( GA::GlobalArray* x0,
                                GA::GlobalArray* X,
                                GA::GlobalArray* Z,
                                GA::GlobalArray* w,
                                double (* kernel)( const double& ) );


} // namespace Kriga

#endif // KRIGA
