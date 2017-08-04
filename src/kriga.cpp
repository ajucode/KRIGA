/*__________________________________________________________________________________________________

  author: Pedro Guarderas
  email: pedro.felipe.guarderas@gmail.com
  date: 08-04-2013
  file: kriging.cpp

  This program is free software; you can redistribute it and/or modify it under the 
  terms of the GNU General Public License as published by the Free Software Foundation; 
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
 */

#include "kriging.hpp"
#include <gsl/gsl_math.h>
#include <string>


Kriga::KrigaError::KrigaError( const Exep e, const GA_array& array ) :
_error( e ),
_garray( array ) {
}

const char* Kriga::KrigaError::what() const throw () {
	std::string message;
	switch ( _error ) {
	case type:
		message = "Problems in type agreement with the following Global Arrays:\n";
		break;
	case dimension:
		message = "Problem in dimension agreement with the following Global Arrays:\n";
		break;
	case name:
		message = "Problem in name agreement with the following Global Arrays:\n";
		break;
	}
	for ( GA::GlobalArray* a : _garray ) {
		message = message + " " +
				std::string( a->inquireName() );
	}
	return message.c_str();
}

std::ostream& operator << ( std::ostream& out, const Kriga::KrigaError& Error  ) {
	out << "Error:" << std::endl;
	for ( GA::GlobalArray* a : Error. _garray ) {
		a->printDistribution();
		out << "--------------------------------------------------" << std::endl;
	}
}

/*__________________________________________________________________________________________________
  Covariance matrix constructor based in Global Arrays
 */
GA::GlobalArray* Kriga::covariance( GA::GlobalArray* X,
		GA::GlobalArray* w,
		double (* kernel)( const double& ) ) {
	int n = 1;
	int typew, typeX;
	int ndimw, ndimX;
	int dimsw[1], dimsX[2];

	X->inquire( &typeX, &ndimX, dimsX );
	w->inquire( &typew, &ndimw, dimsw );

	if ( ndimX == 2 &&  ndimw == 1 ) {
		if ( dimsX[1] == dimsw[0] ) {
			int i, j;
			int liX[2], ljX[2], hiX[2], hjX[2];
			int ld[1], hd[1], lD[2], hD[2];
			double alpha, beta, dot = 0;

			int Dims[2] = { dimsX[0], dimsX[0] };
			int dims[1] = { dimsw[0] };

			GA::GlobalArray* D = GA::SERVICES.createGA( typeX, 2, Dims, (char*)"cov", nullptr );
			GA::GlobalArray* d = GA::SERVICES.createGA( typeX, 1, dims, (char*)"d", nullptr );

			alpha = 1.0;
			beta = -1.0;

			ld[0] = 0; hd[0] = dimsX[1] - 1;

			for ( i = 0; i < dimsX[0]; i++ ) {

				liX[0] = i; liX[1] = 0;
				hiX[0] = i; hiX[1] = dimsX[1] - 1;

				for( j = i; j < dimsX[0]; j++ ) {

					ljX[0] = j; ljX[1] = 0;
					hjX[0] = j; hjX[1] = dimsX[1] - 1;

					d->addPatch( &alpha, X, liX, hiX,
							&beta, X, ljX, hjX,
							ld, hd );
					d->elemMultiply( d, d );
					dot = (* kernel)( sqrt( d->ddot( w ) ) );

					lD[0] = i; lD[1] = j;
					hD[0] = i; hD[1] = j;

					D->put( lD, hD, &dot, &n );

					if ( i != j ) {
						lD[0] = j; lD[1] = i;
						hD[0] = j; hD[1] = i;
						D->put( lD, hD, &dot, &n );
					}
				}
			}

			d->destroy();

			return D;
		}
		else {
			GA_array galist = { X, w };
			throw KrigaError( Kriga::KrigaError::dimension, galist );
			return nullptr;
		}
	}
	else {
		GA_array galist = { X, w };
		throw KrigaError( Kriga::KrigaError::type, galist );
		return nullptr;
	}
}


GA::GlobalArray* Kriga::covariance_prediction( GA::GlobalArray* x0,
		GA::GlobalArray* X,
		GA::GlobalArray* w,
		double (* kernel)( const double& ) ) {
	int n = 1;
	int typex0, typeX, typew;
	int ndimx0, ndimX, ndimw;
	int dimsx0[2], dimsX[2], dimsw[1];

	x0->inquire( &typex0, &ndimx0, dimsx0 );
	X->inquire( &typeX, &ndimX, dimsX );
	w->inquire( &typew, &ndimw, dimsw );

	if ( ndimX == 2 &&  ndimw == 1 && ndimx0 >= 1 && ndimx0 <= 2  ) {
		int i, j;
		int lx0[ndimx0], hx0[ndimx0];
		int lX[2], hX[2];
		int ld[1], hd[1];
		double alpha, beta, dot = 0;


		int dims[1] = { dimsw[0] };
		GA::GlobalArray* D;
		GA::GlobalArray* d = GA::SERVICES.createGA( typeX, 1, dims, (char*)"d", nullptr );

		alpha = 1.0;
		beta = -1.0;

		ld[0] = 0; hd[0] = dimsw[0] - 1;

		if ( ndimx0 == 1 ) {
			int lD[1], hD[1];
			int Dims[1] = { dimsX[0] };

			D = GA::SERVICES.createGA( typeX, 1, Dims, (char*)"cov", nullptr );
			lx0[0] = 0; hx0[0] = dimsx0[0] - 1;

			for ( j = 0; j < dimsX[0]; j++ ) {

				lX[0] = j; lX[1] = 0;
				hX[0] = j; hX[1] = dimsX[1] - 1;

				d->addPatch( &alpha, x0, lx0, hx0,
						&beta, X, lX, hX,
						ld, hd );
				d->elemMultiply( d, d );
				dot = (* kernel)( sqrt( d->ddot( w ) ) );

				lD[0] = j;
				hD[0] = j;

				D->put( lD, hD, &dot, &n );
			}
		}
		else if ( ndimx0 == 2 ) {
			int lD[2], hD[2];
			int Dims[2] = { dimsX[0], dimsx0[0] };

			D = GA::SERVICES.createGA( typeX, 2, Dims, (char*)"cov", nullptr );

			for ( i = 0; i < dimsx0[0]; i++ ) {
				lx0[0] = i; lx0[1] = 0;
				hx0[0] = i; hx0[1] = dimsx0[0] - 1;

				for ( j = 0; j < dimsX[0]; j++ ) {

					lX[0] = j; lX[1] = 0;
					hX[0] = j; hX[1] = dimsX[1] - 1;

					d->addPatch( &alpha, x0, lx0, hx0,
							&beta, X, lX, hX,
							ld, hd );
					d->elemMultiply( d, d );
					dot = (* kernel)( sqrt( d->ddot( w ) ) );

					lD[0] = j; lD[1] = i;
					hD[0] = j; hD[1] = i;

					D->put( lD, hD, &dot, &n );

				}
			}
		}
		d->destroy();

		return D;
	}
	else {
		GA_array galist = { x0, X, w };
		throw KrigaError( Kriga::KrigaError::dimension, galist );
		return nullptr;
	}
}

/*__________________________________________________________________________________________________
  Simple Kriging
 */
Kriga::KrigResult Kriga::simple_kriging( GA::GlobalArray* x0,
		GA::GlobalArray* X,
		GA::GlobalArray* Z,
		GA::GlobalArray* w,
		double (* kernel)( const double& ) ) {

	int typex0, typeX, typeZ;
	int ndimx0, ndimX, ndimZ;
	int dimsx0[2], dimsX[2], dimsZ[2];

	x0->inquire( &typex0, &ndimx0, dimsx0 );
	X->inquire( &typeX, &ndimX, dimsX );
	Z->inquire( &typeZ, &ndimZ, dimsZ );

	if ( typex0 == typeX && typeX == typeZ &&
			ndimX == 2 && ndimZ == 2 && ndimx0 >= 1 && ndimx0 <= 2 &&
			dimsX[0] == dimsZ[0] && dimsx0[1] == dimsX[1] )  {

		int m = dimsX[0] - 1;
		int n = dimsx0[0] - 1;
		int p = dimsZ[1] - 1;

		int dimsKk[2] = { m + 1, n + 1 };
		int dimsMSE[2] = { n + 1, n + 1 };
		int dimsz[2] = { n + 1, p + 1 };

		double alpha = 1.0, beta = 0.0, gamma = -1.0, sigma;

		GA::GlobalArray* Kk = GA::SERVICES.createGA( typeX, 2, dimsKk, (char*)"Kk", nullptr );
		GA::GlobalArray* MSE = GA::SERVICES.createGA( typeX, 2, dimsMSE, (char*)"MSE", nullptr );
		GA::GlobalArray* z = GA::SERVICES.createGA( typeZ, 2, dimsz, (char*)"z", nullptr );

		GA::GlobalArray* K = covariance( X, w, kernel );
		GA::GlobalArray* k = covariance_prediction( x0, X, w, kernel );

		// expensive part
		K->spdInvert();

		Kk->matmulPatch( 't', 'n', &alpha, &beta,
				k, 0, m, 0, n,
				K, 0, m, 0, m,
				0, n, 0, m );

		z->matmulPatch( 'n', 'n', &alpha, &beta,
				Kk, 0, n, 0, m,
				Z, 0, m, 0, p,
				0, n, 0, p);

		MSE->matmulPatch( 'n', 'n', &gamma, &beta,
				Kk, 0, m, 0, m,
				k, 0, m, 0, n,
				0, n, 0, n );

		sigma = (* kernel)( 0.0 );
		MSE->addConstant( &sigma );

		Kk->destroy();
		Kriga::KrigResult Result = { z, k, K, MSE, nullptr, Kriga::KrigType::simple };
		return Result;
	}
	else {
		GA_array galist = { x0, X, Z, w };
		throw KrigaError( Kriga::KrigaError::type, galist );
		Kriga::KrigResult Result = { nullptr, nullptr, nullptr, nullptr, nullptr, Kriga::KrigType::simple };
		return Result;
	}
}

/*__________________________________________________________________________________________________
  Ordinary Kriging
 */
Kriga::KrigResult Kriga::ordinary_kriging( GA::GlobalArray* x0,
		GA::GlobalArray* X,
		GA::GlobalArray* Z,
		GA::GlobalArray* w,
		double (* kernel)( const double& ) ) {

	int typex0, typeX, typeZ;
	int ndimx0, ndimX, ndimZ;
	int dimsx0[2], dimsX[2], dimsZ[2];

	x0->inquire( &typex0, &ndimx0, dimsx0 );
	X->inquire( &typeX, &ndimX, dimsX );
	Z->inquire( &typeZ, &ndimZ, dimsZ );

	if ( typex0 == typeX && typeX == typeZ &&
			ndimX == 2 && ndimZ == 2 && ndimx0 >= 1 && ndimx0 <= 2 &&
			dimsX[0] == dimsZ[0] && dimsx0[1] == dimsX[1] )  {

		int m = dimsX[0] - 1;
		int n = dimsx0[0] - 1;
		int p = dimsZ[1] - 1;

		int dimsKk[2] = { m + 1, n + 1 };
		int dimsMSE[2] = { n + 1, n + 1 };
		int dimsz[2] = { n + 1, p + 1 };
		int dimsone[2] = { n + 1, m + 1 };
		int dimsu[2] = { n + 1, n + 1 };

		double alpha = 1.0, beta = 0.0, gamma = -1.0, sigma, sigmau, delta;

		GA::GlobalArray* Kk = GA::SERVICES.createGA( typeX, 2, dimsKk, (char*)"Kk", nullptr );
		GA::GlobalArray* MSE = GA::SERVICES.createGA( typeX, 2, dimsMSE, (char*)"MSE", nullptr );
		GA::GlobalArray* z = GA::SERVICES.createGA( typeZ, 2, dimsz, (char*)"z", nullptr );
		GA::GlobalArray* u = GA::SERVICES.createGA( typeZ, 2, dimsu, (char*)"u", nullptr );
		GA::GlobalArray* ones = GA::SERVICES.createGA( typeZ, 2, dimsone, (char*)"ones", nullptr );

		GA::GlobalArray* K = covariance( X, w, kernel );
		GA::GlobalArray* k = covariance_prediction( x0, X, w, kernel );

		// expensive part
		K->spdInvert();

		ones->fill( &alpha );

		ones->matmulPatch( 't', 'n', &alpha, &beta,
				ones, 0, n, 0, m,
				K, 0, m, 0, m,
				0, n, 0, m );

		sigmau = 1.0 / ones->ddot( ones );

		delta = 1.0;// - ones->( k );

		u->matmulPatch( 'n', 'n', &sigmau, &beta,
				ones, 0, n, 0, m,
				Z, 0, m, 0, p,
				0, n, 0, p );

		Kk->matmulPatch( 't', 'n', &alpha, &beta,
				k, 0, m, 0, n,
				K, 0, m, 0, m,
				0, n, 0, m );

		z->matmulPatch( 'n', 'n', &alpha, &beta,
				Kk, 0, n, 0, m,
				Z, 0, m, 0, p,
				0, n, 0, p );


		MSE->matmulPatch( 'n', 'n', &gamma, &beta,
				Kk, 0, m, 0, m,
				k, 0, m, 0, n,
				0, n, 0, n );

		Z->add( &alpha, z, &delta, u );

		sigma = (* kernel)( 0.0 ) + delta * delta * sigmau;
		MSE->addConstant( &sigma );

		Kk->destroy();
		ones->destroy();
		Kriga::KrigResult Result = { z, k, K, MSE, nullptr, Kriga::KrigType::ordinary };
		return Result;
	}
	else {
		GA_array galist = { x0, X, Z, w };
		throw KrigaError( Kriga::KrigaError::type, galist );
		Kriga::KrigResult Result = { nullptr, nullptr, nullptr, nullptr, nullptr, Kriga::KrigType::ordinary };
		return Result;
	}
}

/*__________________________________________________________________________________________________
  Universal Kriging
 */
Kriga::KrigResult Kriga::universal_kriging( GA::GlobalArray* x0,
		GA::GlobalArray* X,
		GA::GlobalArray* Z,
		GA::GlobalArray* w,
		double (* kernel)( const double& ) ) {

	int typex0, typeX, typeZ;
	int ndimx0, ndimX, ndimZ;
	int dimsx0[2], dimsX[2], dimsZ[2];

	x0->inquire( &typex0, &ndimx0, dimsx0 );
	X->inquire( &typeX, &ndimX, dimsX );
	Z->inquire( &typeZ, &ndimZ, dimsZ );

	if ( typex0 == typeX && typeX == typeZ &&
			ndimX == 2 && ndimZ == 2 && ndimx0 >= 1 && ndimx0 <= 2 &&
			dimsX[0] == dimsZ[0] && dimsx0[1] == dimsX[1] )  {

		int m = dimsX[0] - 1;
		int n = dimsx0[0] - 1;
		int p = dimsZ[1] - 1;

		int dimsKk[2] = { m + 1, n + 1 };
		int dimsMSE[2] = { n + 1, n + 1 };
		int dimsz[2] = { n + 1, p + 1 };
		int dimsone[2] = { n + 1, m + 1 };
		int dimsu[2] = { n + 1, n + 1 };

		double alpha = 1.0, beta = 0.0, gamma = -1.0, sigma, sigmau, delta;

		GA::GlobalArray* Kk = GA::SERVICES.createGA( typeX, 2, dimsKk, (char*)"Kk", nullptr );
		GA::GlobalArray* MSE = GA::SERVICES.createGA( typeX, 2, dimsMSE, (char*)"MSE", nullptr );
		GA::GlobalArray* z = GA::SERVICES.createGA( typeZ, 2, dimsz, (char*)"z", nullptr );
		GA::GlobalArray* u = GA::SERVICES.createGA( typeZ, 2, dimsu, (char*)"u", nullptr );
		GA::GlobalArray* ones = GA::SERVICES.createGA( typeZ, 2, dimsone, (char*)"ones", nullptr );

		GA::GlobalArray* K = covariance( X, w, kernel );
		GA::GlobalArray* k = covariance_prediction( x0, X, w, kernel );

		// expensive part
		K->spdInvert();

		ones->fill( &alpha );

		ones->matmulPatch( 't', 'n', &alpha, &beta,
				ones, 0, n, 0, m,
				K, 0, m, 0, m,
				0, n, 0, m );

		sigmau = 1.0 / ones->ddot( ones );

		delta = 1.0;// - ones->( k );

		u->matmulPatch( 'n', 'n', &sigmau, &beta,
				ones, 0, n, 0, m,
				Z, 0, m, 0, p,
				0, n, 0, p );

		Kk->matmulPatch( 't', 'n', &alpha, &beta,
				k, 0, m, 0, n,
				K, 0, m, 0, m,
				0, n, 0, m );

		z->matmulPatch( 'n', 'n', &alpha, &beta,
				Kk, 0, n, 0, m,
				Z, 0, m, 0, p,
				0, n, 0, p );


		MSE->matmulPatch( 'n', 'n', &gamma, &beta,
				Kk, 0, m, 0, m,
				k, 0, m, 0, n,
				0, n, 0, n );

		Z->add( &alpha, z, &delta, u );

		sigma = (* kernel)( 0.0 ) + delta * delta * sigmau;
		MSE->addConstant( &sigma );

		Kk->destroy();
		ones->destroy();
		Kriga::KrigResult Result = { z, k, K, MSE, nullptr, Kriga::KrigType::ordinary };
		return Result;
	}
	else {
		GA_array galist = { x0, X, Z, w };
		throw KrigaError( Kriga::KrigaError::type, galist );
		Kriga::KrigResult Result = { nullptr, nullptr, nullptr, nullptr, nullptr, Kriga::KrigType::ordinary };
		return Result;
	}
}
