/*__________________________________________________________________________________________________

  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 02-04-2013
  file: main.cpp

  This program is free software; you can redistribute it and/or modify it under the 
  terms of the GNU General Public License as published by the Free Software Foundation; 
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#include <ga++.h>
#include <iostream>
#include "kriging.hpp"
#include "kernels.hpp"
#include <mpi.h>
#include <type_traits>

using namespace std;

inline double kernel( const double& h ) {
    return krig::kf_matern( h, 2.0, 1.0, 1.0 );
}

/*inline double dist( const int& size, double* x, double* y ) {
  return krig::dst_weighted( size, x, y, 1.0, 2.0 );
  }*/

int main( int argc, char* argv[] ) {

    int N = 2;
    int n = 3;
    int dimX[2] = {2,3};
    int li1[2] = {0,0};
    int hi1[2] = {0,2};
    int li2[2] = {1,0};
    int hi2[2] = {1,2};
    int dimw[1] = {3};
    int dimx0[1] = {3};
  
    double val = 1.0;
    double vx0 = 1.0;
    double r1[] = { 1.0, 0.0, 0.0 };
    double r2[] = { 0.0, 0.0, 1.0 };

    //  int Dim[3] = {2,3,4};
    //  MPI::Init( argc, argv );
    GA::Initialize( argc, argv, 30000, 30000, MT_DBL, 0 );

    GA::GlobalArray* X = GA::SERVICES.createGA( MT_DBL, N, dimX, (char*)"X", nullptr );
    GA::GlobalArray* w = GA::SERVICES.createGA( MT_DBL, 1, dimw, (char*)"w", nullptr );
    GA::GlobalArray* x0 = GA::SERVICES.createGA( MT_DBL, 1, dimx0, (char*)"x0", nullptr );
    GA::GlobalArray* K;

    X->put( li1, hi1, r1, &n ); 
    X->put( li2, hi2, r2, &n ); 
    w->fill( &val );
    x0->fill( &vx0 );

    try {
	K = krig::covariance( X, w, &kernel );
    
    }
    catch ( krig::KrigError& error ) {
	cout << error.what() << endl;
    }

    X->printDistribution();
    X->print();
    
    w->printDistribution();
    w->print();
    
    K->printDistribution();
    K->print();
    
    X->destroy();
    w->destroy();
    K->destroy();

    //  cout << "rank: "<<rank< decltype( r1 ) >::value  << endl;
    //  cout << "extent: " << extent< decltype( r1 ) >::value  << endl;
    GA::Terminate();
    //  MPI::Finalize();

}
