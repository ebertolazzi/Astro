/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2014
  All Rights Reserved.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation;

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
\****************************************************************************/

#include "Utils_mex.hh"
#include "Astro.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#define MEX_ERROR_MESSAGE \
"=====================================================================================\n" \
"\n" \
"USAGE:\n" \
"    [deltaV,deltaV2] = DeltaV(P0,V0,P1,V1,mu);\n" \
"\n"

#define CHECK_IN(N) \
  UTILS_ASSERT( nrhs == N, CMD "Expected {} argument(s), nrhs = {}\n", N, nrhs )

#define CHECK_OUT(N) \
  UTILS_ASSERT( nlhs == N, CMD "Expected {} argument(s), nlhs = {}\n", N, nlhs )

namespace AstroLib {
  // function [deltaV,deltaV2] = DeltaV(mu,P0,V0,P1,V1)

  using namespace std;

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    try {
      UTILS_ASSERT0( nrhs == 5, "DeltaV: Expected 5 arguments" );
      UTILS_ASSERT0( nlhs == 2, "DeltaV: Expected 2 outputs" );

      real_type mu = Utils::mex_get_scalar_value( arg_in_0, "Lambert first argument (mu) must be a scalar\n" );

      mwSize n;
      double const * PP0 = Utils::mex_vector_pointer( arg_in_1, n, "DeltaV second argument must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt("DeltaV second argument must be a 3d vector\n");

      double const * VV0 = Utils::mex_vector_pointer( arg_in_2, n, "DeltaV third argument must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt("DeltaV third argument must be a 3d vector\n");

      double const * PP1 = Utils::mex_vector_pointer( arg_in_3, n, "DeltaV 4th argument must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt("DeltaV 4th argument must be a 3d vector\n");

      double const * VV1 = Utils::mex_vector_pointer( arg_in_4, n, "DeltaV 5th argument must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt("DeltaV 5th argument must be a 3d vector\n");

      dvec3_t P0, V0, P1, V1;
      std::copy_n( PP0, 3, P0.data() );
      std::copy_n( VV0, 3, V0.data() );
      std::copy_n( PP1, 3, P1.data() );
      std::copy_n( VV1, 3, V1.data() );
/*
  typedef struct {
    bool      long_path;   // true se il minimo è sulla "long path"
    bool      left_branch; // true se il minimo è sul "left branch"
    real_type optimal_travel_time;
    real_type period;
    dvec3_t   W1;
    dvec3_t   W2;
  } minimum_DeltaV_extra;
*/

      real_type * DV  = Utils::mex_create_matrix_value( arg_out_0, 3, 1 );
      real_type * DV2 = Utils::mex_create_matrix_value( arg_out_1, 3, 1 );

      DV[0]  = minimum_DeltaV( mu, P0, V0, P1, V1, DV[1], DV[2], nullptr );
      DV2[0] = minimum_DeltaV2( mu, P0, V0, P1, V1, DV2[1], DV2[2] );

    }
    catch ( std::exception & exc ) {
      mexPrintf("Error: %s\n", exc.what() );
    }
    catch (...) {
      mexPrintf("Unknown error\n");
    }
  }
}
