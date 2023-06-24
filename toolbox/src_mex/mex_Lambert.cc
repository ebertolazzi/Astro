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
"    [V1,V2,ok] = Lambert(t1,P1,t2,P2,m,mu);\n" \
"\n"

#define CHECK_IN(N) \
  UTILS_ASSERT( nrhs == N, CMD "Expected {} argument(s), nrhs = {}\n", N, nrhs )

#define CHECK_OUT(N) \
  UTILS_ASSERT( nlhs == N, CMD "Expected {} argument(s), nlhs = {}\n", N, nlhs )

namespace AstroLib {
  // function [V1,V2,ok] = Lambert(t1,P1,t2,P2,m,mu)

  using namespace std;

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    try {
      UTILS_ASSERT0( nrhs == 6, "Lambert: Expected 6 arguments" );
      UTILS_ASSERT0( nlhs == 3, "Lambert: Expected 3 outputs" );

      // Get values of the scalar inputs
      // Check for the proper type of argument
      real_type t1 = Utils::mex_get_scalar_value( arg_in_0, "Lambert first argument (t1) must be a scalar\n" );
      real_type t2 = Utils::mex_get_scalar_value( arg_in_2, "Lambert 3rd argument (t2) must be a scalar\n" );

      mwSize n;
      double const * RR1 = Utils::mex_vector_pointer( arg_in_1, n, "Lambert second argument must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt("Lambert second argument must be a 3d vector\n");

      double const * RR2 = Utils::mex_vector_pointer( arg_in_3, n, "Lambert 4th argument must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt("Lambert 4th argument must be a 3d vector\n");

      integer m = Utils::mex_get_int64( arg_in_4, "Lambert 5th argument (m) must be an integer\n" );
      real_type mu = Utils::mex_get_scalar_value( arg_in_5, "Lambert 6th argument (mu) must be a scalar\n" );

      real_type dt = t2-t1;
      UTILS_ASSERT( dt > 0, "lambert, bad t2-t1 = {}\n", dt );

      dvec3_t R1, R2, V1, V2;
      std::copy_n( RR1, 3, R1.data() );
      std::copy_n( RR2, 3, R2.data() );

      int ok = Lambert( R1, R2, dt, m, mu, V1, V2 );

      double * VV1 = Utils::mex_create_matrix_value( arg_out_0, 3, 1 );
      double * VV2 = Utils::mex_create_matrix_value( arg_out_1, 3, 1 );

      std::copy_n( V1.data(), 3, VV1 );
      std::copy_n( V2.data(), 3, VV2 );

      Utils::mex_set_scalar_int64(arg_out_2,ok);

    }
    catch ( std::exception & exc ) {
      mexPrintf("Error: %s\n", exc.what() );
    }
    catch (...) {
      mexPrintf("Unknown error\n");
    }
  }
}
