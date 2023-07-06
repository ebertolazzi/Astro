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
"    [V1,V2,ok] = Lambert(P1,P2,DT,m,mu);\n" \
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

    #define CMD "[V1,V2,ok] = Lambert(P1,P2,DT,m,mu): "

    try {
      UTILS_ASSERT0( nrhs == 5, CMD "Expected 5 arguments" );
      UTILS_ASSERT0( nlhs == 3, CMD "Expected 3 outputs" );

      mwSize n;
      double const * RR1 = Utils::mex_vector_pointer( arg_in_0, n, CMD "P1 must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt( CMD "P1 must be a 3d vector\n" );

      double const * RR2 = Utils::mex_vector_pointer( arg_in_1, n, CMD "P2 must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt( CMD "P2 must be a 3d vector\n" );

      real_type DT = Utils::mex_get_scalar_value( arg_in_2, CMD "DT must be a scalar\n" );

      integer m = Utils::mex_get_int64( arg_in_3, CMD "m must be an integer\n" );
      real_type mu = Utils::mex_get_scalar_value( arg_in_4, CMD "mu must be a scalar\n" );

      UTILS_ASSERT( DT > 0, "lambert, bad DT = {}\n", DT );

      dvec3_t R1, R2, V1, V2;
      std::copy_n( RR1, 3, R1.data() );
      std::copy_n( RR2, 3, R2.data() );

      int ok = Lambert( R1, R2, DT, m, mu, V1, V2 );

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
