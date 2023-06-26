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
"    [ta,tb,DV] = LambertMinDeltaV(mu,P0,V0,P1,V1,t0,t1,dt);\n" \
"\n"

#define CHECK_IN(N) \
  UTILS_ASSERT( nrhs == N, CMD "Expected {} argument(s), nrhs = {}\n", N, nrhs )

#define CHECK_OUT(N) \
  UTILS_ASSERT( nlhs == N, CMD "Expected {} argument(s), nlhs = {}\n", N, nlhs )

namespace AstroLib {

  using namespace std;

  #define CMD "[ta,tb,DV] = LambertMinDeltaV(mu,P0,V0,P1,V1,t0,t1,dt): "

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    try {
      UTILS_ASSERT0( nrhs == 8, CMD "Expected 8 input arguments" );
      UTILS_ASSERT0( nlhs == 3, CMD "Expected 3 outputs" );

      real_type muS = Utils::mex_get_scalar_value( arg_in_0, CMD "first argument (mu) must be a scalar\n" );

      mwSize n;
      double const * PP0 = Utils::mex_vector_pointer( arg_in_1, n, CMD "second argument must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt(CMD "second argument must be a 3d vector\n");

      double const * VV0 = Utils::mex_vector_pointer( arg_in_2, n, CMD "third argument must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt(CMD "third argument must be a 3d vector\n");

      double const * PP1 = Utils::mex_vector_pointer( arg_in_3, n, CMD "4th argument must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt(CMD "4th argument must be a 3d vector\n");

      double const * VV1 = Utils::mex_vector_pointer( arg_in_4, n, CMD "5th argument must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt(CMD "5th argument must be a 3d vector\n");

      real_type t_begin = Utils::mex_get_scalar_value( arg_in_5, CMD "6th argument must be a scalar\n" );
      real_type t_end   = Utils::mex_get_scalar_value( arg_in_6, CMD "7th argument must be a scalar\n" );
      real_type t_delta = Utils::mex_get_scalar_value( arg_in_7, CMD "8th argument must be a scalar\n" );

      dvec3_t P0, V0, P1, V1;
      std::copy_n( PP0, 3, P0.data() );
      std::copy_n( VV0, 3, V0.data() );
      std::copy_n( PP1, 3, P1.data() );
      std::copy_n( VV1, 3, V1.data() );

      minimum_DeltaV_trip trip;

      // typedef struct {
      //   real_type t_begin;
      //   real_type t_end;
      //   real_type DeltaV1;
      //   real_type DeltaV2;
      // } minimum_DeltaV_trip;

      Astro a_from, a_to;

      bool oka = a_from.setup_using_point_and_velocity( "FROM", P0, V0, muS, t_begin );
      bool okb = a_to.setup_using_point_and_velocity( "TO", P1, V1, muS, t_begin );

      UTILS_ASSERT0( oka, CMD "failed to setup first astro" );
      UTILS_ASSERT0( okb, CMD "failed to setup second astro" );

      bool ok  = minimum_DeltaV( muS, t_begin, t_end, t_delta, a_from, a_to, trip );
      UTILS_ASSERT0( oka, CMD "failed to compute minimum maneuvre" );

      Utils::mex_set_scalar_value( arg_out_0, trip.t_begin );
      Utils::mex_set_scalar_value( arg_out_1, trip.t_end );

      real_type * DV = Utils::mex_create_matrix_value( arg_out_2, 2, 1 );
      DV[0] = trip.DeltaV1;
      DV[1] = trip.DeltaV2;

    }
    catch ( std::exception & exc ) {
      mexPrintf("Error: %s\n", exc.what() );
    }
    catch (...) {
      mexPrintf("Unknown error\n");
    }
  }
}
