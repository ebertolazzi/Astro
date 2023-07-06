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
"    [t_begin,t_end,P1,V1,W1,P2,V2,W2,DV1,DV2] = globalMinimumDeltaV(mu,t0,P0,V0,P1,V1,t_begin,t_end,t_tolerance,maxDV);\n" \
"\n"

#define CHECK_IN(N) \
  UTILS_ASSERT( nrhs == N, CMD "Expected {} argument(s), nrhs = {}\n", N, nrhs )

#define CHECK_OUT(N) \
  UTILS_ASSERT( nlhs == N, CMD "Expected {} argument(s), nlhs = {}\n", N, nlhs )

namespace AstroLib {
  // function [deltaV,deltaV2] = DeltaV(mu,P0,V0,P1,V1)

  using namespace std;

  #define CMD "[t_begin,t_end,P1,V1,W1,P2,V2,W2,DV1,DV2] = globalMinimumDeltaV(mu,t0,P0,V0,P1,V1,t_begin,t_end,t_tolerance): "

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    try {
      UTILS_ASSERT0( nrhs == 10, CMD "Expected 10 arguments" );
      UTILS_ASSERT0( nlhs == 10, CMD "Expected 10 outputs" );

      real_type mu = Utils::mex_get_scalar_value( arg_in_0, CMD "mu must be a scalar\n" );
      real_type t0 = Utils::mex_get_scalar_value( arg_in_1, CMD "mu must be a scalar\n" );

      mwSize n;
      double const * PP0 = Utils::mex_vector_pointer( arg_in_2, n, CMD " P0 must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt( CMD " P0 must be a 3d vector\n" );

      double const * VV0 = Utils::mex_vector_pointer( arg_in_3, n, CMD "V0 must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt( CMD "V0 must be a 3d vector\n" );

      double const * PP1 = Utils::mex_vector_pointer( arg_in_4, n, CMD "P1 must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt( CMD "P1 must be a 3d vector\n" );

      double const * VV1 = Utils::mex_vector_pointer( arg_in_5, n,  CMD "V1 must be a 3d vector\n" );
      if ( n != 3 ) mexErrMsgTxt( CMD "V1 must be a 3d vector\n" );

      real_type t_begin         = Utils::mex_get_scalar_value( arg_in_6, CMD "t_begin must be a scalar\n" );
      real_type t_end           = Utils::mex_get_scalar_value( arg_in_7, CMD "t_end must be a scalar\n" );
      real_type t_tolerance     = Utils::mex_get_scalar_value( arg_in_8, CMD "t_tolerance must be a scalar\n" );
      real_type max_accepted_DV = Utils::mex_get_scalar_value( arg_in_9, CMD "max_accepted_DV must be a scalar\n" );

      dvec3_t P0, V0, P1, V1;
      std::copy_n( PP0, 3, P0.data() ); std::copy_n( VV0, 3, V0.data() );
      std::copy_n( PP1, 3, P1.data() ); std::copy_n( VV1, 3, V1.data() );

      Astro a_from, a_to;
      bool ok1 = a_from.setup_using_point_and_velocity( P0, V0, mu, t0 );
      bool ok2 = a_to.setup_using_point_and_velocity( P1, V1, mu, t0 );

      vector<minimum_DeltaV_trip> trips;
      global_minimum_DeltaV( a_from, a_to, t_begin, t_end, t_tolerance, max_accepted_DV, trips, nullptr );

      integer N = trips.size();

      double * v_t_begin = Utils::mex_create_matrix_value( arg_out_0, N, 1 );
      double * v_t_end   = Utils::mex_create_matrix_value( arg_out_1, N, 1 );
      double * v_P1      = Utils::mex_create_matrix_value( arg_out_2, 3, N );
      double * v_V1      = Utils::mex_create_matrix_value( arg_out_3, 3, N );
      double * v_W1      = Utils::mex_create_matrix_value( arg_out_4, 3, N );
      double * v_P2      = Utils::mex_create_matrix_value( arg_out_5, 3, N );
      double * v_V2      = Utils::mex_create_matrix_value( arg_out_6, 3, N );
      double * v_W2      = Utils::mex_create_matrix_value( arg_out_7, 3, N );
      double * v_DV1     = Utils::mex_create_matrix_value( arg_out_8, N, 1 );
      double * v_DV2     = Utils::mex_create_matrix_value( arg_out_9, N, 1 );

      for ( minimum_DeltaV_trip const & S : trips ) {
        *v_t_begin++ = S.t_begin;
        *v_t_end++   = S.t_end;

        *v_P1++ = S.P1.coeff(0);
        *v_P1++ = S.P1.coeff(1);
        *v_P1++ = S.P1.coeff(2);

        *v_V1++ = S.V1.coeff(0);
        *v_V1++ = S.V1.coeff(1);
        *v_V1++ = S.V1.coeff(2);

        *v_W1++ = S.W1.coeff(0);
        *v_W1++ = S.W1.coeff(1);
        *v_W1++ = S.W1.coeff(2);

        *v_P2++ = S.P2.coeff(0);
        *v_P2++ = S.P2.coeff(1);
        *v_P2++ = S.P2.coeff(2);

        *v_V2++ = S.V2.coeff(0);
        *v_V2++ = S.V2.coeff(1);
        *v_V2++ = S.V2.coeff(2);

        *v_W2++ = S.W2.coeff(0);
        *v_W2++ = S.W2.coeff(1);
        *v_W2++ = S.W2.coeff(2);

        *v_DV1++ = S.DeltaV1;
        *v_DV2++ = S.DeltaV2;
      }
    } catch ( std::exception & exc ) {
      mexPrintf("Error: %s\n", exc.what() );
    } catch (...) {
      mexPrintf("Unknown error\n");
    }
  }
}
