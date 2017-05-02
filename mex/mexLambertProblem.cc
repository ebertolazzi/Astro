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

#include "common_Mex.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#define arg0 prhs[0]
#define arg1 prhs[1]
#define arg2 prhs[2]
#define arg3 prhs[3]
#define arg4 prhs[4]
#define arg5 prhs[5]

#define out0 plhs[0]
#define out1 plhs[1]
#define out2 plhs[2]

#define MEX_ASSERT(COND,MSG)  \
  if ( !(COND) ) {            \
    std::ostringstream ost ;  \
    ost << MSG << '\n' ;      \
    mexErrMsgTxt(ost.str().c_str()) ; \
  }

// function [V1,V2,ok] = Lambert(t1,P1,t2,P2,m,mu)
namespace astro {

  using namespace std ;

  static real_type const oneDay_s = 86400 ;

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    try {
      MEX_ASSERT( nrhs == 6, "lambert: Expected 6 arguments" ) ;
      MEX_ASSERT( nlhs == 3, "lambert: Expected 3 outputs" ) ;

      // Get values of the scalar inputs
      // Check for the proper type of argument
      real_type t1_day ;
      if ( !is_scalar(arg0,t1_day) )
        mexErrMsgTxt("Lambert first argument (t1) must be a scalar\n");

      real_type t2_day ;
      if ( !is_scalar(arg2,t2_day) )
        mexErrMsgTxt("Lambert 3rd argument (t2) must be a scalar\n");

      int_type n;
      if ( !is_vector(arg1,n) && n == 3 )
        mexErrMsgTxt("second argument must be a 3d vector\n");

      if ( !is_vector(arg3,n) && n == 3 )
        mexErrMsgTxt("4th argument must be a 3d vector\n");

      double const * R1_m = mxGetPr(arg1);
      double const * R2_m = mxGetPr(arg3);

      int_type m ;
      if ( !is_integer(arg4,m) )
        mexErrMsgTxt("Lambert 5th argument (m) must be an integer\n");

      real_type mu_m3_over_s2 ;
      if ( !is_scalar(arg5,mu_m3_over_s2) )
        mexErrMsgTxt("Lambert 5th argument (mu) must be a scalar\n");

      out0 = mxCreateDoubleMatrix(3,1,mxREAL);
      if ( out0 == nullptr )
        mexErrMsgTxt("Could not create mxArray[out0]\n");

      out1 = mxCreateDoubleMatrix(3,1,mxREAL);
      if ( out1 == nullptr )
        mexErrMsgTxt("Could not create mxArray[out1]\n");

      double * V1_ms = mxGetPr(out0);
      double * V2_ms = mxGetPr(out1);
      
      real_type dt_day = t2_day-t1_day ;
      ASSERT( dt_day > 0, "lambert, bad t2-t1 = " << dt_day ) ;

      real_type dt_s = oneDay_s * dt_day ;
      int ok = lambert( R1_m, R2_m, dt_s, m, mu_m3_over_s2, V1_ms, V2_ms ) ;
      
      to_mxArray(ok, out2);

    }
    catch ( std::exception & exc ) {
      mexPrintf("Error: %s\n", exc.what() ) ;
    }
    catch (...) {
      mexPrintf("Unknown error\n") ;
    }
  }
}
