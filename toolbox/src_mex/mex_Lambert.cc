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

#include "GenericContainerMatlabInterface.hh"
#include "mex_utils.hh"
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

/*
// redirect stdout, found at
// https://it.mathworks.com/matlabcentral/answers/132527-in-mex-files-where-does-output-to-stdout-and-stderr-go
*/

#ifdef UTILS_OS_LINUX

class mystream : public std::streambuf {
protected:
  virtual
  std::streamsize
  xsputn(const char *s, std::streamsize n) override
  { mexPrintf("%.*s", n, s); mexEvalString("drawnow;"); return n; }

  virtual
  int
  overflow(int c=EOF) override
  { if (c != EOF) { mexPrintf("%.1s", &c); } return 1; }
};

class scoped_redirect_cout {
public:
  scoped_redirect_cout()
  { old_buf = std::cout.rdbuf(); std::cout.rdbuf(&mout); }
  ~scoped_redirect_cout()
  { std::cout.rdbuf(old_buf); }
private:
  mystream mout;
  std::streambuf *old_buf;
};
static scoped_redirect_cout mycout_redirect;

#endif

#define CHECK_IN(N) \
  MEX_ASSERT2( nrhs == N, CMD "Expected {} argument(s), nrhs = {}\n", N, nrhs )

#define CHECK_OUT(N) \
  MEX_ASSERT2( nlhs == N, CMD "Expected {} argument(s), nlhs = {}\n", N, nlhs )

namespace AstroLib {
  // function [V1,V2,ok] = Lambert(t1,P1,t2,P2,m,mu)

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
      real_type t1_day = getScalarValue(
        arg_in_0,
        "Lambert first argument (t1) must be a scalar\n"
      );

      real_type t2_day = getScalarValue(
        arg_in_2,
        "Lambert 3rd argument (t2) must be a scalar\n"
      );

      mwSize n;
      double const * R1_m = getVectorPointer(
        arg_in_1, n, "Lambert second argument must be a 3d vector\n"
      );
      if ( n != 3 )
        mexErrMsgTxt("Lambert second argument must be a 3d vector\n");
        
      double const * R2_m = getVectorPointer(
        arg_in_3, n, "Lambert 4th argument must be a 3d vector\n"
      );
      if ( n != 3 )
        mexErrMsgTxt("Lambert 4th argument must be a 3d vector\n");

      integer m = getInt(
        arg_in_4,
        "Lambert 5th argument (m) must be an integer\n"
      );

      real_type mu_m3_over_s2 = getScalarValue(
        arg_in_5,
        "Lambert 6th argument (mu) must be a scalar\n"
      );

      double * V1_ms = createMatrixValue(arg_out_0,3,1);
      double * V2_ms = createMatrixValue(arg_out_1,3,1);
      
      real_type dt_day = t2_day-t1_day ;
      MEX_ASSERT( dt_day > 0, "lambert, bad t2-t1 = " << dt_day ) ;

      real_type dt_s = oneDay_s * dt_day ;
      int ok = Lambert( R1_m, R2_m, dt_s, m, mu_m3_over_s2, V1_ms, V2_ms ) ;
      
      setScalarInt(arg_out_2,ok);

    }
    catch ( std::exception & exc ) {
      mexPrintf("Error: %s\n", exc.what() ) ;
    }
    catch (...) {
      mexPrintf("Unknown error\n") ;
    }
  }
}
