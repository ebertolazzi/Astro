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

#include "mex_utils.hh"
#include "Astro.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#define MEX_ERROR_MESSAGE \
"=====================================================================================\n" \
"\n" \
"USAGE:\n" \
"    ok = Check_EQ_for_consistency(p,f,g,h,k,L);\n" \
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
  // function ok = Check_EQ_for_consistency(p,f,g,h,k,L)

  using namespace std;

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {
    try {
      MEX_ASSERT( nrhs == 6, "Check_EQ_for_consistency: Expected 6 arguments" );
      MEX_ASSERT( nlhs == 1, "Check_EQ_for_consistency: Expected 1 outputs" );

      // Get values of the scalar inputs
      // Check for the proper type of argument
      real_type p = getScalarValue( arg_in_0, "Lambert first argument (p) must be a scalar\n" );
      real_type f = getScalarValue( arg_in_1, "Lambert 2rd argument (f) must be a scalar\n" );
      real_type g = getScalarValue( arg_in_2, "Lambert 3rd argument (g) must be a scalar\n" );
      real_type h = getScalarValue( arg_in_3, "Lambert 4rd argument (h) must be a scalar\n" );
      real_type k = getScalarValue( arg_in_4, "Lambert 5rd argument (k) must be a scalar\n" );
      real_type L = getScalarValue( arg_in_5, "Lambert 6rd argument (K) must be a scalar\n" );

      bool ok = check_EQ_for_consistency( p, f, g, h, k, L );
      
      setScalarBool(arg_out_0,ok);
    }
    catch ( std::exception & exc ) {
      mexPrintf("Error: %s\n", exc.what() );
    }
    catch (...) {
      mexPrintf("Unknown error\n");
    }
  }
}
