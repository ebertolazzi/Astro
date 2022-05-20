/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2019
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Utils_mex.hh"
#include "Astro.hh"
#include "GenericContainerMatlabInterface.hh"


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

#define CHECK_IN(N)  UTILS_ASSERT( nrhs == N, CMD ": expected {} argument(s), nrhs = {}\n", N, nrhs )
#define CHECK_OUT(N) UTILS_ASSERT( nlhs == N, CMD ": expected {} argument(s), nlhs = {}\n", N, nlhs )

namespace AstroLib {

  /*\
   *                      _____                 _   _
   *  _ __ ___   _____  _|  ___|   _ _ __   ___| |_(_) ___  _ __
   * | '_ ` _ \ / _ \ \/ / |_ | | | | '_ \ / __| __| |/ _ \| '_ \
   * | | | | | |  __/>  <|  _|| |_| | | | | (__| |_| | (_) | | | |
   * |_| |_| |_|\___/_/\_\_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|
   *
  \*/

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //   _ __   _____      __
  //  | '_ \ / _ \ \ /\ / /
  //  | | | |  __/\ V  V /
  //  |_| |_|\___| \_/\_/
  */
  static
  void
  do_new(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_1 "AstroMexWrapper('new'[,'name'])"
    #define CMD MEX_ERROR_MESSAGE_1

    UTILS_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    UTILS_ASSERT( nrhs == 1 || nrhs == 2, CMD ": expected 1 or 2 input, nrhs = {}\n", nrhs );
    Astro * ptr = nullptr;
    if ( nrhs == 1 ) {
      ptr = new Astro();
    } else {
      UTILS_ASSERT0( mxIsChar(arg_in_1), CMD ": second argument must be a string" );
      string name = mxArrayToString(arg_in_1);
      ptr = new Astro(name);
    }
    arg_out_0 = Utils::mex_convert_ptr_to_mx<Astro>(ptr);
    #undef CMD
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //       _      _      _
  //    __| | ___| | ___| |_ ___
  //   / _` |/ _ \ |/ _ \ __/ _ \
  //  | (_| |  __/ |  __/ ||  __/
  //   \__,_|\___|_|\___|\__\___|
  */
  static
  void
  do_delete(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_2 "AstroMexWrapper('delete',obj)"
    #define CMD MEX_ERROR_MESSAGE_2
    CHECK_IN(2);
    CHECK_OUT(0);
    Utils::mex_destroy_object<Astro>( arg_in_1 );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_copy(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define MEX_ERROR_MESSAGE_3 "AstroMexWrapper('copy',obj)"
    #define CMD MEX_ERROR_MESSAGE_3

    CHECK_IN(2);
    CHECK_OUT(1);

    Astro const * ptr = Utils::mex_convert_mx_to_ptr<Astro>(arg_in_1);
    arg_out_0 = Utils::mex_convert_ptr_to_mx<Astro>(new Astro( *ptr ));
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_name(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_4 "res = AstroMexWrapper('do_name',obj)"
    #define CMD MEX_ERROR_MESSAGE_4

    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::to_mxArray(ptr->name(),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_check_for_consistency(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_5 "ok = AstroMexWrapper('check_for_consistency',obj)"
    #define CMD MEX_ERROR_MESSAGE_5

    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::to_mxArray(ptr->check_for_consistency(),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_position(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_6 "[x,y,z] = AstroMexWrapper('position',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_6

    CHECK_IN(3);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    if ( nlhs == 1 ) {
      GC_namespace::mat_real_type xyz;
      xyz.resize(3,sz);
      for ( mwSize i = 0 ; i < sz ; ++i ) ptr->position( t[i], xyz(0,i), xyz(1,i), xyz(2,i) );
      GC_namespace::to_mxArray(xyz,arg_out_0);
    } else if ( nlhs == 3 ) {
      GC_namespace::vec_real_type x(sz), y(sz), z(sz);
      for ( mwSize i = 0 ; i < sz ; ++i ) ptr->position( t[i], x[i], y[i], z[i] );
      GC_namespace::to_mxArray(x,arg_out_0);
      GC_namespace::to_mxArray(y,arg_out_1);
      GC_namespace::to_mxArray(z,arg_out_2);
    } else {
      UTILS_ASSERT( false, CMD ": expected 1 or 3 output, nlhs = {}\n", nlhs );
    }
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_velocity(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_7 "[vx,vy,vz] = AstroMexWrapper('velocity',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_7

    CHECK_IN(3);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    if ( nlhs == 1 ) {
      GC_namespace::mat_real_type vxyz;
      vxyz.resize(3,sz);
      for ( mwSize i = 0 ; i < sz ; ++i ) ptr->velocity( t[i], vxyz(0,i), vxyz(1,i), vxyz(2,i) );
      GC_namespace::to_mxArray(vxyz,arg_out_0);
    } else if ( nlhs == 3 ) {
      GC_namespace::vec_real_type vx(sz), vy(sz), vz(sz);
      for ( mwSize i = 0 ; i < sz ; ++i ) ptr->velocity( t[i], vx[i], vy[i], vz[i] );
      GC_namespace::to_mxArray(vx,arg_out_0);
      GC_namespace::to_mxArray(vy,arg_out_1);
      GC_namespace::to_mxArray(vz,arg_out_2);
    } else {
      UTILS_ASSERT( false, CMD ": expected 1 or 3 output, nlhs = {}\n", nlhs );
    }
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_acceleration(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_8 "[ax,ay,az] = AstroMexWrapper('acceleration',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_8

    CHECK_IN(3);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    if ( nlhs == 1 ) {
      GC_namespace::mat_real_type xyz;
      xyz.resize(3,sz);
      for ( mwSize i = 0 ; i < sz ; ++i ) ptr->acceleration( t[i], xyz(0,i), xyz(1,i), xyz(2,i) );
      GC_namespace::to_mxArray(xyz,arg_out_0);
    } else if ( nlhs == 3 ) {
      GC_namespace::vec_real_type ax(sz), ay(sz), az(sz);
      for ( mwSize i = 0 ; i < sz ; ++i ) ptr->acceleration( t[i], ax[i], ay[i], az[i] );
      GC_namespace::to_mxArray(ax,arg_out_0);
      GC_namespace::to_mxArray(ay,arg_out_1);
      GC_namespace::to_mxArray(az,arg_out_2);
    } else {
      UTILS_ASSERT( false, CMD ": expected 1 or 3 output, nlhs = {}\n", nlhs );
    }
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_jerk(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_9 "[jx,jy,jz] = AstroMexWrapper('jerk',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_9

    CHECK_IN(3);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ":  param t" );
    if ( nlhs == 1 ) {
      GC_namespace::mat_real_type xyz;
      xyz.resize(3,sz);
      for ( mwSize i = 0 ; i < sz ; ++i ) ptr->jerk( t[i], xyz(0,i), xyz(1,i), xyz(2,i) );
      GC_namespace::to_mxArray(xyz,arg_out_0);
    } else if ( nlhs == 3 ) {
      GC_namespace::vec_real_type jx(sz), jy(sz), jz(sz);
      for ( mwSize i = 0 ; i < sz ; ++i ) ptr->jerk( t[i], jx[i], jy[i], jz[i] );
      GC_namespace::to_mxArray(jx,arg_out_0);
      GC_namespace::to_mxArray(jy,arg_out_1);
      GC_namespace::to_mxArray(jz,arg_out_2);
    } else {
      UTILS_ASSERT( false, CMD ": expected 1 or 3 output, nlhs = {}\n", nlhs );
    }
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_setup_Keplerian(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_10 "AstroMexWrapper('setup_Keplerian',obj,[name,t0,a,e,Omega,omega,i,M0,muS]||struct)"
    #define CMD MEX_ERROR_MESSAGE_10

    UTILS_ASSERT( nrhs == 3 || nrhs == 11, CMD ": expected 3 or 11 input, nrhs = {}\n", nrhs );
    CHECK_OUT(0);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    string name = "no-name";
    GC_namespace::real_type t0(0), a(0), e(0), Omega(0), omega(0), i(0), M0(0), muS(0);
    if ( nrhs == 11 ) {
      UTILS_ASSERT0( mxIsChar(arg_in_2), CMD ": param 'name' must be a string" );
      name  = mxArrayToString(arg_in_2);
      t0    = Utils::mex_get_scalar_value( arg_in_3,  CMD ": param t0" );
      a     = Utils::mex_get_scalar_value( arg_in_4,  CMD ": param a" );
      e     = Utils::mex_get_scalar_value( arg_in_5,  CMD ": param e" );
      Omega = Utils::mex_get_scalar_value( arg_in_6,  CMD ": param Omega" );
      omega = Utils::mex_get_scalar_value( arg_in_7,  CMD ": param omega" );
      i     = Utils::mex_get_scalar_value( arg_in_8,  CMD ": param i" );
      M0    = Utils::mex_get_scalar_value( arg_in_9,  CMD ": param M0" );
      muS   = Utils::mex_get_scalar_value( arg_in_10, CMD ": param muS" );
    } else {
      GC_namespace::GenericContainer gc;
      try {
        GC_namespace::mxArray_to_GenericContainer(arg_in_2,gc);
        name = gc("name").get_string(CMD ": param struct field 'name'" );
      } catch ( std::exception const & e ) {
        mexErrMsgTxt( fmt::format( "Astro Error: {}", e.what() ).c_str() );
      } catch (...) {
        mexErrMsgTxt( "mxArray_to_GenericContainer failed" );
      }
      t0    = gc("t0").get_number();
      a     = gc("a").get_number();
      e     = gc("e").get_number();
      Omega = gc("Omega").get_number();
      omega = gc("omega").get_number();
      i     = gc("i").get_number();
      M0    = gc("M0").get_number();
      muS   = gc("muS").get_number();
    }
    ptr->setup_Keplerian( name, t0, a, e, Omega, omega, i, M0, muS );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_setup_Equinoctial(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_11 "AstroMexWrapper('setup_Equinoctial',obj,[name,t0,p,f,g,h,k,retrograde,L0,muS]||struct)"
    #define CMD MEX_ERROR_MESSAGE_11

    UTILS_ASSERT( nrhs == 3 || nrhs == 12, CMD ": expected 3 or 12 input, nrhs = {}\n", nrhs );
    CHECK_OUT(0);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type t0(0), p(0), f(0), g(0), h(0), k(0), L0(0), muS(0);
    bool retrograde = false;
    string name = "";
    if ( nrhs == 12 ) {
      UTILS_ASSERT0( mxIsChar(arg_in_2), CMD ": param 'name' must be a string" );
      name       = mxArrayToString(arg_in_2);
      t0         = Utils::mex_get_scalar_value( arg_in_3,  CMD ": param t0" );
      p          = Utils::mex_get_scalar_value( arg_in_4,  CMD ": param p" );
      f          = Utils::mex_get_scalar_value( arg_in_5,  CMD ": param f" );
      g          = Utils::mex_get_scalar_value( arg_in_6,  CMD ": param g" );
      h          = Utils::mex_get_scalar_value( arg_in_7,  CMD ": param h" );
      k          = Utils::mex_get_scalar_value( arg_in_8,  CMD ": param k" );
      retrograde = Utils::mex_get_bool( arg_in_9, CMD ": param 'retrograde'" );
      L0         = Utils::mex_get_scalar_value( arg_in_10, CMD ": param L0" );
      muS        = Utils::mex_get_scalar_value( arg_in_11, CMD ": param muS" );
    } else {
      GC_namespace::GenericContainer gc;
      try {
        GC_namespace::mxArray_to_GenericContainer(arg_in_2,gc);
        name = gc("name").get_string(CMD ": param struct field 'name'" );
      } catch ( std::exception const & e ) {
        mexErrMsgTxt( fmt::format( "Astro Error: {}", e.what() ).c_str() );
      } catch (...) {
        mexErrMsgTxt( "mxArray_to_GenericContainer failed" );
      }
      GenericContainer const & R = gc("retrograde");
      retrograde = false;
      if ( R.get_type() == GC_namespace::GC_BOOL ) {
        retrograde = R.get_bool();
      } else {
        retrograde = R.get_as_int() != 0;
      }
      name = gc("name").get_string(CMD ": param struct field 'name'" );
      t0   = gc("t0").get_number();
      p    = gc("p").get_number();
      f    = gc("f").get_number();
      g    = gc("g").get_number();
      h    = gc("h").get_number();
      k    = gc("k").get_number();
      L0   = gc("L0").get_number();
      muS  = gc("muS").get_number();
    }
    ptr->setup_Equinoctial( name, t0, p, f, g, h, k, retrograde, L0, muS );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_setup_PV(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_12 "AstroMexWrapper('setup_PV',obj,name,P,V,t0,muS)"
    #define CMD MEX_ERROR_MESSAGE_12

    CHECK_IN(7);
    CHECK_OUT(0);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    UTILS_ASSERT0( mxIsChar(arg_in_2), CMD ": param name must be a string" );
    string n = mxArrayToString(arg_in_2);
    mwSize szP, szV;
    GC_namespace::real_type const * P = Utils::mex_vector_pointer( arg_in_3, szP, CMD ": param P" );
    UTILS_ASSERT( szP == 3, CMD ": size |P| = {} expected 3\n", szP );
    GC_namespace::real_type const * V = Utils::mex_vector_pointer( arg_in_4, szV, CMD ": param V" );
    UTILS_ASSERT( szV == 3, CMD ": size |V| = {} expected 3\n", szV );
    GC_namespace::real_type t0  = Utils::mex_get_scalar_value( arg_in_5, CMD ": param t0" );
    GC_namespace::real_type muS = Utils::mex_get_scalar_value( arg_in_6, CMD ": param muS" );
    ptr->setup_using_point_and_velocity( n, P, V, muS, t0 );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_mean_anomaly(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_13 "res = AstroMexWrapper('mean_anomaly',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_13

    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    GC_namespace::vec_real_type ma;
    ma.reserve(sz);
    for ( mwSize i = 0 ; i < sz ; ++i ) ma.push_back( ptr->mean_anomaly( t[i] ) );
    GC_namespace::to_mxArray(ma,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_true_anomaly(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_14 "res = AstroMexWrapper('true_anomaly',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_14

    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    GC_namespace::vec_real_type ta;
    ta.reserve(sz);
    for ( mwSize i = 0 ; i < sz ; ++i ) ta.push_back( ptr->true_anomaly( t[i] ) );
    GC_namespace::to_mxArray(ta,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_p_orbital(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_15 "res = AstroMexWrapper('p_orbital',obj)"
    #define CMD MEX_ERROR_MESSAGE_15

    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->p_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_f_orbital(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_16 "res = AstroMexWrapper('f_orbital',obj)"
    #define CMD MEX_ERROR_MESSAGE_16

    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->f_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_g_orbital(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_17 "res = AstroMexWrapper('g_orbital',obj)"
    #define CMD MEX_ERROR_MESSAGE_17

    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->g_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_h_orbital(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_18 "res = AstroMexWrapper('h_orbital',obj)"
    #define CMD MEX_ERROR_MESSAGE_18

    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->h_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_k_orbital(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_19 "AstroMexWrapper('k_orbital',obj)"
    #define CMD MEX_ERROR_MESSAGE_19

    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->k_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_L_orbital(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_20 "L = AstroMexWrapper('L_orbital',obj,[t0,dt] or [t])"
    #define CMD MEX_ERROR_MESSAGE_20

    UTILS_ASSERT(
      nrhs == 3 || nrhs == 4, CMD "Expected 3 or 4 argument(s), nrhs = {}\n", nrhs
    );
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    if ( nrhs == 3 ) {
      mwSize sz;
      double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
      GC_namespace::vec_real_type L;
      L.reserve(sz);
      for ( mwSize i = 0 ; i < sz ; ++i ) L.push_back( ptr->L_orbital( t[i] ) );
      GC_namespace::to_mxArray(L,arg_out_0);
    } else {
      GC_namespace::real_type t0 = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t0" );
      mwSize sz;
      double const * dt = Utils::mex_vector_pointer( arg_in_3, sz, CMD ": param dt" );
      GC_namespace::vec_real_type L;
      L.reserve(sz);
      for ( mwSize i = 0 ; i < sz ; ++i ) L.push_back( ptr->L_orbital( t0, dt[i] ) );
      GC_namespace::to_mxArray(L,arg_out_0);
    }
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_L_orbital_D(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_21 "DLDt = AstroMexWrapper('L_orbital_D',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_21

    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    GC_namespace::vec_real_type L;
    L.reserve(sz);
    for ( mwSize i = 0 ; i < sz ; ++i ) L.push_back( ptr->L_orbital_D( t[i] ) );
    GC_namespace::to_mxArray(L,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_L_orbital_DD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_22 "D2LD2t = AstroMexWrapper('L_orbital_DD',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_22

    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    GC_namespace::vec_real_type L;
    L.reserve(sz);
    for ( mwSize i = 0 ; i < sz ; ++i ) L.push_back( ptr->L_orbital_DD( t[i] ) );
    GC_namespace::to_mxArray(L,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_latitude_of_periapsis(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_23 "res = AstroMexWrapper('latitude_of_periapsis',obj)"
    #define CMD MEX_ERROR_MESSAGE_23

    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::to_mxArray(ptr->latitude_of_periapsis(),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_latitude_of_apoapsis(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_24 "res = AstroMexWrapper('latitude_of_apoapsis',obj)"
    #define CMD MEX_ERROR_MESSAGE_24

    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::to_mxArray(ptr->latitude_of_apoapsis(),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_info(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_25 "info = AstroMexWrapper('info',obj)"
    #define CMD MEX_ERROR_MESSAGE_25

    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    
    GC_namespace::GenericContainer gc;
    GC_namespace::GenericContainer & EQ = gc["Equinoctial"];
    EQ["p"] = ptr->p_orbital();
    EQ["f"] = ptr->f_orbital();
    EQ["g"] = ptr->g_orbital();
    EQ["h"] = ptr->h_orbital();
    EQ["k"] = ptr->k_orbital();
    EQ["retrograde"] = ptr->retrograde();
    GC_namespace::GenericContainer & K = gc["Kepler"];
    K["a"]     = ptr->a_orbital();
    K["e"]     = ptr->e_orbital();
    K["i"]     = ptr->i_orbital();
    K["Omega"] = ptr->Omega_orbital();
    K["omega"] = ptr->omega_orbital();

    gc["t0"] = ptr->t0_orbital();
    gc["M0"] = ptr->M0_orbital();
    
    GC_namespace::GenericContainer_to_mxArray(gc,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_print_info(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_26 "AstroMexWrapper('print_info',obj)"
    #define CMD MEX_ERROR_MESSAGE_26

    CHECK_IN(2);
    CHECK_OUT(0);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mexPrintf("%s\n",ptr->info().c_str());
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_number_of_revolution(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_27 "n = AstroMexWrapper('number_of_revolution',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_27

    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type t = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t" );
    GC_namespace::to_mxArray(ptr->number_of_revolution(t),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_time_from_L_angle(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_28 "time = AstroMexWrapper('time_from_L_angle',obj,t,L)"
    #define CMD MEX_ERROR_MESSAGE_28

    CHECK_IN(4);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type t = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t" );
    GC_namespace::real_type L = Utils::mex_get_scalar_value( arg_in_3, CMD ": param L" );
    GC_namespace::to_mxArray(ptr->time_from_L_angle(t,L),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_absolute_position(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_29 "p = AstroMexWrapper('absolute_position',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_29

    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type t = Utils::mex_get_scalar_value( arg_in_2, CMD " param t" );
    GC_namespace::to_mxArray(ptr->absolute_position(t),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_absolute_velocity(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_30 "v = AstroMexWrapper('absolute_velocity',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_30

    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type t = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t" );
    GC_namespace::to_mxArray(ptr->absolute_velocity(t),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_period(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_31 "period = AstroMexWrapper('period',obj)"
    #define CMD MEX_ERROR_MESSAGE_31
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::to_mxArray(ptr->period(),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_apoapsis(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_32 "a = AstroMexWrapper('apoapsis',obj)"
    #define CMD MEX_ERROR_MESSAGE_32
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::to_mxArray(ptr->apoapsis(),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_periapsis(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_33 "p = AstroMexWrapper('periapsis',obj)"
    #define CMD MEX_ERROR_MESSAGE_33
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::to_mxArray(ptr->periapsis(),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_muS(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_34 "muS = AstroMexWrapper('muS',obj)"
    #define CMD MEX_ERROR_MESSAGE_34
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::to_mxArray(ptr->get_muS(),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_orbit_energy(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_35 "E = AstroMexWrapper('orbit_energy',obj)"
    #define CMD MEX_ERROR_MESSAGE_35
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::to_mxArray(ptr->orbit_energy(),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_radius_by_L(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_36 "R = AstroMexWrapper('radius_by_L',obj,L)"
    #define CMD MEX_ERROR_MESSAGE_36
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * L = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param L" );
    GC_namespace::vec_real_type ray;
    ray.reserve(sz);
    for ( mwSize i = 0 ; i < sz ; ++i ) ray.push_back( ptr->radius_by_L( L[i] ) );
    GC_namespace::to_mxArray(ray,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_radius_by_L_D(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_37 "DRDL = AstroMexWrapper('radius_by_L_D',obj,L)"
    #define CMD MEX_ERROR_MESSAGE_37
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * L = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param L" );
    GC_namespace::vec_real_type ray;
    ray.reserve(sz);
    for ( mwSize i = 0 ; i < sz ; ++i ) ray.push_back( ptr->radius_by_L_D( L[i] ) );
    GC_namespace::to_mxArray(ray,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_radius_by_L_DD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_38 "D2RD2L = AstroMexWrapper('radius_by_L_DD',obj,L)"
    #define CMD MEX_ERROR_MESSAGE_38
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * L = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param L" );
    GC_namespace::vec_real_type ray;
    ray.reserve(sz);
    for ( mwSize i = 0 ; i < sz ; ++i ) ray.push_back( ptr->radius_by_L_DD( L[i] ) );
    GC_namespace::to_mxArray(ray,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_L_orbital_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_39 "L_grad = AstroMexWrapper('L_orbital_EQ_gradient',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_39
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    GC_namespace::mat_real_type G;
    G.resize(sz,6);
    for ( mwSize i = 0 ; i < sz ; ++i ) {
      real_type grad[6];
      ptr->L_orbital_EQ_gradient( t[i], grad );
      G(i,0) = grad[0]; G(i,1) = grad[1]; G(i,2) = grad[2];
      G(i,3) = grad[3]; G(i,4) = grad[4]; G(i,5) = grad[5];
    }
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_M0(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_40 "M0 = AstroMexWrapper('M0',obj)"
    #define CMD MEX_ERROR_MESSAGE_40
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->M0_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_M0_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_41 "M0_grad = AstroMexWrapper('M0_EQ_gradient',obj)"
    #define CMD MEX_ERROR_MESSAGE_41
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::vec_real_type G;
    G.resize(6);
    real_type grad[6];
    ptr->M0_EQ_gradient( grad );
    G[0] = grad[0]; G[1] = grad[1]; G[2] = grad[2];
    G[3] = grad[3]; G[4] = grad[4]; G[5] = grad[5];
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_theta0(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_42 "theta0 = AstroMexWrapper('theta0',obj)"
    #define CMD MEX_ERROR_MESSAGE_42
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->theta0_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_theta0_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_43 "theta0_grad = AstroMexWrapper('theta0_EQ_gradient',obj)"
    #define CMD MEX_ERROR_MESSAGE_43
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::vec_real_type G;
    G.resize(6);
    real_type grad[6];
    ptr->theta0_EQ_gradient( grad );
    G[0] = grad[0]; G[1] = grad[1]; G[2] = grad[2];
    G[3] = grad[3]; G[4] = grad[4]; G[5] = grad[5];
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_E0_angle(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_44 "E0 = AstroMexWrapper('E0_angle',obj)"
    #define CMD MEX_ERROR_MESSAGE_44
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->E0_angle() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_E0_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_45 "E0_grad = AstroMexWrapper('E0_EQ_gradient',obj)"
    #define CMD MEX_ERROR_MESSAGE_45
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::vec_real_type G;
    G.resize(6);
    real_type grad[6];
    ptr->E0_EQ_gradient( grad );
    G[0] = grad[0]; G[1] = grad[1]; G[2] = grad[2];
    G[3] = grad[3]; G[4] = grad[4]; G[5] = grad[5];
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_H0_angle(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_46 "H0 = AstroMexWrapper('H0_angle',obj)"
    #define CMD MEX_ERROR_MESSAGE_46
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->H0_angle() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_H0_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_47 "H0_grad = AstroMexWrapper('H0_EQ_gradient',obj)"
    #define CMD MEX_ERROR_MESSAGE_47
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::vec_real_type G;
    G.resize(6);
    real_type grad[6];
    ptr->H0_EQ_gradient( grad );
    G[0] = grad[0]; G[1] = grad[1]; G[2] = grad[2];
    G[3] = grad[3]; G[4] = grad[4]; G[5] = grad[5];
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_E_angle(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_48 "E = AstroMexWrapper('E_angle',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_48
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD " param t" );
    GC_namespace::vec_real_type E;
    E.reserve(sz);
    for ( mwSize i = 0 ; i < sz ; ++i ) E.push_back( ptr->E_angle( t[i] ) );
    GC_namespace::to_mxArray(E,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_E_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_49 "E_grad = AstroMexWrapper('E_EQ_gradient',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_49
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    GC_namespace::mat_real_type G;
    G.resize(sz,6);
    for ( mwSize i = 0 ; i < sz ; ++i ) {
      real_type grad[6];
      ptr->E_EQ_gradient( t[i], grad );
      G(i,0) = grad[0]; G(i,1) = grad[1]; G(i,2) = grad[2];
      G(i,3) = grad[3]; G(i,4) = grad[4]; G(i,5) = grad[5];
    }
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_H_angle(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_50 "H = AstroMexWrapper('H_angle',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_50
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    GC_namespace::vec_real_type H;
    H.reserve(sz);
    for ( mwSize i = 0 ; i < sz ; ++i ) H.push_back( ptr->H_angle( t[i] ) );
    GC_namespace::to_mxArray(H,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_H_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_51 "H_grad = AstroMexWrapper('H_EQ_gradient',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_51
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    GC_namespace::mat_real_type G;
    G.resize(sz,6);
    for ( mwSize i = 0 ; i < sz ; ++i ) {
      real_type grad[6];
      ptr->H_EQ_gradient( t[i], grad );
      G(i,0) = grad[0]; G(i,1) = grad[1]; G(i,2) = grad[2];
      G(i,3) = grad[3]; G(i,4) = grad[4]; G(i,5) = grad[5];
    }
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_radius_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_52 "R_grad = AstroMexWrapper('radius_EQ_gradient',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_52
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    GC_namespace::mat_real_type G;
    G.resize(sz,6);
    for ( mwSize i = 0 ; i < sz ; ++i ) {
      real_type grad[6];
      ptr->radius_EQ_gradient( t[i], grad );
      G(i,0) = grad[0]; G(i,1) = grad[1]; G(i,2) = grad[2];
      G(i,3) = grad[3]; G(i,4) = grad[4]; G(i,5) = grad[5];
    }
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_position_EQ_jacobian(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_53 "P_grad = AstroMexWrapper('position_EQ_jacobian',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_53
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    double t = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t" );
    GC_namespace::mat_real_type G;
    G.resize(3,6);
    real_type JP[3][6];
    ptr->position_EQ_jacobian( t, JP );
    for ( integer i = 0; i < 3; ++i )
      for ( integer j = 0; j < 6; ++j )
        G(i,j) = JP[i][j];
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_position0_EQ_jacobian(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_54 "P0_grad = AstroMexWrapper('position0_EQ_jacobian',obj)"
    #define CMD MEX_ERROR_MESSAGE_54
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::mat_real_type G;
    G.resize(3,6);
    real_type JP[3][6];
    ptr->position0_EQ_jacobian( JP, ptr->L0_orbital() );
    for ( integer i = 0; i < 3; ++i )
      for ( integer j = 0; j < 6; ++j )
        G(i,j) = JP[i][j];
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_position_EQ_jacobian_FD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_55 "P_grad = AstroMexWrapper('position_EQ_jacobian_FD',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_55
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    double t = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t" );
    GC_namespace::mat_real_type G;
    G.resize(3,6);
    real_type JP[3][6];
    ptr->position_EQ_jacobian_FD( t, JP );
    for ( integer i = 0; i < 3; ++i )
      for ( integer j = 0; j < 6; ++j )
        G(i,j) = JP[i][j];
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_absolute_velocity_by_angle(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_56 "V = AstroMexWrapper('absolute_velocity_by_angle',obj,L)"
    #define CMD MEX_ERROR_MESSAGE_56
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type L = Utils::mex_get_scalar_value( arg_in_2, CMD ": param L" );
    GC_namespace::to_mxArray(ptr->absolute_velocity_by_angle(L),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_absolute_velocity_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_57 "V_grad = AstroMexWrapper('absolute_velocity_EQ_gradient',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_57
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * t = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param t" );
    GC_namespace::mat_real_type G;
    G.resize(sz,6);
    for ( mwSize i = 0 ; i < sz ; ++i ) {
      real_type grad[6];
      ptr->absolute_velocity_EQ_gradient( t[i], grad );
      G(i,0) = grad[0]; G(i,1) = grad[1]; G(i,2) = grad[2];
      G(i,3) = grad[3]; G(i,4) = grad[4]; G(i,5) = grad[5];
    }
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_velocity0_EQ_jacobian(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_58 "V0_grad = AstroMexWrapper('velocity0_EQ_jacobian',obj)"
    #define CMD MEX_ERROR_MESSAGE_58
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::mat_real_type G;
    G.resize(3,6);
    real_type JP[3][6];
    ptr->velocity0_EQ_jacobian( JP, ptr->L0_orbital() );
    for ( integer i = 0; i < 3; ++i )
      for ( integer j = 0; j < 6; ++j )
        G(i,j) = JP[i][j];
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_velocity_EQ_jacobian(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_59 "V_jacobian = AstroMexWrapper('velocity_EQ_jacobian',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_59
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    double t = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t" );
    GC_namespace::mat_real_type G;
    G.resize(3,6);
    real_type JP[3][6];
    ptr->velocity_EQ_jacobian( t, JP );
    for ( integer i = 0; i < 3; ++i )
      for ( integer j = 0; j < 6; ++j )
        G(i,j) = JP[i][j];
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_velocity_EQ_jacobian_FD(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_60 "V_jacobian_FD = AstroMexWrapper('velocity_EQ_jacobian_FD',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_60
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    double t = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t" );
    GC_namespace::mat_real_type G;
    G.resize(3,6);
    real_type JP[3][6];
    ptr->velocity_EQ_jacobian_FD( t, JP );
    for ( integer i = 0; i < 3; ++i )
      for ( integer j = 0; j < 6; ++j )
        G(i,j) = JP[i][j];
    GC_namespace::to_mxArray(G,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_L_from_true_anomaly(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_61 "L = AstroMexWrapper('L_from_true_anomaly',obj,nu)"
    #define CMD MEX_ERROR_MESSAGE_61
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mwSize sz;
    double const * nu = Utils::mex_vector_pointer( arg_in_2, sz, CMD ": param nu" );
    GC_namespace::vec_real_type L;
    L.reserve(sz);
    for ( mwSize i = 0 ; i < sz ; ++i ) L.push_back( ptr->L_from_true_anomaly( nu[i] ) );
    GC_namespace::to_mxArray(L,arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_make_retrograde(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_62 "AstroMexWrapper('make_retrograde',obj)"
    #define CMD MEX_ERROR_MESSAGE_62
    CHECK_IN(2);
    CHECK_OUT(0);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    ptr->make_retrograde();
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_make_not_retrograde(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_63 "AstroMexWrapper('make_not_retrograde',obj)"
    #define CMD MEX_ERROR_MESSAGE_63
    CHECK_IN(2);
    CHECK_OUT(0);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    ptr->make_not_retrograde();
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_retrograde(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_64 "res = AstroMexWrapper('retrograde',obj)"
    #define CMD MEX_ERROR_MESSAGE_64
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    Utils::mex_set_scalar_value( arg_out_0, ptr->retrograde() ? 1 : 0 );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_normal(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_65 "N = AstroMexWrapper('normal',obj)"
    #define CMD MEX_ERROR_MESSAGE_65
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type * n = Utils::mex_create_matrix_value( arg_out_0, 3, 1 );
    ptr->normal(n);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_local_frame(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_66 "frame = AstroMexWrapper('local_frame',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_66
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type t = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t" );
    GC_namespace::real_type * frame = Utils::mex_create_matrix_value( arg_out_0, 3, 3 );
    GC_namespace::real_type M[3][3];
    ptr->local_frame(t,M);
    for ( mwSize i=0; i < 3; ++i )
      for ( mwSize j=0; j < 3; ++j )
        frame[i+j*3] = M[i][j];
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_local_frame_by_L(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_67 "frame = AstroMexWrapper('local_frame_by_L',obj,L)"
    #define CMD MEX_ERROR_MESSAGE_67
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type L = Utils::mex_get_scalar_value( arg_in_2, CMD ": param L" );
    GC_namespace::real_type * frame = Utils::mex_create_matrix_value( arg_out_0, 3, 3 );
    GC_namespace::real_type M[3][3];
    ptr->local_frame_by_L(L,M);
    for ( mwSize i=0; i < 3; ++i )
      for ( mwSize j=0; j < 3; ++j )
        frame[i+j*3] = M[i][j];
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_ellipse_frame(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_68 "frame = AstroMexWrapper('ellipse_frame',obj)"
    #define CMD MEX_ERROR_MESSAGE_68
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type * frame = Utils::mex_create_matrix_value( arg_out_0, 3, 3 );
    GC_namespace::real_type M[3][3];
    ptr->ellipse_frame(M);
    for ( mwSize i=0; i < 3; ++i )
      for ( mwSize j=0; j < 3; ++j )
        frame[i+j*3] = M[i][j];
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_eval_E(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_69 "E = AstroMexWrapper('eval_E',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_69
    CHECK_IN(3);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type t = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t" );
    UTILS_ASSERT(
      nrhs > 0 && nrhs < 4,
      CMD "Expected 1,2,3 or 4  argument(s), nrhs = {}\n", nrhs
    );

    GC_namespace::real_type E[4];
    ptr->eval_E(t,E,nrhs-1);
    Utils::mex_set_scalar_value( arg_out_0, E[0] );
    if ( nlhs > 1 ) Utils::mex_set_scalar_value( arg_out_1, E[1] );
    if ( nlhs > 2 ) Utils::mex_set_scalar_value( arg_out_2, E[2] );
    if ( nlhs > 3 ) Utils::mex_set_scalar_value( arg_out_3, E[3] );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_eval_L(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_70 "L = AstroMexWrapper('eval_L',obj,t)"
    #define CMD MEX_ERROR_MESSAGE_70
    CHECK_IN(3);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type t = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t" );
    UTILS_ASSERT(
      nrhs > 0 && nrhs < 4,
      CMD "Expected 1,2,3 or 4  argument(s), nrhs = {}\n", nrhs
    );

    GC_namespace::real_type L[4];
    ptr->eval_L(t,L,nrhs-1);
    Utils::mex_set_scalar_value( arg_out_0, L[0] );
    if ( nlhs > 1 ) Utils::mex_set_scalar_value( arg_out_1, L[1] );
    if ( nlhs > 2 ) Utils::mex_set_scalar_value( arg_out_2, L[2] );
    if ( nlhs > 3 ) Utils::mex_set_scalar_value( arg_out_3, L[3] );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static std::map<std::string,DO_CMD> cmd_to_fun = {
    {"new",do_new},
    {"delete",do_delete},
    {"copy",do_copy},
    {"name",do_name},
    {"check_for_consistency",do_check_for_consistency},
    {"position",do_position},
    {"velocity",do_velocity},
    {"acceleration",do_acceleration},
    {"jerk",do_jerk},
    {"setup_Keplerian",do_setup_Keplerian},
    {"setup_Equinoctial",do_setup_Equinoctial},
    {"setup_PV",do_setup_PV},
    {"mean_anomaly",do_mean_anomaly},
    {"true_anomaly",do_true_anomaly},
    {"p_orbital",do_p_orbital},
    {"f_orbital",do_f_orbital},
    {"g_orbital",do_g_orbital},
    {"h_orbital",do_h_orbital},
    {"k_orbital",do_k_orbital},
    {"L_orbital",do_L_orbital},
    {"L_orbital_D",do_L_orbital_D},
    {"L_orbital_DD",do_L_orbital_DD},
    {"latitude_of_periapsis",do_latitude_of_periapsis},
    {"latitude_of_apoapsis",do_latitude_of_apoapsis},
    {"info",do_info},
    {"print_info",do_print_info},
    {"number_of_revolution",do_number_of_revolution},
    {"time_from_L_angle",do_time_from_L_angle},
    {"absolute_position",do_absolute_position},
    {"absolute_velocity",do_absolute_velocity},
    {"muS",do_muS},
    {"period",do_period},
    {"apoapsis",do_apoapsis},
    {"periapsis",do_periapsis},
    {"orbit_energy",do_orbit_energy},
    {"radius_by_L",do_radius_by_L},
    {"radius_by_L_D",do_radius_by_L_D},
    {"radius_by_L_DD",do_radius_by_L_DD},
    {"L_orbital_EQ_gradient",do_L_orbital_EQ_gradient},
    {"M0",do_M0},
    {"M0_EQ_gradient",do_M0_EQ_gradient},
    {"theta0",do_theta0},
    {"theta0_EQ_gradient",do_theta0_EQ_gradient},
    {"E0_angle",do_E0_angle},
    {"E0_EQ_gradient",do_E0_EQ_gradient},
    {"E_angle",do_E_angle},
    {"E_EQ_gradient",do_E_EQ_gradient},
    {"H0_angle",do_H0_angle},
    {"H0_EQ_gradient",do_H0_EQ_gradient},
    {"H_angle",do_H_angle},
    {"H_EQ_gradient",do_H_EQ_gradient},
    {"radius_EQ_gradient",do_radius_EQ_gradient},
    {"position_EQ_jacobian",do_position_EQ_jacobian},
    {"position_EQ_jacobian_FD",do_position_EQ_jacobian_FD},
    {"position0_EQ_jacobian",do_position0_EQ_jacobian},
    {"velocity_EQ_jacobian",do_velocity_EQ_jacobian},
    {"velocity_EQ_jacobian_FD",do_velocity_EQ_jacobian_FD},
    {"velocity0_EQ_jacobian",do_velocity0_EQ_jacobian},
    {"absolute_velocity_by_angle",do_absolute_velocity_by_angle},
    {"absolute_velocity_EQ_gradient",do_absolute_velocity_EQ_gradient},
    {"L_from_true_anomaly",do_L_from_true_anomaly},
    {"make_retrograde",do_make_retrograde},
    {"make_not_retrograde",do_make_not_retrograde},
    {"retrograde",do_retrograde},
    {"normal",do_normal},
    {"local_frame",do_local_frame},
    {"local_frame_by_L",do_local_frame_by_L},
    {"ellipse_frame",do_ellipse_frame},
    {"eval_E",do_eval_E},
    {"eval_L",do_eval_L}
  };

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  #define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"AstroMexWrapper\n" \
"\n" \
"USAGE:\n" \
"\n" \
MEX_ERROR_MESSAGE_1 "\n" \
MEX_ERROR_MESSAGE_2 "\n" \
MEX_ERROR_MESSAGE_3 "\n" \
MEX_ERROR_MESSAGE_4 "\n" \
MEX_ERROR_MESSAGE_5 "\n" \
MEX_ERROR_MESSAGE_6 "\n" \
MEX_ERROR_MESSAGE_7 "\n" \
MEX_ERROR_MESSAGE_8 "\n" \
MEX_ERROR_MESSAGE_9 "\n" \
MEX_ERROR_MESSAGE_10 "\n" \
MEX_ERROR_MESSAGE_11 "\n" \
MEX_ERROR_MESSAGE_12 "\n" \
MEX_ERROR_MESSAGE_13 "\n" \
MEX_ERROR_MESSAGE_14 "\n" \
MEX_ERROR_MESSAGE_15 "\n" \
MEX_ERROR_MESSAGE_16 "\n" \
MEX_ERROR_MESSAGE_17 "\n" \
MEX_ERROR_MESSAGE_18 "\n" \
MEX_ERROR_MESSAGE_19 "\n" \
MEX_ERROR_MESSAGE_20 "\n" \
MEX_ERROR_MESSAGE_21 "\n" \
MEX_ERROR_MESSAGE_22 "\n" \
MEX_ERROR_MESSAGE_23 "\n" \
MEX_ERROR_MESSAGE_24 "\n" \
MEX_ERROR_MESSAGE_25 "\n" \
MEX_ERROR_MESSAGE_26 "\n" \
MEX_ERROR_MESSAGE_27 "\n" \
MEX_ERROR_MESSAGE_28 "\n" \
MEX_ERROR_MESSAGE_29 "\n" \
MEX_ERROR_MESSAGE_30 "\n" \
MEX_ERROR_MESSAGE_31 "\n" \
MEX_ERROR_MESSAGE_32 "\n" \
MEX_ERROR_MESSAGE_33 "\n" \
MEX_ERROR_MESSAGE_34 "\n" \
MEX_ERROR_MESSAGE_35 "\n" \
MEX_ERROR_MESSAGE_36 "\n" \
MEX_ERROR_MESSAGE_37 "\n" \
MEX_ERROR_MESSAGE_38 "\n" \
MEX_ERROR_MESSAGE_39 "\n" \
MEX_ERROR_MESSAGE_40 "\n" \
MEX_ERROR_MESSAGE_41 "\n" \
MEX_ERROR_MESSAGE_42 "\n" \
MEX_ERROR_MESSAGE_43 "\n" \
MEX_ERROR_MESSAGE_44 "\n" \
MEX_ERROR_MESSAGE_45 "\n" \
MEX_ERROR_MESSAGE_46 "\n" \
MEX_ERROR_MESSAGE_47 "\n" \
MEX_ERROR_MESSAGE_48 "\n" \
MEX_ERROR_MESSAGE_49 "\n" \
MEX_ERROR_MESSAGE_50 "\n" \
MEX_ERROR_MESSAGE_51 "\n" \
MEX_ERROR_MESSAGE_52 "\n" \
MEX_ERROR_MESSAGE_53 "\n" \
MEX_ERROR_MESSAGE_54 "\n" \
MEX_ERROR_MESSAGE_55 "\n" \
MEX_ERROR_MESSAGE_56 "\n" \
MEX_ERROR_MESSAGE_57 "\n" \
MEX_ERROR_MESSAGE_58 "\n" \
MEX_ERROR_MESSAGE_59 "\n" \
MEX_ERROR_MESSAGE_60 "\n" \
MEX_ERROR_MESSAGE_61 "\n" \
MEX_ERROR_MESSAGE_62 "\n" \
MEX_ERROR_MESSAGE_63 "\n" \
MEX_ERROR_MESSAGE_64 "\n" \
MEX_ERROR_MESSAGE_65 "\n" \
MEX_ERROR_MESSAGE_66 "\n" \
MEX_ERROR_MESSAGE_67 "\n" \
MEX_ERROR_MESSAGE_68 "\n" \
MEX_ERROR_MESSAGE_69 "\n" \
MEX_ERROR_MESSAGE_70 "\n" \
"\n" \
"%======================================================================%\n" \
"%                                                                      %\n" \
"%  Autor: Enrico Bertolazzi                                            %\n" \
"%         Department of Industrial Engineering                         %\n" \
"%         University of Trento                                         %\n" \
"%         enrico.bertolazzi@unitn.it                                   %\n" \
"%                                                                      %\n" \
"%======================================================================%\n"

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  extern "C"
  void
  mexFunction(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    // the first argument must be a string
    if ( nrhs == 0 ) {
      mexErrMsgTxt(MEX_ERROR_MESSAGE);
      return;
    }

    try {
      UTILS_ASSERT0( mxIsChar(arg_in_0), "AstroMexWrapper: first argument must be a string" );
      string cmd = mxArrayToString(arg_in_0);
      std::map<std::string,DO_CMD>::const_iterator it = cmd_to_fun.find(cmd);
      if ( it != cmd_to_fun.end() ) {
        it->second( nlhs, plhs, nrhs, prhs );
      } else {
        mexErrMsgTxt( fmt::format( "AstroMexWrapper Error: key {} unknown", cmd ).c_str() );
      }
    } catch ( std::exception const & e ) {
      mexErrMsgTxt( fmt::format( "AstroMexWrapper Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt( "AstroMexWrapper failed" );
    }
  }
}
