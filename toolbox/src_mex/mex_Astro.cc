/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2019
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "mex_utils.hh"
#include "Astro.hh"

#define MEX_ERROR_MESSAGE \
"=====================================================================================\n" \
"mex_Astro: error\n" \
"\n" \
"USAGE:\n" \
"  - Constructors:\n" \
"    OBJ = AstroMexWrapper( 'new' );\n" \
"\n" \
"  On output:\n" \
"    OBJ = pointer to the internal object\n" \
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

  /*\
   |  ____    _  _____  _
   | |  _ \  / \|_   _|/ \
   | | | | |/ _ \ | | / _ \
   | | |_| / ___ \| |/ ___ \
   | |____/_/   \_\_/_/   \_\
   |
  \*/

  static
  inline
  Astro *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<Astro>(mx_id);
  }

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
    #define CMD "AstroMexWrapper('new'[,'name']): "
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 1 || nrhs == 2, CMD "expected 1 or 2 input, nrhs = {}\n", nrhs );
    Astro * ptr = nullptr;
    if ( nrhs == 1 ) {
      ptr = new Astro();
    } else {
      MEX_ASSERT( mxIsChar(arg_in_1), CMD "second argument must be a string" );
      string name = mxArrayToString(arg_in_1);
      ptr = new Astro(name);
    }
    arg_out_0 = convertPtr2Mat<Astro>(ptr);
    #undef CMD
  }

  static
  void
  do_new_handle(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('new_handle',old_handle,new_handle): "
    MEX_ASSERT2( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 3, CMD "expected 2 input, nrhs = {}\n", nrhs );
    Astro * ptr = DATA_GET( arg_in_2 );
    Mat2Ptr_change_ptr<Astro>( arg_in_1, ptr );
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
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('delete',obj): "
    CHECK_IN(2);
    CHECK_OUT(0);
    destroyObject<Astro>( arg_in_1 );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_copy(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('copy',obj,obj1): "
    CHECK_IN(3);
    CHECK_OUT(0);
    Astro * ptr  = DATA_GET( arg_in_1 );
    Astro * ptr1 = DATA_GET( arg_in_2 );
    ptr1->setup(
      ptr->name(), ptr->t0_orbital(), ptr->get_K(),
      ptr->M0_orbital(), ptr->get_muS()
    );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_name(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('do_name',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('check_for_consistency',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('position',obj,t): "
    CHECK_IN(3);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
      MEX_ASSERT2( false, CMD "expected 1 or 3 output, nlhs = {}\n", nlhs );
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
    #define CMD "AstroMexWrapper('velocity',obj,t): "
    CHECK_IN(3);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
      MEX_ASSERT2( false, CMD "expected 1 or 3 output, nlhs = {}\n", nlhs );
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
    #define CMD "AstroMexWrapper('acceleration',obj,t): "
    CHECK_IN(3);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
      MEX_ASSERT2( false, CMD "expected 1 or 3 output, nlhs = {}\n", nlhs );
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
    #define CMD "AstroMexWrapper('jerk',obj,t): "
    CHECK_IN(3);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
      MEX_ASSERT2( false, CMD "expected 1 or 3 output, nlhs = {}\n", nlhs );
    }
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_setup_Keplerian(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('setup_Keplerian',obj,[name,t0,a,e,Omega,omega,i,M0,muS]||struct): "
    MEX_ASSERT2( nrhs == 3 || nrhs == 11, CMD "expected 3 or 11 input, nrhs = {}\n", nrhs );
    CHECK_OUT(0);
    Astro * ptr = DATA_GET( arg_in_1 );
    string name = "no-name";
    GC_namespace::real_type t0(0), a(0), e(0), Omega(0), omega(0), i(0), M0(0), muS(0);
    if ( nrhs == 11 ) {
      MEX_ASSERT( mxIsChar(arg_in_2), CMD " param 'name' must be a string" );
      name  = mxArrayToString(arg_in_2);
      t0    = getScalarValue( arg_in_3,  CMD " param t0" );
      a     = getScalarValue( arg_in_4,  CMD " param a" );
      e     = getScalarValue( arg_in_5,  CMD " param e" );
      Omega = getScalarValue( arg_in_6,  CMD " param Omega" );
      omega = getScalarValue( arg_in_7,  CMD " param omega" );
      i     = getScalarValue( arg_in_8,  CMD " param i" );
      M0    = getScalarValue( arg_in_9,  CMD " param M0" );
      muS   = getScalarValue( arg_in_10, CMD " param muS" );
    } else {
      GC_namespace::GenericContainer gc;
      try {
        GC_namespace::mxArray_to_GenericContainer(arg_in_2,gc);
        name = gc("name").get_string(CMD " param struct field 'name': " );
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
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('setup_Equinoctial',obj,[name,t0,p,f,g,h,k,retrograde,L0,muS]||struct): "
    MEX_ASSERT2( nrhs == 3 || nrhs == 12, CMD "expected 3 or 12 input, nrhs = {}\n", nrhs );
    CHECK_OUT(0);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type t0(0), p(0), f(0), g(0), h(0), k(0), L0(0), muS(0);
    bool retrograde = false;
    string name = "";
    if ( nrhs == 12 ) {
      MEX_ASSERT( mxIsChar(arg_in_2), CMD " param 'name' must be a string" );
      name       = mxArrayToString(arg_in_2);
      t0         = getScalarValue( arg_in_3,  CMD " param t0" );
      p          = getScalarValue( arg_in_4,  CMD " param p" );
      f          = getScalarValue( arg_in_5,  CMD " param f" );
      g          = getScalarValue( arg_in_6,  CMD " param g" );
      h          = getScalarValue( arg_in_7,  CMD " param h" );
      k          = getScalarValue( arg_in_8,  CMD " param k" );
      retrograde = getBool( arg_in_9, CMD " param 'retrograde'" );
      L0         = getScalarValue( arg_in_10, CMD " param L0" );
      muS        = getScalarValue( arg_in_11, CMD " param muS" );
    } else {
      GC_namespace::GenericContainer gc;
      try {
        GC_namespace::mxArray_to_GenericContainer(arg_in_2,gc);
        name = gc("name").get_string(CMD " param struct field 'name': " );
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
      name = gc("name").get_string(CMD " param struct field 'name'" );
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
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('setup_PV',obj,name,P,V,t0,muS): "
    CHECK_IN(7);
    CHECK_OUT(0);
    Astro * ptr = DATA_GET( arg_in_1 );
    MEX_ASSERT( mxIsChar(arg_in_2), CMD " param name must be a string" );
    string n = mxArrayToString(arg_in_2);
    mwSize szP, szV;
    GC_namespace::real_type const * P = getVectorPointer( arg_in_3, szP, CMD " param P" );
    MEX_ASSERT2( szP == 3, CMD "size |P| = {} expected 3\n", szP );
    GC_namespace::real_type const * V = getVectorPointer( arg_in_4, szV, CMD " param V" );
    MEX_ASSERT2( szV == 3, CMD "size |V| = {} expected 3\n", szV );
    GC_namespace::real_type t0  = getScalarValue( arg_in_5, CMD " param t0" );
    GC_namespace::real_type muS = getScalarValue( arg_in_6, CMD " param muS" );
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
    #define CMD "AstroMexWrapper('mean_anomaly',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
    #define CMD "AstroMexWrapper('true_anomaly',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
    #define CMD "AstroMexWrapper('p_orbital',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    setScalarValue( arg_out_0, ptr->p_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_f_orbital(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('f_orbital',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    setScalarValue( arg_out_0, ptr->f_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_g_orbital(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('g_orbital',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    setScalarValue( arg_out_0, ptr->g_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_h_orbital(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('h_orbital',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    setScalarValue( arg_out_0, ptr->h_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_k_orbital(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('k_orbital',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    setScalarValue( arg_out_0, ptr->k_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_L_orbital(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('L_orbital',obj,[t0,dt|t]): "
    MEX_ASSERT2(
      nrhs == 3 || nrhs == 4, CMD "Expected 3 or 4 argument(s), nrhs = {}\n", nrhs
    );
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    if ( nrhs == 3 ) {
      mwSize sz;
      double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
      GC_namespace::vec_real_type L;
      L.reserve(sz);
      for ( mwSize i = 0 ; i < sz ; ++i ) L.push_back( ptr->L_orbital( t[i] ) );
      GC_namespace::to_mxArray(L,arg_out_0);
    } else {
      GC_namespace::real_type t0 = getScalarValue( arg_in_2, CMD " param t0" );
      mwSize sz;
      double const * dt = getVectorPointer( arg_in_3, sz, CMD " param dt" );
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
    #define CMD "AstroMexWrapper('L_orbital_D',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
    #define CMD "AstroMexWrapper('L_orbital_DD',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
    #define CMD "AstroMexWrapper('latitude_of_periapsis',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('latitude_of_apoapsis',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('info',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    
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
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('print_info',obj): "
    CHECK_IN(2);
    CHECK_OUT(0);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('number_of_revolution',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type t = getScalarValue( arg_in_2, CMD " param t" );
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
    #define CMD "AstroMexWrapper('time_from_L_angle',obj,t,L): "
    CHECK_IN(4);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type t = getScalarValue( arg_in_2, CMD " param t" );
    GC_namespace::real_type L = getScalarValue( arg_in_3, CMD " param L" );
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
    #define CMD "AstroMexWrapper('absolute_position',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type t = getScalarValue( arg_in_2, CMD " param t" );
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
    #define CMD "AstroMexWrapper('absolute_velocity',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type t = getScalarValue( arg_in_2, CMD " param t" );
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
    #define CMD "AstroMexWrapper('period',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('apoapsis',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('periapsis',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('muS',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('orbit_energy',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('radius_by_L',obj,L): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * L = getVectorPointer( arg_in_2, sz, CMD " param L" );
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
    #define CMD "AstroMexWrapper('radius_by_L_D',obj,L): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * L = getVectorPointer( arg_in_2, sz, CMD " param L" );
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
    #define CMD "AstroMexWrapper('radius_by_L_DD',obj,L): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * L = getVectorPointer( arg_in_2, sz, CMD " param L" );
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
    #define CMD "AstroMexWrapper('L_orbital_EQ_gradient',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
    #define CMD "AstroMexWrapper('M0',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    setScalarValue( arg_out_0, ptr->M0_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_M0_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('M0_EQ_gradient',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('theta0',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    setScalarValue( arg_out_0, ptr->theta0_orbital() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_theta0_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('theta0_EQ_gradient',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('E0_angle',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    setScalarValue( arg_out_0, ptr->E0_angle() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_E0_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('E0_EQ_gradient',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('H0_angle',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    setScalarValue( arg_out_0, ptr->H0_angle() );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_H0_EQ_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('H0_EQ_gradient',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('E_angle',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
    #define CMD "AstroMexWrapper('E_EQ_gradient',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
    #define CMD "AstroMexWrapper('H_angle',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
    #define CMD "AstroMexWrapper('H_EQ_gradient',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
    #define CMD "AstroMexWrapper('radius_EQ_gradient',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
    #define CMD "AstroMexWrapper('position_EQ_jacobian',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    double t = getScalarValue( arg_in_2, CMD " param t" );
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
    #define CMD "AstroMexWrapper('position0_EQ_jacobian',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('position_EQ_jacobian_FD',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    double t = getScalarValue( arg_in_2, CMD " param t" );
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
    #define CMD "AstroMexWrapper('absolute_velocity_by_angle',obj,L): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type L = getScalarValue( arg_in_2, CMD " param L" );
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
    #define CMD "AstroMexWrapper('absolute_velocity_EQ_gradient',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * t = getVectorPointer( arg_in_2, sz, CMD " param t" );
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
    #define CMD "AstroMexWrapper('velocity0_EQ_jacobian',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('velocity_EQ_jacobian',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    double t = getScalarValue( arg_in_2, CMD " param t" );
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
    #define CMD "AstroMexWrapper('velocity_EQ_jacobian_FD',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    double t = getScalarValue( arg_in_2, CMD " param t" );
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
    #define CMD "AstroMexWrapper('L_from_true_anomaly',obj,nu): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    mwSize sz;
    double const * nu = getVectorPointer( arg_in_2, sz, CMD " param nu" );
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
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('make_retrograde',obj): "
    CHECK_IN(2);
    CHECK_OUT(0);
    Astro * ptr = DATA_GET( arg_in_1 );
    ptr->make_retrograde();
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_make_not_retrograde(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('make_not_retrograde',obj): "
    CHECK_IN(2);
    CHECK_OUT(0);
    Astro * ptr = DATA_GET( arg_in_1 );
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
    #define CMD "AstroMexWrapper('retrograde',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    setScalarValue( arg_out_0, ptr->retrograde() ? 1 : 0 );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_normal(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('normal',obj): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type * n = createMatrixValue( arg_out_0, 3, 1 );
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
    #define CMD "AstroMexWrapper('local_frame',obj,t): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type t = getScalarValue( arg_in_2, CMD " param t" );
    GC_namespace::real_type * frame = createMatrixValue( arg_out_0, 3, 3 );
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
    #define CMD "AstroMexWrapper('local_frame_by_L',obj,L): "
    CHECK_IN(3);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type L = getScalarValue( arg_in_2, CMD " param L" );
    GC_namespace::real_type * frame = createMatrixValue( arg_out_0, 3, 3 );
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
    #define CMD "AstroMexWrapper('ellipse_frame',obj): "
    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type * frame = createMatrixValue( arg_out_0, 3, 3 );
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
    #define CMD "AstroMexWrapper('eval_E',obj,t): "
    CHECK_IN(3);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type t = getScalarValue( arg_in_2, CMD " param t" );
    MEX_ASSERT2(
      nrhs > 0 && nrhs < 4,
      CMD "Expected 1,2,3 or 4  argument(s), nrhs = {}\n", nrhs
    );

    GC_namespace::real_type E[4];
    ptr->eval_E(t,E,nrhs-1);
    setScalarValue( arg_out_0, E[0] );
    if ( nrhs > 1 ) setScalarValue( arg_out_1, E[1] );
    if ( nrhs > 2 ) setScalarValue( arg_out_2, E[2] );
    if ( nrhs > 3 ) setScalarValue( arg_out_3, E[3] );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_eval_L(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "AstroMexWrapper('eval_L',obj,t): "
    CHECK_IN(3);
    Astro * ptr = DATA_GET( arg_in_1 );
    GC_namespace::real_type t = getScalarValue( arg_in_2, CMD " param t" );
    MEX_ASSERT2(
      nrhs > 0 && nrhs < 4,
      CMD "Expected 1,2,3 or 4  argument(s), nrhs = {}\n", nrhs
    );

    GC_namespace::real_type L[4];
    ptr->eval_L(t,L,nrhs-1);
    setScalarValue( arg_out_0, L[0] );
    if ( nrhs > 1 ) setScalarValue( arg_out_1, L[1] );
    if ( nrhs > 2 ) setScalarValue( arg_out_2, L[2] );
    if ( nrhs > 3 ) setScalarValue( arg_out_3, L[3] );
    #undef CMD
  }

  #if 0
  integer
  Lambert(
    real_type const R1[3],
    real_type const R2[3],
    real_type       tf_in, // tempo di volo
    integer         m_in,  // numero di rivoluzioni
    real_type       muC,
    real_type       V1[3],
    real_type       V2[3]
  );

  integer
  Lambert_Lancaster_Blanchard(
    real_type const r1vec[3],
    real_type const r2vec[3],
    real_type       tf_in,
    integer         m_in,
    real_type       muC,
    real_type       V1[3],
    real_type       V2[3]
  );

  void
  Lambert_minmax_distances(
    real_type const r1vec[3],
    real_type       r1,
    real_type const r2vec[3],
    real_type       r2,
    real_type       dth,
    real_type       a,
    real_type const V1[3],
    real_type const V2[3],
    integer         m,
    real_type       muC,
    real_type &     minimum_distance,
    real_type &     maximum_distance
  );

  #endif

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static std::map<std::string,DO_CMD> cmd_to_fun = {
    {"new",do_new},
    {"new_handle",do_new_handle},
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
      MEX_ASSERT( mxIsChar(arg_in_0), "First argument must be a string" );
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
