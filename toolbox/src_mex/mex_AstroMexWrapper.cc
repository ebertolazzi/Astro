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
    #define MEX_ERROR_MESSAGE_001 "AstroMexWrapper('new'[,'name'])"
    #define CMD MEX_ERROR_MESSAGE_001

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
    #define MEX_ERROR_MESSAGE_002 "AstroMexWrapper('delete',obj)"
    #define CMD MEX_ERROR_MESSAGE_002
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

    #define MEX_ERROR_MESSAGE_003 "AstroMexWrapper('copy',obj)"
    #define CMD MEX_ERROR_MESSAGE_003

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
    #define MEX_ERROR_MESSAGE_004 "res = AstroMexWrapper('do_name',obj)"
    #define CMD MEX_ERROR_MESSAGE_004

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
    #define MEX_ERROR_MESSAGE_005 "ok = AstroMexWrapper('check_for_consistency',obj)"
    #define CMD MEX_ERROR_MESSAGE_005

    CHECK_IN(2);
    CHECK_OUT(1);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::to_mxArray(ptr->check_for_consistency(),arg_out_0);
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_setup_Keplerian(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_006 "AstroMexWrapper('setup_Keplerian',obj,[name,t0,a,e,Omega,omega,i,M0,muS]||struct)"
    #define CMD MEX_ERROR_MESSAGE_006

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
    #define MEX_ERROR_MESSAGE_007 "AstroMexWrapper('setup_Equinoctial',obj,[name,t0,p,f,g,h,k,retrograde,L0,muS]||struct)"
    #define CMD MEX_ERROR_MESSAGE_007

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
      if ( R.get_type() == GC_namespace::GC_type::BOOL ) {
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
    #define MEX_ERROR_MESSAGE_008 "AstroMexWrapper('setup_PV',obj,name,P,V,t0,muS)"
    #define CMD MEX_ERROR_MESSAGE_008

    CHECK_IN(7);
    CHECK_OUT(0);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    UTILS_ASSERT0( mxIsChar(arg_in_2), CMD ": param name must be a string" );
    string n = mxArrayToString(arg_in_2);
    mwSize szP, szV;
    GC_namespace::real_type const * PP = Utils::mex_vector_pointer( arg_in_3, szP, CMD ": param P" );
    UTILS_ASSERT( szP == 3, CMD ": size |P| = {} expected 3\n", szP );
    GC_namespace::real_type const * VV = Utils::mex_vector_pointer( arg_in_4, szV, CMD ": param V" );
    UTILS_ASSERT( szV == 3, CMD ": size |V| = {} expected 3\n", szV );
    GC_namespace::real_type t0  = Utils::mex_get_scalar_value( arg_in_5, CMD ": param t0" );
    GC_namespace::real_type muS = Utils::mex_get_scalar_value( arg_in_6, CMD ": param muS" );
    dvec3_t P, V;
    P << PP[0], PP[1], PP[2];
    V << VV[0], VV[1], VV[2];
    ptr->setup_using_point_and_velocity( n, P, V, muS, t0 );
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_get_PV(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_009 "AstroMexWrapper('setup_PV',obj,t0)"
    #define CMD MEX_ERROR_MESSAGE_009

    CHECK_IN(3);
    CHECK_OUT(2);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    GC_namespace::real_type t0 = Utils::mex_get_scalar_value( arg_in_2, CMD ": param t0" );

    double * P = Utils::mex_create_matrix_value( arg_out_0, 3, 1 );
    double * V = Utils::mex_create_matrix_value( arg_out_1, 3, 1 );

    ptr->position( t0, P );
    ptr->velocity( t0, V );

    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  void
  do_info(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define MEX_ERROR_MESSAGE_010 "info = AstroMexWrapper('info',obj)"
    #define CMD MEX_ERROR_MESSAGE_010

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
    #define MEX_ERROR_MESSAGE_011 "AstroMexWrapper('print_info',obj)"
    #define CMD MEX_ERROR_MESSAGE_011

    CHECK_IN(2);
    CHECK_OUT(0);
    Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
    mexPrintf("%s\n",ptr->info().c_str());
    #undef CMD
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #include "mex_AstroPosition.cxx"
  #include "mex_AstroOrbital.cxx"
  #include "mex_AstroJacobian.cxx"

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
    {"get_PV",do_get_PV},
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
    {"radius",do_radius},
    {"absolute_velocity",do_absolute_velocity},
    {"muS",do_muS},
    {"period",do_period},
    {"apoapsis",do_apoapsis},
    {"periapsis",do_periapsis},
    {"orbit_energy",do_orbit_energy},
    {"radius_by_L",do_radius_by_L},
    {"radius_by_L_D",do_radius_by_L_D},
    {"radius_by_L_DD",do_radius_by_L_DD},
    {"M0",do_M0},
    {"theta0",do_theta0},
    {"E0_angle",do_E0_angle},
    {"E_angle",do_E_angle},
    {"H0_angle",do_H0_angle},
    {"H_angle",do_H_angle},
    {"position_EQ_jacobian",do_position_EQ_jacobian},
    {"position_EQ_jacobian_FD",do_position_EQ_jacobian_FD},
    {"position0_EQ_jacobian",do_position0_EQ_jacobian},
    {"velocity_EQ_jacobian",do_velocity_EQ_jacobian},
    {"velocity_EQ_jacobian_FD",do_velocity_EQ_jacobian_FD},
    {"velocity0_EQ_jacobian",do_velocity0_EQ_jacobian},
    {"absolute_velocity_by_angle",do_absolute_velocity_by_angle},

    {"L_orbital_EQ_gradient",do_L_orbital_EQ_gradient},
    {"M0_EQ_gradient",do_M0_EQ_gradient},
    {"theta0_EQ_gradient",do_theta0_EQ_gradient},
    {"E0_EQ_gradient",do_E0_EQ_gradient},
    {"E_EQ_gradient",do_E_EQ_gradient},
    {"H0_EQ_gradient",do_H0_EQ_gradient},
    {"H_EQ_gradient",do_H_EQ_gradient},
    {"radius_EQ_gradient",do_radius_EQ_gradient},
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
MEX_ERROR_MESSAGE_001 "\n" \
MEX_ERROR_MESSAGE_002 "\n" \
MEX_ERROR_MESSAGE_003 "\n" \
MEX_ERROR_MESSAGE_004 "\n" \
MEX_ERROR_MESSAGE_005 "\n" \
MEX_ERROR_MESSAGE_006 "\n" \
MEX_ERROR_MESSAGE_007 "\n" \
MEX_ERROR_MESSAGE_008 "\n" \
MEX_ERROR_MESSAGE_009 "\n" \
MEX_ERROR_MESSAGE_010 "\n" \
MEX_ERROR_MESSAGE_011 "\n" \
MEX_ERROR_MESSAGE_101 "\n" \
MEX_ERROR_MESSAGE_102 "\n" \
MEX_ERROR_MESSAGE_103 "\n" \
MEX_ERROR_MESSAGE_104 "\n" \
MEX_ERROR_MESSAGE_105 "\n" \
MEX_ERROR_MESSAGE_106 "\n" \
MEX_ERROR_MESSAGE_107 "\n" \
MEX_ERROR_MESSAGE_201 "\n" \
MEX_ERROR_MESSAGE_202 "\n" \
MEX_ERROR_MESSAGE_203 "\n" \
MEX_ERROR_MESSAGE_204 "\n" \
MEX_ERROR_MESSAGE_205 "\n" \
MEX_ERROR_MESSAGE_206 "\n" \
MEX_ERROR_MESSAGE_207 "\n" \
MEX_ERROR_MESSAGE_208 "\n" \
MEX_ERROR_MESSAGE_209 "\n" \
MEX_ERROR_MESSAGE_210 "\n" \
MEX_ERROR_MESSAGE_211 "\n" \
MEX_ERROR_MESSAGE_212 "\n" \
MEX_ERROR_MESSAGE_213 "\n" \
MEX_ERROR_MESSAGE_214 "\n" \
MEX_ERROR_MESSAGE_215 "\n" \
MEX_ERROR_MESSAGE_216 "\n" \
MEX_ERROR_MESSAGE_217 "\n" \
MEX_ERROR_MESSAGE_218 "\n" \
MEX_ERROR_MESSAGE_219 "\n" \
MEX_ERROR_MESSAGE_220 "\n" \
MEX_ERROR_MESSAGE_221 "\n" \
MEX_ERROR_MESSAGE_222 "\n" \
MEX_ERROR_MESSAGE_223 "\n" \
MEX_ERROR_MESSAGE_224 "\n" \
MEX_ERROR_MESSAGE_225 "\n" \
MEX_ERROR_MESSAGE_226 "\n" \
MEX_ERROR_MESSAGE_227 "\n" \
MEX_ERROR_MESSAGE_228 "\n" \
MEX_ERROR_MESSAGE_229 "\n" \
MEX_ERROR_MESSAGE_230 "\n" \
MEX_ERROR_MESSAGE_231 "\n" \
MEX_ERROR_MESSAGE_232 "\n" \
MEX_ERROR_MESSAGE_233 "\n" \
MEX_ERROR_MESSAGE_234 "\n" \
MEX_ERROR_MESSAGE_235 "\n" \
MEX_ERROR_MESSAGE_236 "\n" \
MEX_ERROR_MESSAGE_237 "\n" \
MEX_ERROR_MESSAGE_238 "\n" \
MEX_ERROR_MESSAGE_301 "\n" \
MEX_ERROR_MESSAGE_302 "\n" \
MEX_ERROR_MESSAGE_303 "\n" \
MEX_ERROR_MESSAGE_304 "\n" \
MEX_ERROR_MESSAGE_305 "\n" \
MEX_ERROR_MESSAGE_306 "\n" \
MEX_ERROR_MESSAGE_307 "\n" \
MEX_ERROR_MESSAGE_308 "\n" \
MEX_ERROR_MESSAGE_309 "\n" \
MEX_ERROR_MESSAGE_310 "\n" \
MEX_ERROR_MESSAGE_311 "\n" \
MEX_ERROR_MESSAGE_312 "\n" \
MEX_ERROR_MESSAGE_313 "\n" \
MEX_ERROR_MESSAGE_314 "\n" \
MEX_ERROR_MESSAGE_315 "\n" \
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
