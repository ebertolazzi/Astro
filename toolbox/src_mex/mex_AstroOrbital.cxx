// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static
void
do_true_anomaly(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_201 "res = AstroMexWrapper('true_anomaly',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_201
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
  #define MEX_ERROR_MESSAGE_202 "res = AstroMexWrapper('p_orbital',obj)"
  #define CMD MEX_ERROR_MESSAGE_202
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
  #define MEX_ERROR_MESSAGE_203 "res = AstroMexWrapper('f_orbital',obj)"
  #define CMD MEX_ERROR_MESSAGE_203
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
  #define MEX_ERROR_MESSAGE_204 "res = AstroMexWrapper('g_orbital',obj)"
  #define CMD MEX_ERROR_MESSAGE_204
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
  #define MEX_ERROR_MESSAGE_205 "res = AstroMexWrapper('h_orbital',obj)"
  #define CMD MEX_ERROR_MESSAGE_205
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
  #define MEX_ERROR_MESSAGE_206 "AstroMexWrapper('k_orbital',obj)"
  #define CMD MEX_ERROR_MESSAGE_206
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
  #define MEX_ERROR_MESSAGE_207 "L = AstroMexWrapper('L_orbital',obj,[t0,dt] or [t])"
  #define CMD MEX_ERROR_MESSAGE_207
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
  #define MEX_ERROR_MESSAGE_208 "DLDt = AstroMexWrapper('L_orbital_D',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_208
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
  #define MEX_ERROR_MESSAGE_209 "D2LD2t = AstroMexWrapper('L_orbital_DD',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_209
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
  #define MEX_ERROR_MESSAGE_210 "res = AstroMexWrapper('latitude_of_periapsis',obj)"
  #define CMD MEX_ERROR_MESSAGE_210
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
  #define MEX_ERROR_MESSAGE_211 "res = AstroMexWrapper('latitude_of_apoapsis',obj)"
  #define CMD MEX_ERROR_MESSAGE_211
  CHECK_IN(2);
  CHECK_OUT(1);
  Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
  GC_namespace::to_mxArray(ptr->latitude_of_apoapsis(),arg_out_0);
  #undef CMD
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static
void
do_number_of_revolution(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_212 "n = AstroMexWrapper('number_of_revolution',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_212
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
  #define MEX_ERROR_MESSAGE_213 "time = AstroMexWrapper('time_from_L_angle',obj,t,L)"
  #define CMD MEX_ERROR_MESSAGE_213
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
do_period(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_214 "period = AstroMexWrapper('period',obj)"
  #define CMD MEX_ERROR_MESSAGE_214
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
  #define MEX_ERROR_MESSAGE_215 "a = AstroMexWrapper('apoapsis',obj)"
  #define CMD MEX_ERROR_MESSAGE_215
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
  #define MEX_ERROR_MESSAGE_216 "p = AstroMexWrapper('periapsis',obj)"
  #define CMD MEX_ERROR_MESSAGE_216
  CHECK_IN(2);
  CHECK_OUT(1);
  Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
  GC_namespace::to_mxArray(ptr->periapsis(),arg_out_0);
  #undef CMD
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static
void
do_mean_anomaly(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_217 "res = AstroMexWrapper('mean_anomaly',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_217
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
do_muS(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_218 "muS = AstroMexWrapper('muS',obj)"
  #define CMD MEX_ERROR_MESSAGE_218
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
  #define MEX_ERROR_MESSAGE_219 "E = AstroMexWrapper('orbit_energy',obj)"
  #define CMD MEX_ERROR_MESSAGE_219
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
  #define MEX_ERROR_MESSAGE_220 "R = AstroMexWrapper('radius_by_L',obj,L)"
  #define CMD MEX_ERROR_MESSAGE_220
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
  #define MEX_ERROR_MESSAGE_221 "DRDL = AstroMexWrapper('radius_by_L_D',obj,L)"
  #define CMD MEX_ERROR_MESSAGE_221
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
  #define MEX_ERROR_MESSAGE_222 "D2RD2L = AstroMexWrapper('radius_by_L_DD',obj,L)"
  #define CMD MEX_ERROR_MESSAGE_222
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
do_M0(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_223 "M0 = AstroMexWrapper('M0',obj)"
  #define CMD MEX_ERROR_MESSAGE_223
  CHECK_IN(2);
  CHECK_OUT(1);
  Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
  Utils::mex_set_scalar_value( arg_out_0, ptr->M0_orbital() );
  #undef CMD
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static
void
do_theta0(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_224 "theta0 = AstroMexWrapper('theta0',obj)"
  #define CMD MEX_ERROR_MESSAGE_224
  CHECK_IN(2);
  CHECK_OUT(1);
  Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
  Utils::mex_set_scalar_value( arg_out_0, ptr->theta0_orbital() );
  #undef CMD
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static
void
do_E0_angle(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_225 "E0 = AstroMexWrapper('E0_angle',obj)"
  #define CMD MEX_ERROR_MESSAGE_225
  CHECK_IN(2);
  CHECK_OUT(1);
  Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
  Utils::mex_set_scalar_value( arg_out_0, ptr->E0_angle() );
  #undef CMD
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static
void
do_H0_angle(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_226 "H0 = AstroMexWrapper('H0_angle',obj)"
  #define CMD MEX_ERROR_MESSAGE_226
  CHECK_IN(2);
  CHECK_OUT(1);
  Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
  Utils::mex_set_scalar_value( arg_out_0, ptr->H0_angle() );
  #undef CMD
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static
void
do_E_angle(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_227 "E = AstroMexWrapper('E_angle',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_227
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
do_H_angle(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_228 "H = AstroMexWrapper('H_angle',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_228
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
do_L_from_true_anomaly(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_229 "L = AstroMexWrapper('L_from_true_anomaly',obj,nu)"
  #define CMD MEX_ERROR_MESSAGE_229
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
  #define MEX_ERROR_MESSAGE_230 "AstroMexWrapper('make_retrograde',obj)"
  #define CMD MEX_ERROR_MESSAGE_230
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
  #define MEX_ERROR_MESSAGE_231 "AstroMexWrapper('make_not_retrograde',obj)"
  #define CMD MEX_ERROR_MESSAGE_231
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
  #define MEX_ERROR_MESSAGE_232 "res = AstroMexWrapper('retrograde',obj)"
  #define CMD MEX_ERROR_MESSAGE_232
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
  #define MEX_ERROR_MESSAGE_233 "N = AstroMexWrapper('normal',obj)"
  #define CMD MEX_ERROR_MESSAGE_233
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
  #define MEX_ERROR_MESSAGE_234 "frame = AstroMexWrapper('local_frame',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_234
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
  #define MEX_ERROR_MESSAGE_235 "frame = AstroMexWrapper('local_frame_by_L',obj,L)"
  #define CMD MEX_ERROR_MESSAGE_235
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
  #define MEX_ERROR_MESSAGE_236 "frame = AstroMexWrapper('ellipse_frame',obj)"
  #define CMD MEX_ERROR_MESSAGE_236
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
  #define MEX_ERROR_MESSAGE_237 "E = AstroMexWrapper('eval_E',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_237
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
  #define MEX_ERROR_MESSAGE_238 "L = AstroMexWrapper('eval_L',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_238
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
