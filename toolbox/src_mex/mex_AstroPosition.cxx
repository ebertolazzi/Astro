
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static
void
do_position(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_101 "xyz|[x,y,z] = AstroMexWrapper('position',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_101
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
  #define MEX_ERROR_MESSAGE_102 "vxyz|[vx,vy,vz] = AstroMexWrapper('velocity',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_102
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
  #define MEX_ERROR_MESSAGE_103 "axyz|[ax,ay,az] = AstroMexWrapper('acceleration',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_103
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
  #define MEX_ERROR_MESSAGE_104 "jxyz|[jx,jy,jz] = AstroMexWrapper('jerk',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_104
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
do_radius(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_105 "radius = AstroMexWrapper('radius',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_105
  CHECK_IN(3);
  CHECK_OUT(1);
  Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
  GC_namespace::real_type t = Utils::mex_get_scalar_value( arg_in_2, CMD " param t" );
  GC_namespace::to_mxArray( ptr->radius(t), arg_out_0 );
  #undef CMD
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static
void
do_absolute_velocity(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_106 "v = AstroMexWrapper('absolute_velocity',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_106
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
do_absolute_velocity_by_angle(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_107 "V = AstroMexWrapper('absolute_velocity_by_angle',obj,L)"
  #define CMD MEX_ERROR_MESSAGE_107
  CHECK_IN(3);
  CHECK_OUT(1);
  Astro * ptr = Utils::mex_convert_mx_to_ptr<Astro>( arg_in_1 );
  GC_namespace::real_type L = Utils::mex_get_scalar_value( arg_in_2, CMD ": param L" );
  GC_namespace::to_mxArray(ptr->absolute_velocity_by_angle(L),arg_out_0);
  #undef CMD
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
