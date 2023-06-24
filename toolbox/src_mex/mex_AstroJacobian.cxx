// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static
void
do_L_orbital_EQ_gradient(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_301 "L_grad = AstroMexWrapper('L_orbital_EQ_gradient',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_301
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
do_M0_EQ_gradient(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_302 "M0_grad = AstroMexWrapper('M0_EQ_gradient',obj)"
  #define CMD MEX_ERROR_MESSAGE_302
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
do_theta0_EQ_gradient(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_303 "theta0_grad = AstroMexWrapper('theta0_EQ_gradient',obj)"
  #define CMD MEX_ERROR_MESSAGE_303
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
do_E0_EQ_gradient(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_304 "E0_grad = AstroMexWrapper('E0_EQ_gradient',obj)"
  #define CMD MEX_ERROR_MESSAGE_304
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
do_H0_EQ_gradient(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_305 "H0_grad = AstroMexWrapper('H0_EQ_gradient',obj)"
  #define CMD MEX_ERROR_MESSAGE_305
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
do_E_EQ_gradient(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_306 "E_grad = AstroMexWrapper('E_EQ_gradient',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_306
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
do_H_EQ_gradient(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_307 "H_grad = AstroMexWrapper('H_EQ_gradient',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_307
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
  #define MEX_ERROR_MESSAGE_308 "R_grad = AstroMexWrapper('radius_EQ_gradient',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_308
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
  #define MEX_ERROR_MESSAGE_309 "P_grad = AstroMexWrapper('position_EQ_jacobian',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_309
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
  #define MEX_ERROR_MESSAGE_310 "P0_grad = AstroMexWrapper('position0_EQ_jacobian',obj)"
  #define CMD MEX_ERROR_MESSAGE_310
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
  #define MEX_ERROR_MESSAGE_311 "P_grad = AstroMexWrapper('position_EQ_jacobian_FD',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_311
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
do_absolute_velocity_EQ_gradient(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {
  #define MEX_ERROR_MESSAGE_312 "V_grad = AstroMexWrapper('absolute_velocity_EQ_gradient',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_312
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
  #define MEX_ERROR_MESSAGE_313 "V0_grad = AstroMexWrapper('velocity0_EQ_jacobian',obj)"
  #define CMD MEX_ERROR_MESSAGE_313
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
  #define MEX_ERROR_MESSAGE_314 "V_jacobian = AstroMexWrapper('velocity_EQ_jacobian',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_314
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
  #define MEX_ERROR_MESSAGE_315 "V_jacobian_FD = AstroMexWrapper('velocity_EQ_jacobian_FD',obj,t)"
  #define CMD MEX_ERROR_MESSAGE_315
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
