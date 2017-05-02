#include "common_Mex.hh"

namespace astro {

  using namespace std ;

  bool
  is_matrix( mxArray const * mx,
             int_type  & nrow,
             int_type  & ncol ) {
    bool ok = mxIsDouble(mx) &&
              !mxIsEmpty(mx) &&
              !mxIsComplex(mx) &&
              mxIsNumeric(mx) ;
    if ( ok ) {
      nrow = int_type(mxGetM(mx));
      ncol = int_type(mxGetN(mx));
    }
    return ok ;
  }

  bool
  is_vector( mxArray const * mx, int_type & dim ) {
    int_type nrow, ncol ;
    bool ok = is_matrix( mx, nrow, ncol ) ;
    ok  = min( nrow, ncol ) == 1 ;
    dim = max( nrow, ncol ) ;
    return ok ;
  }

  bool
  is_scalar( mxArray const * mx, real_type & d ) {
    int_type nrow, ncol ;
    bool ok = is_matrix( mx, nrow, ncol ) ;
    if ( ok ) ok = nrow == 1 && ncol == 1 ;
    if ( ok ) d = mxGetScalar(mx) ;
    return ok;
  }

  bool
  is_integer( mxArray const * mx, int_type & d ) {
    int_type nrow, ncol ;
    bool ok = is_matrix( mx, nrow, ncol ) ;
    if ( ok ) ok = nrow == 1 && ncol == 1 ;
    if ( ok ) {
      double res  = mxGetScalar(mx) ;
      double ires = floor(res);
      ok = std::abs(res-ires) < 1e-9 ;
      if ( ok ) d = int_type(ires);
    }
    return ok;
  }
  // ===========================================================================

  void
  to_mxArray( bool const & val_in, mxArray * & mx ) {
    mxLogical val = val_in ? 1 : 0 ;
    mx = mxCreateLogicalScalar(val) ;
  }

  void
  to_mxArray( int_type const & val, mxArray * & mx ) {
    mwSize dims[2] = {1,1} ;
    mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL) ;
    *(mwSize *)mxGetData(mx) = mwSize(val) ;
  }

  void
  to_mxArray( real_type const & val, mxArray * & mx ) {
    mx = mxCreateDoubleScalar(val) ;
  }

  void
  to_mxArray( std::string const & val, mxArray * & mx ) {
    mx = mxCreateString( val.c_str() ) ;
  }

  void
  to_mxArray( std::vector<int_type> const & val, mxArray * & mx ) {
    mwSize dims[2] = { mwSize(val.size()), 1 } ;
    mx = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL) ;
    int32_t * ptr = (int32_t*)mxGetData(mx) ;
    for ( mwSize i = 0 ; i < dims[0] ; ++i ) ptr[i] = int32_t(val[i]) ;
  }

  void
  to_mxArray( std::vector<real_type> const & val, mxArray * & mx ) {
    mwSize dims[2] = { mwSize(val.size()), 1 } ;
    mx = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL) ;
    double * ptr = mxGetPr(mx) ;
    for ( mwSize i = 0 ; i < dims[0] ; ++i ) ptr[i] = val[i] ;
  }

  void
  to_mxArray( std::vector<std::string> const & val, mxArray * & mx ) {
    mwSize dims[2] = { mwSize(val.size()), 1 } ;
    mx = mxCreateCellMatrix(dims[0], dims[1]) ;
    for( mwSize i = 0 ; i < dims[0] ; ++i )
      mxSetCell(mx,i,mxCreateString( val[i].c_str()) );
  }
}
