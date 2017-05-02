#include "mex.h"
#include "astro.hh"

#include <string>
#include <vector>
#include <algorithm>

#ifdef OLD_VERSION
#define mwSize  int
#define mwIndex int
#endif

namespace astro {

  using namespace std ;

  bool
  is_matrix( mxArray const * mx,
             int_type  & nrow,
             int_type  & ncol ) ;

  bool
  is_vector( mxArray const * mx, int_type & dim ) ;

  bool
  is_scalar( mxArray const * mx, real_type & d ) ;

  bool
  is_integer( mxArray const * mx, int_type & d ) ;

  // ===========================================================================

  void
  to_mxArray( bool const & val_in, mxArray * & mx ) ;

  void
  to_mxArray( int_type const & val, mxArray * & mx ) ;

  void
  to_mxArray( real_type const & val, mxArray * & mx ) ;

  void
  to_mxArray( std::string const & val, mxArray * & mx ) ;

  void
  to_mxArray( std::vector<int_type> const & val, mxArray * & mx );

  void
  to_mxArray( std::vector<real_type> const & val, mxArray * & mx ) ;

  void
  to_mxArray( std::vector<std::string> const & val, mxArray * & mx );

}
