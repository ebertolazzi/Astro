/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: Utils_NelderMead.cc
///

#include "Utils_NelderMead.hh"

//#define UTILS_DEBUG

#ifdef UTILS_DEBUG
  #define REF(A,...)   A.coeffRef(__VA_ARGS__)
  #define VALUE(A,...) A.coeff(__VA_ARGS__)
#else
  #define REF(A,...)   A(__VA_ARGS__)
  #define VALUE(A,...) A(__VA_ARGS__)
#endif

namespace Utils {

  // =================================================================
  // set_tolerance
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::set_tolerance( Real tol ) {
    UTILS_ASSERT(
      tol > 0,
      "NelderMead::set_tolerance({}) argument must be >0\n", tol
    );
    m_tolerance = tol;
  }

  // =================================================================
  // set_max_iterations
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::set_max_iterations( integer mit ) {
    UTILS_ASSERT(
      mit > 0,
      "NelderMead::set_max_iterations({}) argument must be >0\n", mit
    );
    m_max_iteration = mit;
  }

  // =================================================================
  // set_max_fun_evaluation
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::set_max_fun_evaluation( integer mfev ) {
    UTILS_ASSERT(
      mfev > 0,
      "NelderMead::set_max_fun_evaluation({}) argument must be >0\n", mfev
    );
    m_max_fun_evaluation = mfev;
  }

  // =================================================================
  // allocate
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::allocate( integer n ) {
    m_dim           = n;
    m_r_dim         = Real(1)/Real(n);
    m_dim_factorial = 1;
    for ( integer i = 2; i <= n+1; ++i ) m_dim_factorial *= i;

    m_dim_regular_simplex_volume = (sqrt(Real(n+1))/Real(m_dim_factorial))/std::pow(Real(2),Real(n)/2);

    integer ntot = (3*n+10)*n+2;
    m_base_value.reallocate( ntot );

    size_t N  = size_t(n);
    size_t N1 = size_t(n+1);

    new (&this->m_f)      MapVec( m_base_value( N1 ),    n+1      );
    new (&this->m_p)      MapMat( m_base_value( N*N1 ),  n,   n+1 );
    new (&this->m_p_work) MapMat( m_base_value( N*N ),   n,   n   );
    new (&this->m_dist)   MapMat( m_base_value( N1*N1 ), n+1, n+1 );

    new (&this->m_psum)   MapVec( m_base_value( N ), n );
    new (&this->m_f_work) MapVec( m_base_value( N ), n );
    new (&this->m_grad)   MapVec( m_base_value( N ), n );
    new (&this->m_pr)     MapVec( m_base_value( N ), n );
    new (&this->m_pe)     MapVec( m_base_value( N ), n );
    new (&this->m_pc)     MapVec( m_base_value( N ), n );

    m_base_value.must_be_empty("NelderMead<Real>::allocate");
  }

  // =================================================================
  // setup
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::setup( integer dim, NMFunc & fun, Console const * console ) {
    // HJPatternSearch The constructor initialize the solver
    // parameters and check the inputs when the class is instanciated.
    m_fun     = fun;
    m_console = console;
    allocate( dim );
  }

  // =================================================================
  // replace_point
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::replace_point(
    Real           fj,
    MapVec const & pj,
    integer        jpos
  ) {
    REF(m_f,jpos)           = fj;
    m_psum                 += pj - m_p.col(jpos);
    m_p.col(jpos).noalias() = pj;
    this->dist_update( jpos );
  }

  // =================================================================
  // spendley
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::spendley( Real const X0[], Real delta ) {
    Real t1 = sqrt( Real(m_dim+1) ) - Real(1);
    Real t2 = Real(m_dim)*sqrt(Real(2));
    Real p  = delta * (m_dim+t1)/t2;
    Real q  = delta * t1/t2;
    for ( integer i = 0; i < m_dim; ++i ) REF(m_p,i,0) = X0[i];
    REF(m_f,0) = this->eval_function( m_p.col(0).data() );
    for ( integer i = 0; i < m_dim; ++i ) {
      m_p.col(i+1)    = m_p.col(0).array()+p;
      REF(m_p,i,i+1) += q;
      REF(m_f,i+1)    = this->eval_function( m_p.col(i+1).data() );
    }
  }

  // =================================================================
  // diamond
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::diamond( Real const X0[], Real delta ) {
    for ( integer i = 0; i < m_dim; ++i ) m_p.row(i).fill(X0[i]);
    REF(m_f,0) = this->eval_function( X0 );
    for ( integer i = 0; i < m_dim; ++i ) {
      Real & vi = REF(m_p,i,i+1);
      Real * cl = m_p.col(i+1).data();
      vi = X0[i]+delta; Real fp = this->eval_function( cl );
      vi = X0[i]-delta; Real fm = this->eval_function( cl );
      if ( fp < fm ) {
        REF(m_f,i+1) = fp;
        vi = X0[i]+delta;
      } else {
        REF(m_f,i+1) = fm;
      }
    }
  }

  // =================================================================
  // grad_update
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::grad_update( integer j0 ) {
    integer j = 0;
    for ( integer i = 0; i <= m_dim; ++i  ) {
      if ( i != j0 ) {
        m_p_work.col(j).noalias() = m_p.col(i)-m_p.col(j0);
        REF(m_f_work,j)           = VALUE(m_f,i)-VALUE(m_f,j0);
        ++j;
      }
    }
    m_lu.compute( m_p_work.transpose() );
    m_simplex_volume = std::abs(m_lu.determinant()) / m_dim_factorial;
    m_grad.noalias() = m_lu.solve( m_f_work );
  }

  // =================================================================
  // dist_init
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::dist_init() {
    for ( integer i = 0; i < m_dim; ++i  ) {
      REF(m_dist,i,i) = 0;
      for ( integer j = i+1; j <= m_dim; ++j ) {
        REF(m_dist,i,j) = ( m_p.col(i) - m_p.col(j) ).norm();
        REF(m_dist,j,i) = VALUE(m_dist,i,j);
      }
    }
    m_diameter = m_dist.maxCoeff();
    grad_update( 0 );
  }

  // =================================================================
  // dist_update
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::dist_update( integer j ) {
    for ( integer i = 0; i <= m_dim; ++i  ) {
      if ( i != j ) {
        REF(m_dist,i,j) = ( m_p.col(i) - m_p.col(j) ).norm();
        REF(m_dist,j,i) = VALUE(m_dist,i,j);
      }
    }
    m_diameter = m_dist.maxCoeff();
    grad_update( m_low );
  }

  // =================================================================
  // shrink
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::shrink() {
    Real c1 = 1-m_sigma;
    Real c2 = m_sigma;
    for ( integer i = 0; i <= m_dim; ++i ) {
      if ( i != m_low ) {
        m_p.col(i) = c1 * m_p.col(i) + c2 * m_p.col(m_low);
        REF(m_f,i) = this->eval_function( m_p.col(i).data() );
      }
    }
    m_psum.noalias() = m_p.col(m_dim);
    for ( integer i = 0 ; i < m_dim; ++i ) m_psum.noalias() += m_p.col(i);
    m_dist           *= c1;
    m_diameter       *= c1;
    m_simplex_volume *= std::pow(c1,m_dim);
  }

  // =================================================================
  // extrapolate
  // =================================================================

  template <typename Real>
  Real
  NelderMead<Real>::extrapolate( Real alpha, integer j, MapVec & pe ) const {
    //
    // Extrapolates by a factor alpha through the face of the simplex
    // across from the high point.
    //
    // pface   = (sum(p)-pfrom)/N
    // pcenter = sum(p)/(N+1);
    // pe      = pface + alpha*(pface-pfrom)
    //         = pface*(1+alpha)-pfrom*alpha
    //         = (sum(p)-pfrom)/N*(1+alpha)-pfrom*alpha
    //         = sum(p)/N*(1+alpha)-pfrom*(alpha+(1+alpha)/N)
    //         = sum(p)/(N+1) * (1+alpha)*((N+1)/N) - pfrom*[1/N+alpha*(N+1)/N]
    //         = pcenter * (1+alpha)*((N+1)/N) - pfrom*(1/N-(N+1)/N+(1+alpha)*(N+1)/N)
    //         = pcenter * (1+alpha)*((N+1)/N) - pfrom*(-1+(1+alpha)*(N+1)/N)
    //
    // let beta = (1+alpha)*((N+1)/N);
    //
    // pe = pcenter * beta + pfrom * (1-beta)
    //
    Real beta1   = (1+alpha)/m_dim;
    Real beta    = (m_dim+1)*beta1;
    pe.noalias() = beta1 * m_psum + (1-beta) * m_p.col(j);
    return this->eval_function( pe.data() );
  }

  // =================================================================
  // info
  // =================================================================

  template <typename Real>
  string
  NelderMead<Real>::info() const {
    string res = "";
    if ( m_verbose > 0 ) {
      Real Vd = std::pow(m_simplex_volume,m_r_dim);
      string const line = "-------------------------------------------------------------------------\n";
      if ( m_verbose > 1 ) res += line;
      res += fmt::format(
        "#it={:5} #f={:5} {:15} f(x)={:12} grad(x)={:12} V^(1/d)/D={} [{}/{}]\n",
        m_iteration_count, m_fun_evaluation_count,
        to_string(m_which_step),
        fmt::format("{:.6}",VALUE(m_f,m_low)),
        fmt::format("{:.6}",m_grad.template lpNorm<Eigen::Infinity>()),
        fmt::format("{:.4}",Vd/m_diameter),
        fmt::format("{:.4}",Vd),
        fmt::format("{:.4}",m_diameter)
      );
    }
    return res;
  }

  // =================================================================
  // message
  // =================================================================

  template <typename Real>
  void
  NelderMead<Real>::message( Real rtol ) const {
    if ( m_console == nullptr ) return;
    Real Vd = std::pow(m_simplex_volume,m_r_dim);
    string msg = fmt::format(
      "#it={:<4} #f={:<5} {:<12} f(x)={:<12} grad(x)={:<12}"
      " |err|={:<10} V^(1/d)/D={:<10} [{}/{}]\n",
      m_iteration_count, m_fun_evaluation_count, to_string(m_which_step),
      fmt::format("{:.6}",VALUE(m_f,m_low)),
      fmt::format("{:.6}",m_grad.template lpNorm<Eigen::Infinity>()),
      fmt::format("{:.4}",rtol),
      fmt::format("{:.4}",Vd/m_diameter),
      fmt::format("{:.4}",Vd),
      fmt::format("{:.4}",m_diameter)
    );
    /*
    msg += fmt::format(
     "x={:8} y={:8}\n",
     fmt::format("{:.8}",m_p(0,m_low)),
     fmt::format("{:.8}",m_p(1,m_low))
    );
    */
    m_console->message( msg, 2 );
  }

  // =================================================================
  // search
  // =================================================================

  template <typename Real>
  bool
  NelderMead<Real>::search() {

    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    this->dist_init();
    m_psum.noalias() = m_p.col(m_dim);
    for ( integer i = 0 ; i < m_dim; ++i ) m_psum.noalias() += m_p.col(i);
    m_which_step = NM_move::INIT;

    for (; m_iteration_count <= m_max_iteration; ++m_iteration_count ) {

      if ( m_fun_evaluation_count >= m_max_fun_evaluation ) return false;

      /*
      // Determine which point is the highest (worst),
      // next-highest, and lowest (best).
      */

      m_low   = 0;
      m_0high = 1;
      m_high  = 2;
      if ( VALUE(m_f,m_low)   > VALUE(m_f,m_0high) ) std::swap( m_low,   m_0high );
      if ( VALUE(m_f,m_low)   > VALUE(m_f,m_high)  ) std::swap( m_low,   m_high  );
      if ( VALUE(m_f,m_0high) > VALUE(m_f,m_high)  ) std::swap( m_0high, m_high  );
      for ( integer i = 3 ; i <= m_dim; ++i ) {
        if ( VALUE(m_f,i) < VALUE(m_f,m_low) ) {
          m_low = i; // new minima
        } else if ( VALUE(m_f,i) > VALUE(m_f,m_high) ) {
          m_0high = m_high;
          m_high  = i;
        } else if ( VALUE(m_f,i) > VALUE(m_f,m_0high) ) {
          m_0high = i;
        }
      }

      Real f0  = VALUE(m_f,m_low);
      Real fn  = VALUE(m_f,m_0high);
      Real fn1 = VALUE(m_f,m_high);

      /*
      // Compute the fractional range from highest to lowest
      // and return if satisfactory.
      */
      Real rtol = abs( fn1 - f0 ) / ( abs(f0) + m_tolerance );

      if ( m_verbose > 0 ) message( rtol );
      if ( rtol < m_tolerance || m_diameter < m_tolerance ) return true;

      Real ratio = std::pow(m_simplex_volume/m_dim_regular_simplex_volume,m_r_dim)/m_diameter;
      if ( ratio <= m_volume_tolerance ) {
        this->diamond( m_p.col(m_low).data(), m_diameter*m_volume_tolerance );
        this->dist_init();
        m_psum.noalias() = m_p.col(m_dim);
        for ( integer i = 0 ; i < m_dim; ++i ) m_psum.noalias() += m_p.col(i);
        m_which_step = NM_move::RESTART;
        continue;
      }

      //
      // Begin a new iteration.
      // First extrapolate by a factor m_alpha (default: −1 a reflection)
      // through the face of the simplex across from the high point,
      // i.e., reflect the simplex from the high point.
      //
      Real fr = this->reflect( m_pr );
      if ( fr < f0 ) {
        Real fe = this->expand( m_pe );
        if ( fe < fr ) {
          m_which_step = NM_move::EXPAND_FE;
          this->replace_point( fe, m_pe, m_high );
        } else {
          m_which_step = NM_move::EXPAND_FR;
          this->replace_point( fr, m_pr, m_high );
        }
      } else if ( fr < fn ) {
        m_which_step = NM_move::REFLECT;
        this->replace_point( fr, m_pr, m_high );
      } else if ( fr < fn1 ) {
        Real fc = this->outside( m_pc );
        if ( fc < fr ) {
          m_which_step = NM_move::CONTRACT_O;
          this->replace_point( fc, m_pc, m_high );
        } else {
          m_which_step = NM_move::REFLECT;
          this->replace_point( fr, m_pr, m_high );
        }
      } else {
        Real fc = this->inside( m_pc );
        if ( fc < fn1 ) {
          m_which_step = NM_move::CONTRACT_I;
          this->replace_point( fc, m_pc, m_high );
        } else {
          m_which_step = NM_move::SHRINK;
          this->shrink();
        }
      }
    }
    return false;
  }

  // =================================================================
  // run
  // =================================================================

  template <typename Real>
  bool
  NelderMead<Real>::run( Real const X0[], Real delta ) {
    this->diamond( X0, delta );
    return this->search();
  }

  // =================================================================

  template class NelderMead<double>;
  template class NelderMead<float>;

}

///
/// eof: Utils_NelderMead.cc
///
