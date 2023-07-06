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
/// file: Utils_Poly.cc
///

#if defined(__llvm__) || defined(__clang__)
#pragma clang diagnostic ignored "-Wdeprecated-copy-with-dtor"
#endif

#include "Utils_Poly.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#define UTILS_POLY_INTANTIATE( REAL ) \
namespace Utils { \
  template class Poly<REAL>; \
  template class Sturm<REAL>; \
\
  template Poly<REAL> operator + ( Poly<REAL> const & a, Poly<REAL> const & b ); \
  template Poly<REAL> operator + ( Poly<REAL> const & a, REAL b ); \
  template Poly<REAL> operator + ( REAL a, Poly<REAL> const & b ); \
  template Poly<REAL> operator - ( Poly<REAL> const & a, Poly<REAL> const & b ); \
  template Poly<REAL> operator - ( Poly<REAL> const & a, REAL b ); \
  template Poly<REAL> operator - ( REAL a, Poly<REAL> const & b ); \
  template Poly<REAL> operator * ( Poly<REAL> const & a, Poly<REAL> const & b ); \
  template Poly<REAL> operator * ( REAL a, Poly<REAL> const & b ); \
  template Poly<REAL> operator * ( Poly<REAL> const & a, REAL b ); \
  template void divide( Poly<REAL> const & p, Poly<REAL> const & q, Poly<REAL> & M, Poly<REAL> & R ); \
  template void GCD( Poly<REAL> const & p, Poly<REAL> const & q, Poly<REAL> & g, REAL epsi ); \
}
#endif

namespace Utils {

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real> &
  Poly<Real>::set_scalar( Real a ) {
    this->resize(1);
    this->coeffRef(0) = a;
    m_order = 1;
    return *this;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real> &
  Poly<Real>::set_monomial( Real a ) {
    this->resize(2);
    this->coeffRef(0) = a;
    this->coeffRef(1) = 1;
    m_order           = 2;
    return *this;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Real
  Poly<Real>::eval( Real x ) const {
    // Calcolo il polinomio usando il metodo di Horner
    Integer n = m_order-1;
    Real res = this->coeff(n);
    while ( n-- > 0 ) res = res*x+this->coeff(n);
    return res;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Real
  Poly<Real>::eval_D( Real x ) const {
    // Calcolo il polinomio usando il metodo di Horner
    Integer n = m_order-1;
    Real Dp = this->coeff(n)*n;
    while ( --n > 0 ) Dp = Dp*x+this->coeff(n)*n;
    return Dp;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  void
  Poly<Real>::eval( Real x, Real & p, Real & Dp ) const {
    // Calcolo il polinomio usando il metodo di Horner
    Integer n = m_order-1;
    p  = this->coeff(n);
    Dp = this->coeff(n)*n;
    while ( --n > 0 ) {
      p  = p*x+this->coeff(n);
      Dp = Dp*x+this->coeff(n)*n;
    }
    p = p*x+this->coeff(0);
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  void
  Poly<Real>::derivative( Poly<Real> & der ) const {
    der.resize( m_order-1 ); // nuovo polinomio contenente il risultato
    for( Integer i = 1; i < m_order; ++i )
      der.coeffRef(i-1) = i * this->coeff(i);
    der.m_order = m_order-1;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  void
  Poly<Real>::integral( Poly<Real> & itg ) const {
    itg.resize( m_order+1 ); // nuovo polinomio contenente il risultato
    itg.coeffRef(0) = 0;
    for ( Integer i = 1; i <= m_order; ++i )
      itg.coeffRef(i) = this->coeff(i-1)/i;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Real
  Poly<Real>::normalize() {
    // search max module coeff
    Real S = this->cwiseAbs().maxCoeff();
    if ( S > 0 ) this->to_eigen() /= S;
    return S;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  void
  Poly<Real>::purge( Real epsi ) {
    if ( m_order > 0 ) {
      Real MX = this->cwiseAbs().maxCoeff();
      if ( MX < 1 ) MX = 1;
      Real EPS = epsi*MX;
      for ( Integer i = 0; i < m_order; ++i ) {
        Real & ai = this->coeffRef(i);
        if ( std::abs( ai ) <= EPS ) ai = 0;
      }
    }
    adjust_degree();
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  void
  Poly<Real>::adjust_degree() {
    while ( m_order > 0 && this->coeff(m_order-1) == 0 ) --m_order;
    this->conservativeResize( m_order );
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  typename Poly<Real>::Integer
  Poly<Real>::sign_variations() const {
    Integer sign_var  = 0;
    Integer last_sign = 0;
    for ( Integer i=0; i < m_order; ++i ) {
      Real v = this->coeff(i);
      if ( v > 0 ) {
        if ( last_sign == -1 ) ++sign_var;
        last_sign = 1;
      } else if ( v < 0 ) {
        if ( last_sign == 1 ) ++sign_var;
        last_sign = -1;
      }
    }
    return sign_var;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real> &
  Poly<Real>::operator = ( Poly<Real> const & b ) {
    this->resize( b.m_order );
    this->to_eigen().noalias() = b.to_eigen();
    m_order = b.m_order;
    return *this;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real> &
  Poly<Real>::operator += ( Poly<Real> const & b ) {
    Integer max_order = std::max( m_order, b.m_order );
    Integer min_order = std::min( m_order, b.m_order );

    // ridimensiona vettore coefficienti senza distruggere il contenuto
    this->conservativeResize( max_order );

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    this->head( min_order ).noalias() += b.head(min_order);
    Integer n_tail = b.m_order - m_order;
    if ( n_tail > 0 ) this->tail( n_tail ).noalias() = b.tail(n_tail);
    m_order = max_order;
    return *this;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real> &
  Poly<Real>::operator += ( Real b ) {
    if ( m_order > 0 ) this->coeffRef(0) += b;
    else {
      this->resize(1);
      this->coeffRef(0) = b;
      m_order = 1;
    }
    return *this;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real> &
  Poly<Real>::operator -= ( Poly<Real> const & b ) {
    Integer max_order = std::max( m_order, b.m_order );
    Integer min_order = std::min( m_order, b.m_order );

    // ridimensiona vettore coefficienti senza distruggere il contenuto
    this->conservativeResize( max_order );

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    this->head( min_order ).noalias() -= b.head(min_order);
    Integer n_tail = b.m_order - m_order;
    if ( n_tail > 0 ) this->tail( n_tail ).noalias() = -b.tail(n_tail);
    m_order = max_order;
    return *this;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real> &
  Poly<Real>::operator -= ( Real b ) {
    if ( m_order > 0 ) this->coeffRef(0) -= b;
    else {
      this->resize(1);
      this->coeffRef(0) = -b;
      m_order = 1;
    }
    return *this;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real> &
  Poly<Real>::operator *= ( Poly<Real> const & b ) {
    dvec_t a(this->to_eigen()); // fa una copia dei coefficienti del vettore
    Integer new_order = m_order + b.m_order - 1;
    this->resize( m_order + b.m_order - 1 ); // nuovo polinomio contenente il risultato
    this->setZero();
    for( Integer i=0; i<m_order; ++i )
      for( Integer j=0; j<b.m_order; ++j )
        this->coeffRef(i+j) += a.coeff(i) * b.coeff(j);
    m_order = new_order;
    return *this;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real> &
  Poly<Real>::operator *= ( Real b ) {
    this->to_eigen() *= b;
    return *this;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real>
  operator + ( Poly<Real> const & a, Poly<Real> const & b ) {
    using Integer = typename Poly<Real>::Integer;
    Integer max_order = std::max( a.order(), b.order() );
    Integer min_order = std::min( a.order(), b.order() );
    Poly<Real> sum( max_order ); // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    sum.head( min_order ).noalias() = a.head(min_order) + b.head(min_order);
    Integer n_tail = max_order - min_order;
    if ( n_tail > 0 ) {
      if ( a.order() > b.order() ) sum.tail( n_tail ).noalias() = a.tail(n_tail);
      else                         sum.tail( n_tail ).noalias() = b.tail(n_tail);
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real>
  operator + ( Poly<Real> const & a, Real b ) {
    using Integer = typename Poly<Real>::Integer;
    Integer max_order = std::max( a.order(), 1 );
    Poly<Real> sum( max_order ); // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    if ( a.order() > 0 ) {
      sum.coeffRef(0) = a.coeff(0) + b;
      if ( a.order() > 1 ) sum.tail( a.order()-1 ).noalias() = a.tail(a.order()-1);
    } else {
      sum.coeffRef(0) = b;
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real>
  operator + ( Real a, Poly<Real> const & b ) {
    using Integer = typename Poly<Real>::Integer;
    Integer max_order = std::max( b.order(), 1 );
    Poly<Real> sum( max_order ); // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    if ( b.order() > 0 ) {
      sum.coeffRef(0) = a + b.coeff(0);
      if ( b.order() > 1 ) sum.tail( b.order()-1 ).noalias() = b.tail(b.order()-1);
    } else {
      sum.coeffRef(0) = a;
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real>
  operator - ( Poly<Real> const & a, Poly<Real> const & b ) {
    using Integer = typename Poly<Real>::Integer;
    Integer max_order = std::max( a.order(), b.order() );
    Integer min_order = std::min( a.order(), b.order() );
    Poly<Real> sum( max_order ); // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    sum.head( min_order ).noalias() = a.head(min_order) - b.head(min_order);
    Integer n_tail = max_order - min_order;
    if ( n_tail > 0 ) {
      if ( a.order() > b.order() ) sum.tail( n_tail ).noalias() = a.tail(n_tail);
      else                         sum.tail( n_tail ).noalias() = -b.tail(n_tail);
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real>
  operator - ( Poly<Real> const & a, Real b ) {
    using Integer = typename Poly<Real>::Integer;
    Integer max_order = std::max( a.order(), 1 );
    Poly<Real> sum( max_order ); // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    if ( a.order() > 0 ) {
      sum.coeffRef(0) = a.coeff(0) - b;
      if ( a.order() > 1 ) sum.tail( a.order()-1 ).noalias() = a.tail(a.order()-1);
    } else {
      sum.coeffRef(0) = -b;
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real>
  operator - ( Real a, Poly<Real> const & b ) {
    using Integer = typename Poly<Real>::Integer;
    Integer max_order = std::max( b.order(), 1 );
    Poly<Real> sum( max_order ); // nuovo polinomio contenente la somma

    // somma i coefficienti fino al grado comune ad entrambi i polinomi
    if ( b.order() > 0 ) {
      sum.coeffRef(0) = a - b.coeff(0);
      if ( b.order() > 1 ) sum.tail( b.order()-1 ).noalias() = -b.tail(b.order()-1);
    } else {
      sum.coeffRef(0) = a;
    }
    return sum;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real>
  operator * ( Poly<Real> const & a, Poly<Real> const & b ) {
    using Integer = typename Poly<Real>::Integer;
    Poly<Real> prd( a.order() + b.order() - 1 ); // nuovo polinomio contenente il risultato
    for( Integer i = 0; i < a.order(); ++i )
      for( Integer j = 0; j < b.order(); ++j )
        prd.coeffRef(i+j) += a.coeff(i) * b.coeff(j);
    return prd;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real>
  operator * ( Real a, Poly<Real> const & b ) {
    Poly<Real> prd( b.order() ); // nuovo polinomio contenente il risultato
    prd.noalias() = a*b.to_eigen();
    return prd;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  Poly<Real>
  operator * ( Poly<Real> const & a, Real b ) {
    Poly<Real> prd( a.order() ); // nuovo polinomio contenente il risultato
    prd.noalias() = a.to_eigen()*b;
    return prd;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */

  //!
  //! Divide the polynomial \f$ p(x) \f$ by \f$ q(x) \f$
  //! with remainder, i.e. \f$ p(x) = s(x)q(x)+r(x) \f$.
  //!
  template <typename Real>
  void
  divide(
    Poly<Real> const & p,
    Poly<Real> const & q,
    Poly<Real>       & M,
    Poly<Real>       & R
  ) {

    // static Real epsi = pow(std::numeric_limits<Real>::epsilon(),Real(0.75));

    using Integer = typename Poly<Real>::Integer;

    Poly<Real> P(p), Q(q);
    //
    //  scale polynomials
    //
    //  P(x) = p(x) / scaleP, Q(x) = q(x) / scaleQ
    //
    Real scaleP = P.normalize();
    Real scaleQ = Q.normalize();

    //
    // P(x) = Q(x) * M(x) + R(x)
    //
    R = P;
    Real lcQ = Q.leading_coeff();
    Integer dd = R.order() - Q.order();
    if ( dd < 0 ) {
      // P = Q +R
      M.set_scalar(1);
      R = P-Q;
    } else {
      Integer R_degree = R.degree();
      M.set_order(dd+1);

      UTILS_ASSERT0(
        !is_zero(lcQ),
        "Poly::divide(p,q,M,R), leading coefficient of q(x) is 0!"
      );

      while ( dd >= 0 && R_degree >= 0 ) {
        Real lcR = R(R_degree);
        Real bf  = lcR/lcQ;
        M.coeffRef(dd) = bf;
        R.segment(dd,Q.degree()).noalias() -= bf*Q.head(Q.degree());
        R.coeffRef(R_degree) = 0;
        --R_degree;
        --dd;
      }

      // don not purge remainder
      // this can be done externally
      // R.purge(epsi);
      R.adjust_degree();
    }

    // scale back polinomials
    //
    // P(x) = Q(x) * M(x) + R(x)
    // p(x) / scaleP = q(x) / scaleQ * M(x) + R(x)
    // p(x) = q(x) * (scaleP/scaleQ) * M(x) + scaleP*R(x)
    //
    M *= scaleP/scaleQ;
    R *= scaleP;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  //!
  //! Given \f$ p(x) \f$ and \f$ q(x) \f$ compute G.C.D, i.e.
  //! the polinomial \f$ g(x) \f$ such that \f$ q(x) | p(x) \f$
  //! and \f$ g(x) | q(x) \f$ and if another polynomial
  //! \f$ h(x) \f$ is such that \f$ h(x) | p(x) \f$
  //! and \f$ h(x) | q(x) \f$ then \f$ h(x) | g(x) \f$.
  //!
  template <typename Real>
  void
  GCD(
    Poly<Real> const & p,
    Poly<Real> const & q,
    Poly<Real>       & g,
    Real               epsi
  ) {
    if ( q.order() > 0 ) {
      Poly<Real> M, R;
      divide( p, q, M, R );
      R.purge( epsi );
      GCD( q, R, g, epsi );
    } else {
      g = p;
    }
    g.normalize();
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  void
  Sturm<Real>::build( Poly_t const & P ) {
    m_intervals.clear();
    Poly_t DP, M, R;
    P.derivative( DP );
    m_sturm.clear();
    m_sturm.reserve(P.order());
    m_sturm.emplace_back(P);  m_sturm.back().adjust_degree();
    m_sturm.emplace_back(DP); m_sturm.back().adjust_degree();
    Integer ns = 1;
    while ( true ) {
      divide( m_sturm[ns-1], m_sturm[ns], M, R );
      if ( R.order() <= 0 ) break;
      m_sturm.emplace_back(-R);
      ++ns;
    }
    // divide by GCD
    for ( Integer i = 0; i < ns; ++i ) {
      divide( m_sturm[i], m_sturm.back(), M, R );
      M.normalize();
      m_sturm[i] = M;
    }
    m_sturm.back().set_scalar(1);
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  typename Sturm<Real>::Integer
  Sturm<Real>::sign_variations( Real x, bool & on_root ) const {
    Integer sign_var  = 0;
    Integer last_sign = 0;
    Integer npoly     = Integer(m_sturm.size());
    Real    v         = m_sturm[0].eval(x);
    on_root = false;
    if      ( v > 0 ) last_sign = 1;
    else if ( v < 0 ) last_sign = -1;
    else {
      on_root   = true;
      last_sign = 0;
    }
    for ( Integer i = 1; i < npoly; ++i ) {
      v = m_sturm[i].eval(x);
      if ( v > 0 ) {
        if ( last_sign == -1 ) ++sign_var;
        last_sign = 1;
      } else if ( v < 0 ) {
        if ( last_sign == 1 ) ++sign_var;
        last_sign = -1;
      }
    }
    return sign_var;
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  typename Sturm<Real>::Integer
  Sturm<Real>::separate_roots( Real a_in, Real b_in ) {
    using std::abs;
    using std::max;

    m_intervals.clear();
    m_intervals.reserve( m_sturm.size() );

    Interval I0, I1;
    m_a = I0.a = a_in;
    m_b = I0.b = b_in;

    I0.va = sign_variations( I0.a, I0.a_on_root );
    I0.vb = sign_variations( I0.b, I0.b_on_root );

    Integer n_roots = std::abs( I0.va - I0.vb );

    if ( n_roots <= 1 ) {
      if ( n_roots == 1 && !I0.a_on_root && !I0.b_on_root ) {
        m_intervals.push_back(I0);
      }
      if ( I0.a_on_root ) {
        I1.a  = I1.b  = I0.a;
        I1.va = I1.vb = I0.va;
        I1.a_on_root = I1.b_on_root = true;
        m_intervals.push_back(I1);
      }
      if ( I0.b_on_root ) {
        I1.a  = I1.b  = I0.b;
        I1.va = I1.vb = I0.vb;
        I1.a_on_root = I1.b_on_root = true;
        m_intervals.push_back(I1);
      }
      return m_intervals.size();
    }

    // search intervals
    vector<Interval> I_stack; I_stack.clear(); I_stack.reserve(m_sturm.size());
    I_stack.push_back(I0);
    while ( I_stack.size() > 0 ) {
      I0 = I_stack.back(); I_stack.pop_back();
      // controllo se una sola radice
      n_roots = std::abs( I0.va - I0.vb );
      if ( n_roots <= 1 ) {
        if ( I0.a_on_root ) {
          I0.b         = I0.a;
          I0.vb        = I0.va;
          I0.b_on_root = true;
          m_intervals.push_back(I0);
        } else if ( I0.b_on_root ) {
          I0.a         = I0.b;
          I0.va        = I0.vb;
          I0.a_on_root = true;
          m_intervals.push_back(I0);
        } else if ( n_roots == 1 ) {
          m_intervals.push_back(I0);
        }
      } else if ( abs(I0.b-I0.a) <= 10*machine_eps<Real>()*max(Real(1),max(abs(I0.b),abs(I0.a))) ) {
        I1.a  = I1.b  = I0.a;
        I1.va = I1.vb = 0;
        I1.a_on_root = I1.b_on_root = true;
        I_stack.push_back(I1);
      } else {
        Real    c  = (I0.a+I0.b)/2;
        bool    c_on_root;
        Integer vc = sign_variations( c, c_on_root );
        // check interval [a,c]
        if ( I0.va != vc || c_on_root || I0.a_on_root ) {
          if ( c < I1.b ) { // check if it is a true reduction
            I1.a = I0.a; I1.va = I0.va; I1.a_on_root = I0.a_on_root;
            I1.b = c;    I1.vb = vc;    I1.b_on_root = c_on_root;
            I_stack.push_back(I1);
          } else if ( c_on_root ) {
            I1.a         = I1.b         = c;
            I1.a_on_root = I1.b_on_root = true;
            I1.va        = I1.vb        = 0;
            I_stack.push_back(I1);
          } else if ( I0.a_on_root ) {
            I1.a         = I1.b         = I0.a;
            I1.a_on_root = I1.b_on_root = true;
            I1.va        = I1.vb        = 0;
            I_stack.push_back(I1);
          }
        }
        // check interval [c,b]
        if ( I0.vb != vc || I0.b_on_root ) {
          if ( c > I0.a ) {
            I1.a = c;    I1.va = vc;    I1.a_on_root = c_on_root;
            I1.b = I0.b; I1.vb = I0.vb; I1.b_on_root = I0.b_on_root;
            I_stack.push_back(I1);
          } else if ( I0.b_on_root ) {
            I1.a         = I1.b         = I0.b;
            I1.a_on_root = I1.b_on_root = true;
            I1.va        = I1.vb        = 0;
            I_stack.push_back(I1);
          }
        }
      }
    }
    // sort intervals
    std::sort(
      m_intervals.begin(),
      m_intervals.end(),
      []( Interval const & Sa, Interval const & Sb ) { return Sa.a < Sb.a; }
    );
    return m_intervals.size();
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  typename Sturm<Real>::Integer
  Sturm<Real>::separate_roots() {
    // Cauchy's bounds for roots
    Real an  = m_sturm[0].leading_coeff();
    Real bnd = 1+m_sturm[0].cwiseAbs().maxCoeff()/abs(an);
    return separate_roots( -bnd, bnd );
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Real>
  void
  Sturm<Real>::refine_roots() {
    m_fun.setup( &m_sturm[0] );
    m_roots.resize( m_intervals.size() );
    Integer n = 0;
    for ( auto & I : m_intervals ) {
      Real & r = m_roots.coeffRef(n++);
      if      ( I.a_on_root ) r = I.a;
      else if ( I.b_on_root ) r = I.b;
      else {
        r = m_solver.eval( I.a, I.b, &m_fun );
        if ( !m_solver.converged() )
          fmt::print( "Warning: Sturm<Real>::refine_roots failed at interval N.{}\n", n );
      }
    }
  }

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
UTILS_POLY_INTANTIATE( double );
UTILS_POLY_INTANTIATE( float );
#endif
