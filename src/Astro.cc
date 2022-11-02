/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2007                                                      |
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
 |      Via Mesiano 77, I-38050 Trento, Italy                               |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Astro.hh"

namespace AstroLib {

  using std::abs;
  using std::sqrt;
  using std::tan;
  using std::sin;
  using std::cos;
  using std::tanh;
  using std::sinh;
  using std::cosh;
  using std::atan;
  using std::atan2;

  bool
  check_EQ_for_consistency(
    real_type p,
    real_type f,
    real_type g,
    real_type h,
    real_type k,
    real_type L0
  ) {
    real_type hk = h*h+k*k;
    real_type e  = hypot(f,g);

    bool ok = p > 0 && hk < 1;
    if ( ok && e > 1 ) {
      real_type theta0  = L0-atan2(g,f);
      real_type max_ang = 2*atan( sqrt( (e+1)/(e-1) ) );
      ok = abs(theta0) < max_ang;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Astro::Astro( string const & __name )
  : m_name(__name)
  , m_t0(0)
  , m_M0(0)
  , m_M0_theta0(0)
  , m_M0_e(0)
  , m_Mdot(0)
  , m_Mdot_p(0)
  , m_Mdot_f(0)
  , m_Mdot_g(0)
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Astro::Astro()
  : Astro("noname")
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Astro::Astro( Astro const & ast ) {
    *this = ast;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Astro::~Astro() {
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Astro &
  Astro::operator = ( Astro const & ast ) {
    this->m_name = ast.m_name;
    m_EQ        = ast.m_EQ;
    m_K         = ast.m_K;
    m_t0        = ast.m_t0; // days
    m_M0        = ast.m_M0; // Angle corresponding to time t0
    m_M0_theta0 = ast.m_M0_theta0;
    m_M0_e      = ast.m_M0_e;
    m_L0        = ast.m_L0;
    m_theta0    = ast.m_theta0;
    m_Mdot      = ast.m_Mdot;
    m_Mdot_p    = ast.m_Mdot_p;
    m_Mdot_f    = ast.m_Mdot_f;
    m_Mdot_g    = ast.m_Mdot_g;
    m_muS       = ast.m_muS;
    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Astro::check_for_consistency() const {
    return check_EQ_for_consistency( m_EQ, m_L0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::latitude_of_periapsis() const {
    real_type angle = m_K.omega+m_K.Omega;
    angle_in_range(angle);
    return angle;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::latitude_of_apoapsis() const {
    real_type angle = m_K.omega+m_K.Omega+m_pi;
    angle_in_range(angle);
    return angle;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::M0_Mdot_grad_eval() {
    real_type f    = m_EQ.f;
    real_type g    = m_EQ.g;
    real_type absa = m_K.a > 0 ? m_K.a : -m_K.a;
    real_type e    = m_K.e;

    {
      real_type t1 = sqrt(m_muS);
      real_type t2 = pow(absa, -1.5);
      real_type t3 = 1/(1-e*e);

      m_Mdot = t1 * t2;

      real_type t5 = -3 * m_Mdot * t3;

      m_Mdot_p = -1.5 * t1 * t2 * t3 / absa;
      m_Mdot_f = t5 * f;
      m_Mdot_g = t5 * g;
      if ( m_K.a < 0 ) {
        m_Mdot_p = -m_Mdot_p;
        m_Mdot_f = -m_Mdot_f;
        m_Mdot_g = -m_Mdot_g;
      }
    }

    real_type C0 = cos(m_theta0);
    real_type S0 = sin(m_theta0);

    if ( e <= 1 ) {

      /*
        E(t) - e*sin(E(t)) = M0 + (t-t0)* Mdot
        for t = t0 -->  M0 = E0 - e*sin(E0)
        E0 = 2*arctan( sqrt(1-e)/sqrt(1+e) * tan( theta0/2 ) )
        theta0 = L0 - arctan(g,f);

        M0(theta0,e) == M0(L0 - arctan(g,f),sqrt(f^2+g^2))
       */
      real_type t1 = 1-e*e;
      real_type t2 = sqrt(t1);
      real_type t3 = C0*e;
      real_type t4 = 1/(1+t3); t4 *= t4 * t2;
      m_M0_theta0 = t1 * t4;
      m_M0_e      = -S0 * (t3+2) * t4;

    } else {

      /*
        e*sinh(E(t)) - E(t) = M0 + (t-t0)* Mdot
        for t = t0 -->  M0 = e*sinh(E0) - E0
        E0 = 2*arctanh( sqrt(e-1)/sqrt(1+e) * tan( theta0/2 ) )
        theta0 = L0 - arctan(g,f);
       */

      real_type t1 = e*e - 1;
      real_type t2 = sqrt(t1);
      real_type t3 = C0*e;
      real_type t4 = 1/(1+t3); t4 *= t4 * t2;
      m_M0_theta0 = t1 * t4;
      m_M0_e      = S0 * (t3 + 2) * t4;
    }
  }

  /*
  //            _
  //   ___  ___| |_ _   _ _ __
  //  / __|/ _ \ __| | | | '_ \
  //  \__ \  __/ |_| |_| | |_) |
  //  |___/\___|\__|\__,_| .__/
  //                     |_|
  */

  string
  Astro::info() const {
    return fmt::format(
      "Astro Name: {} [mode: {}]\n"
      "-------------------------------------------\n"
      "t0  = {:<12} Mdot   = {:<12}\n"
      "M0  = {:<12} L(t0)  = {:<12}\n"
      "muS = {:<12} period = {:<12}\n"
      "Equinoctial        Keplerian\n"
      "-------------------------------------------\n"
      "p = {:<12}   e     = {}\n"
      "f = {:<12}   a     = {}\n"
      "g = {:<12}   i     = {} [degree]\n"
      "h = {:<12}   Omega = {} [degree]\n"
      "k = {:<12}   omega = {} [degree]\n"
      "-------------------------------------------\n",
      m_name, (m_EQ.retrograde?"RETROGRADE":"NORMAL"),
      fmt::format("{:.6}",m_t0),
      fmt::format("{:.6}",m_Mdot),
      fmt::format("{:.6}",m_M0),
      fmt::format("{:.6}",m_L0),
      fmt::format("{:.6}",m_muS),
      fmt::format("{:.6}",period()),
      fmt::format("{:.6}",m_EQ.p), fmt::format("{:.6}",m_K.e),
      fmt::format("{:.6}",m_EQ.f), fmt::format("{:.6}",m_K.a),
      fmt::format("{:.6}",m_EQ.g), fmt::format("{:.6}",radiants_to_degrees(m_K.i)),
      fmt::format("{:.6}",m_EQ.h), fmt::format("{:.6}",radiants_to_degrees(m_K.Omega)),
      fmt::format("{:.6}",m_EQ.k), fmt::format("{:.6}",radiants_to_degrees(m_K.omega))
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Astro const &
  Astro::setup_Keplerian(
    string const & n,
    real_type t0,
    real_type a,
    real_type e,
    real_type Omega, // sempre prima Omega "Grande"!!!
    real_type omega,
    real_type i,
    real_type M0,
    real_type muS
  ) {
    m_name = n;

    m_K.a     = a;
    m_K.e     = e;
    m_K.i     = i;
    m_K.Omega = Omega;
    m_K.omega = omega;
    m_t0      = t0;
    m_M0      = M0;
    m_muS     = muS;

    from_Keplerian_to_Equinoctial( m_K, m_EQ );

    if ( e <= 1 ) {
      real_type E_values[4];
      mean_anomaly_to_E( M0, E_values, 0 );
      m_theta0 = E_to_true_anomaly( E_values[0], e );
    } else {
      real_type H_values[4];
      mean_anomaly_to_H( M0, H_values, 0 );
      m_theta0 = H_to_true_anomaly( H_values[0], e );
    }
    m_L0 = m_theta0 + atan2(m_EQ.g,m_EQ.f);

    M0_Mdot_grad_eval();

    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Astro const &
  Astro::setup_Equinoctial(
    string const & n,
    real_type      t0,
    real_type      p,
    real_type      f,
    real_type      g,
    real_type      h,
    real_type      k,
    bool           retrograde,
    real_type      L0,
    real_type      muS
  ) {

    UTILS_ASSERT(
      p > 0,
      "Astro::setup_Equinoctial, agument p = {} must be positive!\n", p
    );

    m_EQ.p = p;
    m_EQ.f = f;
    m_EQ.g = g;
    m_EQ.h = h;
    m_EQ.k = k;

    UTILS_WARNING(
      h*h+k*k <= 1,
      "Astro::setup_Equinoctial, agument h = {} k = {} should be h^2+k^2 < 1\n",
      h, k
    );

    m_EQ.retrograde = retrograde;

    from_equinoctial_to_Keplerian( m_EQ, m_K );

    m_name = n;
    m_t0   = t0;
    m_muS  = muS;
    // Angle corresponding to time t0
    m_L0     = L0;
    m_theta0 = L0-atan2(g,f);
    m_M0     = true_anomaly_to_mean_anomaly( m_theta0, m_K.e );

    M0_Mdot_grad_eval();

    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Astro const &
  Astro::setup_Keplerian( string const & n, GenericContainer const & vars ) {
    return setup_Keplerian(
      n,
      vars("t0")    . get_number(),
      vars("a")     . get_number(),
      vars("e")     . get_number(),
      vars("Omega") . get_number(),
      vars("omega") . get_number(),
      vars("i")     . get_number(),
      vars("M0")    . get_number(),
      vars("muS")   . get_number()
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Astro const &
  Astro::setup_Equinoctial( string const & n, GenericContainer const & vars ) {
    GenericContainer const & R = vars("retrograde");
    bool retrograde = false;
    if ( R.get_type() == GC_namespace::GC_type::BOOL ) {
      retrograde = R.get_bool();
    } else {
      retrograde = R.get_as_int() != 0;
    }

    return setup_Equinoctial(
      n,
      vars("t0").get_number(),
      vars("p").get_number(),
      vars("f").get_number(),
      vars("g").get_number(),
      vars("h").get_number(),
      vars("k").get_number(),
      retrograde,
      vars("L").get_number(),
      vars("muS").get_number()
    );
  }

  bool
  Astro::setup_using_point_and_velocity(
    real_type const P[3],
    real_type const V[3],
    real_type       muS,
    real_type       t0
  ) {
    m_t0  = t0;
    m_muS = muS;
    point_and_velocity_to_Equinoctial_and_Keplerian(
      P, V, m_muS, m_EQ, m_L0, m_K, m_theta0, m_M0
    );
    M0_Mdot_grad_eval();
    return true;
  }

  /*
  //   _                      _             _
  //  (_)_ ____   ____ _ _ __(_) __ _ _ __ | |_ ___
  //  | | '_ \ \ / / _` | '__| |/ _` | '_ \| __/ __|
  //  | | | | \ V / (_| | |  | | (_| | | | | |_\__ \
  //  |_|_| |_|\_/ \__,_|_|  |_|\__,_|_| |_|\__|___/
  */
  real_type
  Astro::orbit_energy() const { return -m_muS/(2*m_K.a); }

  /*
  // Elliptic
  // solve E-e*sin(E)=M
  // solve E'-e*cos(E)E'=Mdot
  // solve E''-e*cos(E)E''=-e*sin(E)(E')^2
  //
  // Hyperbolic
  // solve e*sinh(H)-H=M
  // solve e*cosh(H)H'-H'=Mdot
  // solve e*cosh(H)H''-H''=-e*sinh(H)(H')^2
  */
  void
  Astro::eval_E(
    real_type t,
    real_type E_values[],
    integer   nderiv
  ) const {
    real_type M = m_M0 + (t-m_t0) * m_Mdot;
    if ( m_K.e <= 1 ) {
      mean_anomaly_to_E( M, E_values, nderiv );
    } else {
      mean_anomaly_to_H( M, E_values, nderiv );
    }
  }

  //
  // solve E-e*sin(E)=M
  // solve E'-e*cos(E)E'=Mdot
  // solve E''-e*cos(E)E''=-e*sin(E)(E')^2
  //
  void
  Astro::eval_L(
    real_type t,
    real_type Lvalues[],
    integer   nderiv
  ) const {
    real_type e = m_K.e;
    real_type M = m_M0 + (t-m_t0) * m_Mdot;
    if ( e <= 1 ) { // caso ellittico
      real_type E_values[4];
      mean_anomaly_to_E( M, E_values, nderiv );
      real_type E     = E_values[0];
      real_type E_D   = E_values[1];
      real_type E_DD  = E_values[2];
      real_type E_DDD = E_values[3];

      // da E calcolo theta = 2 * atan( sqrt(1+e)/sqrt(1-e)*tan(E/2))
      real_type theta = E_to_true_anomaly(E,e);

      Lvalues[0] = theta + m_K.omega; // L = theta + omega + Omega
      if ( m_EQ.retrograde ) Lvalues[0] -= m_K.Omega; // L = theta + omega + Omega
      else                   Lvalues[0] += m_K.Omega; // L = theta + omega + Omega

      if ( nderiv < 1 ) return;

      real_type cosE      = cos(E);
      real_type tmp       = 1-e*cosE;
      real_type tmp1      = sqrt(1-e*e);
      real_type dtheta_dE = tmp1/tmp;
      Lvalues[1] = dtheta_dE*E_D;
      if ( nderiv < 2 ) return;

      real_type sinE        = sin(E_D);
      real_type d2theta_d2E = -tmp1*e*sinE/power2(tmp);
      Lvalues[2] = dtheta_dE * E_DD + d2theta_d2E*power2(E_D);
      if ( nderiv < 3 ) return;

      real_type d3theta_d3E = e*tmp1*(cosE*(cosE*e+1)-2*e)/power3(tmp);
      Lvalues[3] = d3theta_d3E*E_DDD+3*d2theta_d2E*E_D*power2(E_DD)+dtheta_dE*E_DDD;

    } else {

      real_type H_values[4];
      mean_anomaly_to_H( M, H_values, nderiv );
      real_type H     = H_values[0];
      real_type H_D   = H_values[1];
      real_type H_DD  = H_values[2];
      real_type H_DDD = H_values[3];

      // da F calcolo theta = 2 * atan( sqrt(1+e)/sqrt(1-e)*tanh(F/2))
      real_type theta = H_to_true_anomaly( H, e );

      Lvalues[0] = theta + m_K.omega; // L = theta + omega + Omega
      if ( m_EQ.retrograde ) Lvalues[0] -= m_K.Omega; // L = theta + omega + Omega
      else                   Lvalues[0] += m_K.Omega; // L = theta + omega + Omega

      if ( nderiv < 1 ) return;

      real_type cosh_H    = cosh(H);
      real_type tmp       = e*cosh_H-1;
      real_type tmp1      = sqrt(e*e-1);
      real_type dtheta_dH = tmp1/tmp;
      Lvalues[1] = dtheta_dH*H_D;
      if ( nderiv < 2 ) return;

      real_type tmp2        = tmp1*e/power2(tmp);
      real_type sinh_H      = sinh(H);
      real_type d2theta_d2H = -tmp2*sinh_H;
      Lvalues[2] = dtheta_dH * H_values[2] + d2theta_d2H*power2(H_D);
      if ( nderiv < 3 ) return;

      real_type d3theta_d3H = ((cosh_H*e+1)*cosh_H-2*e)*tmp2/tmp;
      Lvalues[3] = d3theta_d3H*H_DDD+3*d2theta_d2H*H_values[1]*power2(H_DD)+dtheta_dH*H_values[3];

    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // E = Eccentric Anomaly
  // M = Mean Anomaly
  //
  // solve E-e*sin(E)=M
  // solve E'-e*cos(E)E'=Mdot
  // solve E''-e*cos(E)E''=-e*sin(E)(E')^2
  //
  // http://en.wikipedia.org/wiki/Eccentric_anomaly
  //

  void
  Astro::mean_anomaly_to_E(
    real_type M,
    real_type Evalues[],
    integer   nderiv
  ) const {

    real_type e = m_K.e;

    angle_in_range(M);

    real_type & E = Evalues[0];
    real_type dE = 0, t1 = (1-e)*(1+e);
    E = M;
    for ( integer k = 0; k < 100; ++k ) { // risolvo con Newton
      real_type SE = e*sin(E);
      real_type CE = e*cos(E);
      dE = (E-SE-M)*(1+CE)/(t1+SE*SE);
      //dE = (E-e*sin(E)-M) / (1-e*cos(E));
      if      ( dE >  0.25 ) dE =  0.25;
      else if ( dE < -0.25 ) dE = -0.25; // scalo se passi troppo grandi
      E -= dE; // per la convergenza quando e ~1
      if ( abs( dE ) < 1E-12 ) break;
    }

    // 1- e^2 cos^2 = 1-e^2 + e^2(1-cos^2) = 1-e^2+e^2*sin^2;

    UTILS_ASSERT(
      abs( dE ) < 1E-10,
      "Astro::mean_anomaly_to_E, do not converge:"
      "\nE  = {}"
      "\ndE = {}"
      "\nM  = {}"
      "\ne  = {}\n"
      "{}\n",
      E, dE, M, e, info()
    );

    if ( nderiv < 1 ) return;

    real_type & E_D = Evalues[1];
    real_type cos_E = cos(E);
    real_type bf = 1 - e*cos_E;
    E_D = m_Mdot/bf;
    if ( nderiv < 2 ) return;

    real_type & E_DD = Evalues[2];
    real_type sin_E = sin(E);
    E_DD = -e*sin_E*power2(E_D)/bf;
    if ( nderiv < 3 ) return;

    Evalues[3] = e*E_D*(cos_E*power2(E_D)+3*sin_E*E_DD)/bf;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // F = Hyperbolic Eccentric Anomaly
  // M = Mean Anomaly
  //
  // solve e*sinh(H)-H=M
  // solve e*cosh(H)H'-H'=Mdot
  // solve e*cosh(H)H''-H''=-e*sinh(H)(H')^2
  // solve e*cosh(H)H'''-H'''=-3*e*sinh(H) H' H'' -e*cosh(H)(H')^3
  //

  void
  Astro::mean_anomaly_to_H(
    real_type M,
    real_type Hvalues[],
    integer   nderiv
  ) const {

    real_type e = m_K.e;
    real_type & H = Hvalues[0];
    // M non va messo "in range"!!!

    real_type absM = M > 0 ? M : -M;
    real_type dH   = 0;

    H = 5*e-2.5 > absM ? pow( 6*absM/e, 1./3.) : log( 2*absM/e);
    for ( integer k = 0; k < 100; ++k ) { // risolvo con Newton
      dH = (e*sinh(H)-H-absM) / (e*cosh(H)-1);
      H -= dH;
      if ( abs( dH ) < 1E-12 ) break;
    }
    UTILS_ASSERT(
      abs( dH ) < 1E-10,
      "Astro::mean_anomaly_to_H, do not converge:"
      "\nH  = {}"
      "\ndE = {}"
      "\nM  = {}"
      "\ne  = {}\n"
      "{}\n",
      H, dH, M, e, info()
    );
    if ( M < 0 ) H = -H;

    if ( nderiv < 1 ) return;

    real_type & H_D  = Hvalues[1];
    real_type cosh_H = cosh(H);
    real_type bf     = e*cosh_H-1;
    H_D = m_Mdot/bf;
    if ( nderiv < 2 ) return;

    real_type & H_DD = Hvalues[2];
    real_type sinh_H = sinh(H);
    H_DD = -e*sinh_H*power2(H_D)/bf;
    if ( nderiv < 3 ) return;

    Hvalues[3] = -e*H_D*(cosh_H*power2(H_D)+3*sinh_H*H_DD)/bf;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::L_orbital( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    return L[0];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::L_orbital_D( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 1 );
    return L[1];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::L_orbital_DD( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 2 );
    return L[2];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::L_orbital( real_type _t0, real_type dt ) const {
    if ( m_K.e < 1 ) {
      real_type L0 = L_orbital( _t0 );
      real_type dL = L_orbital( _t0 + dt ) - L0;
      real_type pf = dt / period(); // frazione di periodo da aggiungere togliere
      // cerco rivoluzioni angolo
      if ( dt > 0 ) {
        while ( dL < 0 ) dL += m_2pi;
        while ( pf > 1 ) { dL += m_2pi; pf -= 1; }
      } else {
        while ( dL > 0  ) dL -= m_2pi;
        while ( pf < -1 ) { dL -= m_2pi; pf += 1; }
      }
      return L0 + dL;
    } else {
      return L_orbital( _t0 + dt );
    }
  }

  //
  // calcolo tempo t rispetto all'angolo, tbase <= t < tbase + period
  //
  real_type
  Astro::time_from_L_angle( real_type tbase, real_type L ) const {
    UTILS_ASSERT0(
      m_K.e < 1,
      "Astro::time_from_L_angle(tbase,L), can be used only for elliptic orbits\n"
    );
    real_type theta = L - m_K.omega;
    if ( m_EQ.retrograde ) theta += m_K.Omega;
    else                   theta -= m_K.Omega;

    angle_in_range(theta);
    //real_type E     = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
    //real_type E     = 2*atan2(sqrt(1-e)*sin(theta/2),sqrt(1+e)*cos(theta/2));
    //real_type M     = E-e*sin(E);

    //real_type sinE  = sqrt(1-e*e)*(sin(theta)/(1+e*cos(theta)));
    //if ( sinE > 1 ) sinE = 1; else if ( sinE < -1 ) sinE = -1;
    //real_type M     = asin(sinE)-e*sinE;
    real_type E    = 2*atan2(sqrt(1-m_K.e)*sin(theta/2),sqrt(1+m_K.e)*cos(theta/2));
    real_type M    = E-m_K.e*sin(E);
    real_type tres = m_t0 + (M-m_M0)/m_Mdot;
    real_type per  = this->period(); // m_2pi/Mdot
    UTILS_ASSERT( per > 0, "Bad period = {}\n", per );
    if ( tres > tbase ) {
      integer numPeriod = integer( (tres-tbase)/per ); // conta i periodi da tbase a tres
      tres -= numPeriod * per;
    } else {
      integer numPeriod = integer( (tbase-tres)/per ); // conta i periodi da tbase a tres
      tres += (1+numPeriod) * per;
    }
    return tres;
  }

  //
  //  solve E-e*sin(E)=M
  //  theta = 2*atan( sqrt( (1+e)/(1-e) ) * tan(E/2) )
  //
  real_type
  Astro::true_anomaly( real_type t ) const {
    real_type E;
    eval_E( t, &E, 0 );
    return eccentric_anomaly_to_true_anomaly(E,m_K.e);
  }

  /*
  //                   _ _   _
  //   _ __   ___  ___(_) |_(_) ___  _ __
  //  | '_ \ / _ \/ __| | __| |/ _ \| '_ \
  //  | |_) | (_) \__ \ | |_| | (_) | | | |
  //  | .__/ \___/|___/_|\__|_|\___/|_| |_|
  //  |_|
  */

  void
  Astro::position_by_L(
    real_type   L,
    real_type & x,
    real_type & y,
    real_type & z
  ) const {
    real_type p    = m_EQ.p;
    real_type f    = m_EQ.f;
    real_type g    = m_EQ.g;
    real_type h    = m_EQ.h;
    real_type k    = m_EQ.k;
    real_type cosL = cos(L);
    real_type sinL = sin(L);
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type hk   = h*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;

    real_type I = m_EQ.retrograde ? -1 : 1;

    x = (1+h2-k2)*X+2*I*hk*Y;
    y = I*(1-h2+k2)*Y+2*hk*X;
    z = 2*(h*Y-I*k*X);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::position(
    real_type   t,
    real_type & x,
    real_type & y,
    real_type & z
  ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    position_by_L( L[0], x, y, z );
  }

  /*
  //             _            _ _
  //  __   _____| | ___   ___(_) |_ _   _
  //  \ \ / / _ \ |/ _ \ / __| | __| | | |
  //   \ V /  __/ | (_) | (__| | |_| |_| |
  //    \_/ \___|_|\___/ \___|_|\__|\__, |
  //                                |___/
  */

  void
  Astro::velocity_by_L(
    real_type   L,
    real_type & vx,
    real_type & vy,
    real_type & vz
  ) const {
    real_type p     = m_EQ.p;
    real_type f     = m_EQ.f;
    real_type g     = m_EQ.g;
    real_type h     = m_EQ.h;
    real_type k     = m_EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type hk    = h*k;
    real_type bf    = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cos(L)+f);
    real_type sinLg = bf * (sin(L)+g);

    real_type I = m_EQ.retrograde ? -1 : 1;

    vx = I*2*hk*cosLf - (1+h2-k2)*sinLg;
    vy = I*(1-h2+k2)*cosLf - 2*hk*sinLg;
    vz = 2 * ( h*cosLf + I*k*sinLg );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::velocity(
    real_type   t,
    real_type & vx,
    real_type & vy,
    real_type & vz
  ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    velocity_by_L( L[0], vx, vy, vz );
  }

  /*
  //                      _                _   _
  //    __ _  ___ ___ ___| | ___ _ __ __ _| |_(_) ___  _ __
  //   / _` |/ __/ __/ _ \ |/ _ \ '__/ _` | __| |/ _ \| '_ \
  //  | (_| | (_| (_|  __/ |  __/ | | (_| | |_| | (_) | | | |
  //   \__,_|\___\___\___|_|\___|_|  \__,_|\__|_|\___/|_| |_|
  */

  void
  Astro::acceleration(
    real_type   t,
    real_type & ax,
    real_type & ay,
    real_type & az
  ) const {
    real_type L[4];
    eval_L( t, L, 1 );
    real_type p      = m_EQ.p;
    real_type h      = m_EQ.h;
    real_type k      = m_EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosL_D = -bf*sin(L[0])*L[1];
    real_type sinL_D =  bf*cos(L[0])*L[1];

    real_type I = m_EQ.retrograde ? -1 : 1;

    ax = I*2*hk*cosL_D - (1+h2-k2)*sinL_D;
    ay = I*(1-h2+k2)*cosL_D - 2*hk*sinL_D;
    az = 2 * ( h*cosL_D + I*k*sinL_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::jerk(
    real_type   t,
    real_type & jx,
    real_type & jy,
    real_type & jz
  ) const {
    real_type L[4];
    eval_L( t, L, 2 );
    real_type p      = m_EQ.p;
    real_type h      = m_EQ.h;
    real_type k      = m_EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosL_2 = bf * ( -sin(L[0])*L[2] - cos(L[0]) * power2(L[1]) );
    real_type sinL_2 = bf * (  cos(L[0])*L[2] - sin(L[0]) * power2(L[1]) );

    real_type I = m_EQ.retrograde ? -1 : 1;

    jx = I*2*hk*cosL_2 - (1+h2-k2)*sinL_2;
    jy = I*(1-h2+k2)*cosL_2 - 2*hk*sinL_2;
    jz = 2 * ( h*cosL_2 + I*k*sinL_2 );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::x_position( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type p    = m_EQ.p;
    real_type f    = m_EQ.f;
    real_type g    = m_EQ.g;
    real_type h    = m_EQ.h;
    real_type k    = m_EQ.k;
    real_type cosL = cos(L[0]);
    real_type sinL = sin(L[0]);
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type hk   = h*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;

    real_type I = m_EQ.retrograde ? -1 : 1;

    return (1+h2-k2)*X+2*I*hk*Y;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::y_position( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type cosL = cos(L[0]);
    real_type sinL = sin(L[0]);
    real_type p    = m_EQ.p;
    real_type f    = m_EQ.f;
    real_type g    = m_EQ.g;
    real_type h    = m_EQ.h;
    real_type k    = m_EQ.k;
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type hk   = h*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;

    real_type I = m_EQ.retrograde ? -1 : 1;

    return I*(1-h2+k2)*Y+2*hk*X;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::z_position( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type cosL = cos(L[0]);
    real_type sinL = sin(L[0]);
    real_type p    = m_EQ.p;
    real_type f    = m_EQ.f;
    real_type g    = m_EQ.g;
    real_type h    = m_EQ.h;
    real_type k    = m_EQ.k;
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;

    real_type I = m_EQ.retrograde ? -1 : 1;

    return 2*(h*Y-I*k*X);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::x_velocity( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type p     = m_EQ.p;
    real_type f     = m_EQ.f;
    real_type g     = m_EQ.g;
    real_type h     = m_EQ.h;
    real_type k     = m_EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type hk    = h*k;
    real_type bf    = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cos(L[0])+f);
    real_type sinLg = bf * (sin(L[0])+g);

    real_type I = m_EQ.retrograde ? -1 : 1;

    return I*2*hk*cosLf - (1+h2-k2)*sinLg;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::y_velocity( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type p     = m_EQ.p;
    real_type f     = m_EQ.f;
    real_type g     = m_EQ.g;
    real_type h     = m_EQ.h;
    real_type k     = m_EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type hk    = h*k;
    real_type bf    = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cos(L[0])+f);
    real_type sinLg = bf * (sin(L[0])+g);

    real_type I = m_EQ.retrograde ? -1 : 1;

    return I*(1-h2+k2)*cosLf - 2*hk*sinLg;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::z_velocity( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type p     = m_EQ.p;
    real_type f     = m_EQ.f;
    real_type g     = m_EQ.g;
    real_type h     = m_EQ.h;
    real_type k     = m_EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type bf    = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cos(L[0])+f);
    real_type sinLg = bf * (sin(L[0])+g);

    real_type I = m_EQ.retrograde ? -1 : 1;

    return 2 * ( h*cosLf + I*k*sinLg );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::x_acceleration( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 1 );
    real_type p      = m_EQ.p;
    real_type h      = m_EQ.h;
    real_type k      = m_EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosL_D = -bf*sin(L[0])*L[1];
    real_type sinL_D =  bf*cos(L[0])*L[1];

    real_type I = m_EQ.retrograde ? -1 : 1;

    return I*2*hk*cosL_D - (1+h2-k2)*sinL_D;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::y_acceleration( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 1 );
    real_type p      = m_EQ.p;
    real_type h      = m_EQ.h;
    real_type k      = m_EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosL_D = -bf*sin(L[0])*L[1];
    real_type sinL_D =  bf*cos(L[0])*L[1];

    real_type I = m_EQ.retrograde ? -1 : 1;

    return I*(1-h2+k2)*cosL_D - 2*hk*sinL_D;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::z_acceleration( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 1 );
    real_type p      = m_EQ.p;
    real_type h      = m_EQ.h;
    real_type k      = m_EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type bf     = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosL_D = -bf*sin(L[0])*L[1];
    real_type sinL_D =  bf*cos(L[0])*L[1];

    real_type I = m_EQ.retrograde ? -1 : 1;

    return 2 * ( h*cosL_D + I*k*sinL_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::x_jerk( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 2 );
    real_type p      = m_EQ.p;
    real_type h      = m_EQ.h;
    real_type k      = m_EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosL   = cos(L[0]);
    real_type sinL   = sin(L[0]);
    real_type cosL_2 = bf * ( -sinL*L[2] - cosL * power2(L[1]) );
    real_type sinL_2 = bf * (  cosL*L[2] - sinL * power2(L[1]) );

    real_type I = m_EQ.retrograde ? -1 : 1;

    return I*2*hk*cosL_2 - (1+h2-k2)*sinL_2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::y_jerk( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 2 );
    real_type p      = m_EQ.p;
    real_type h      = m_EQ.h;
    real_type k      = m_EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosL   = cos(L[0]);
    real_type sinL   = sin(L[0]);
    real_type cosL_2 = bf * ( -sinL*L[2] - cosL * power2(L[1]) );
    real_type sinL_2 = bf * (  cosL*L[2] - sinL * power2(L[1]) );

    real_type I = m_EQ.retrograde ? -1 : 1;

    return I*(1-h2+k2)*cosL_2 - 2*hk*sinL_2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::z_jerk( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 2 );
    real_type p      = m_EQ.p;
    real_type h      = m_EQ.h;
    real_type k      = m_EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type bf     = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosL   = cos(L[0]);
    real_type sinL   = sin(L[0]);
    real_type cosL_2 = bf * ( -sinL*L[2] - cosL * power2(L[1]) );
    real_type sinL_2 = bf * (  cosL*L[2] - sinL * power2(L[1]) );

    real_type I = m_EQ.retrograde ? -1 : 1;

    return 2 * ( h*cosL_2 + I*k*sinL_2 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::radius_by_L( real_type L ) const {
    real_type p = m_EQ.p;
    real_type f = m_EQ.f;
    real_type g = m_EQ.g;
    return p/(1+f*cos(L)+g*sin(L));
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::radius_by_L_D( real_type L ) const {
    real_type p    = m_EQ.p;
    real_type f    = m_EQ.f;
    real_type g    = m_EQ.g;
    real_type sinL = sin(L);
    real_type cosL = cos(L);
    return p*(f*sinL-g*cosL)/power2(1+f*cosL+g*sinL);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::radius_by_L_DD( real_type L ) const {
    real_type p     = m_EQ.p;
    real_type f     = m_EQ.f;
    real_type g     = m_EQ.g;
    real_type cosL  = cos(L);
    real_type sinL  = sin(L);
    real_type fcosL = f*cosL;
    real_type gsinL = g*sinL;
    real_type tmp   = 1.0+fcosL+gsinL;
    return p*(2*power2(g*cosL-f*sinL)/tmp+(fcosL+gsinL) )/power2(tmp);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // normale al piano dell'ellisse
  void
  Astro::normal( real_type N[3] ) const {
    real_type h  = m_EQ.h;
    real_type k  = m_EQ.k;
    real_type k2 = k*k;
    real_type h2 = h*h;
    real_type tmp = 1/(k2+h2+1);
    N[0] =  2*tmp*k;
    N[1] = -2*tmp*h;
    N[2] =  tmp*(1-k2-h2);
    if ( m_EQ.retrograde ) N[2] = -N[2];
  }

  /*
  // terna solidare al satellite
  // colonna 0 = vettore dal satellite verso il centro (fuoco ellisse)
  // colonna 1 = vettore normale a n0 e n2 in derezione opposta alla velocità
  // colonna 2 = normale al piano
  // NB: la matrice inversa è la trasposta
  */
  void
  Astro::local_frame_by_L( real_type L, real_type M[3][3] ) const {
    real_type p    = m_EQ.p;
    real_type f    = m_EQ.f;
    real_type g    = m_EQ.g;
    real_type h    = m_EQ.h;
    real_type k    = m_EQ.k;
    real_type cosL = cos(L);
    real_type sinL = sin(L);
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type hk   = h*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;

    real_type I = m_EQ.retrograde ? -1 : 1;

    real_type x = (1+h2-k2)*X+2*I*hk*Y;
    real_type y = I*(1-h2+k2)*Y+2*hk*X;
    real_type z = 2*(h*Y-I*k*X);

    real_type bf1   = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosLf = bf1 * (cos(L)+f);
    real_type sinLg = bf1 * (sin(L)+g);

    real_type vx = I*2*hk*cosLf - (1+h2-k2)*sinLg;
    real_type vy = I*(1-h2+k2)*cosLf - 2*hk*sinLg;
    real_type vz = 2 * ( h*cosLf + I*k*sinLg );

    real_type r = sqrt(x*x+y*y+z*z);
    x /= r;
    y /= r;
    z /= r;

    real_type r_cross_v_x = y*vz-z*vy;
    real_type r_cross_v_y = z*vx-x*vz;
    real_type r_cross_v_z = x*vy-y*vx;
    real_type r_cross_v_z_2 = sqrt( r_cross_v_x*r_cross_v_x +
                                    r_cross_v_y*r_cross_v_y +
                                    r_cross_v_z*r_cross_v_z );

    r_cross_v_x /= r_cross_v_z_2;
    r_cross_v_y /= r_cross_v_z_2;
    r_cross_v_z /= r_cross_v_z_2;

    real_type nx = r_cross_v_y*z-r_cross_v_z*y;
    real_type ny = r_cross_v_z*x-r_cross_v_x*z;
    real_type nz = r_cross_v_x*y-r_cross_v_y*x;
    real_type n  = sqrt(nx*nx+ny*ny+nz*nz);

    M[0][0] = x;
    M[0][1] = y;
    M[0][2] = z;

    M[1][0] = nx/n;
    M[1][1] = ny/n;
    M[1][2] = nz/n;

    M[2][0] = r_cross_v_x;
    M[2][1] = r_cross_v_y;
    M[2][2] = r_cross_v_z;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::local_frame( real_type t, real_type M[3][3] ) const {
    real_type L = L_orbital(t);
    local_frame_by_L( L, M );
  }

  /*
  // terna del piano dell'ellisse
  // colonna 0 = asse principale ellisse (x)
  // colonna 1 = asse normale all'asse principale (y)
  // colonna 2 = normale al piano ellisse
  */
  void
  Astro::ellipse_frame( real_type M[3][3] ) const {
    real_type h = m_EQ.h;
    real_type k = m_EQ.k;

    real_type h2 = h*h;
    real_type k2 = k*k;
    real_type hk = h*k;
    real_type bf = 1+h2+k2;

    real_type I = m_EQ.retrograde ? -1 : 1;

    M[0][0] = ( 1-k2+h2 ) / bf;
    M[0][1] = ( 2*hk ) / bf;
    M[0][2] = ( -2*I*k ) / bf;

    M[1][0] = ( 2*I*hk ) / bf;
    M[1][1] = ( (1+k2-h2)*I ) / bf;
    M[1][2] = ( 2*h ) / bf;

    M[2][0] = ( 2*k ) / bf;
    M[2][1] = ( -2*h ) / bf;
    M[2][2] = ( I*(1-h2-k2) ) / bf;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::make_retrograde() {
    if ( m_EQ.retrograde ) return;
    real_type dO = m_K.omega-m_K.Omega;
    m_EQ.f = m_K.e*cos(dO);
    m_EQ.g = m_K.e*sin(dO);
    real_type tg = tan(m_K.i/2);
    m_EQ.h = cos(m_K.Omega)/tg;
    m_EQ.k = sin(m_K.Omega)/tg;
    m_EQ.retrograde = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::make_not_retrograde() {
    if ( !m_EQ.retrograde ) return;
    m_EQ.f = m_K.e*cos(m_K.omega+m_K.Omega);
    m_EQ.g = m_K.e*sin(m_K.omega+m_K.Omega);
    real_type tg = tan(m_K.i/2);
    m_EQ.h = tg*cos(m_K.Omega);
    m_EQ.k = tg*sin(m_K.Omega);
    m_EQ.retrograde = false;
  }

  /*
  //      _                 _     _
  //     | | __ _  ___ ___ | |__ (_) __ _ _ __  ___
  //  _  | |/ _` |/ __/ _ \| '_ \| |/ _` | '_ \/ __|
  // | |_| | (_| | (_| (_) | |_) | | (_| | | | \__ \
  //  \___/ \__,_|\___\___/|_.__/|_|\__,_|_| |_|___/
  */

  void
  Astro::theta0_EQ_gradient( real_type grad[6] ) const {
    real_type f  = m_EQ.f;
    real_type g  = m_EQ.g;
    real_type e  = m_K.e;
    real_type e2 = e*e;

    // theta0 = L0 - arctan(g,f);

    grad[0] = 0;
    if ( e2 > 0 ) {
      grad[1] = g/e2;
      grad[2] = -f/e2;
    } else {
      grad[1] = 0;
      grad[2] = 0;
    }
    grad[3] = 0;
    grad[4] = 0;
    grad[5] = 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::M0_EQ_gradient( real_type grad[6] ) const {
    real_type f = m_EQ.f;
    real_type g = m_EQ.g;
    real_type e = m_K.e;

    /*\

      E(t) - e*sin(E(t))  = M0 + (t-t0)* Mdot
      e*sinh(H(t)) - H(t) = M0 + (t-t0)* Mdot

      t = t0

      E0 - e*sin(E0)  = M0
      e*sinh(H0) - H0 = M0

      E0 = 2*arctan( sqrt(1-e)/sqrt(1+e) * tan( theta0/2 ) )
      H0 = 2*arctanh( sqrt(e-1)/sqrt(1+e) * tan( theta0/2 ) )

      theta0 = L0 - arctan(g,f);

      M0 = M0( theta0, e ) = M0( L0 - arctan(g,f), sqrt( f^2 + g^2 ) )

    \*/

    real_type theta0_f=0, theta0_g=0, e_f=0, e_g=0;
    if ( e > 0 ) {
      // do not work for e == 0
      e_f      = f/e;
      e_g      = g/e;
      theta0_f = e_g/e;
      theta0_g = -e_f/e;
    }

    grad[0] = 0;
    grad[1] = m_M0_theta0*theta0_f+m_M0_e*e_f;
    grad[2] = m_M0_theta0*theta0_g+m_M0_e*e_g;
    grad[3] = 0;
    grad[4] = 0;
    grad[5] = m_M0_theta0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::E0_EQ_gradient( real_type grad[6] ) const {
    real_type E0_values[4];
    mean_anomaly_to_E( m_M0, E0_values, 0 );
    real_type E0 = E0_values[0];

    //
    // M = E-e*sin(E) = E - sqrt(f^2+g^2) * E
    //
    // dM0      dMdot                            dE            d sqrt(f^2+g^2)
    // --- + DT ----- = (1-sqrt(f^2+g^2)*cos(E)) --- -  sin(E) ---------------
    // d{}       d{}                             d{}                 d{}
    //

    //real_type tmp   = 1-e*cos(E);
    //real_type sinEe = sin(E)/e;

    real_type grad_theta0[6];
    theta0_EQ_gradient( grad_theta0 );

    real_type e  = m_K.e;
    real_type ep = sqrt(1+e);
    real_type em = sqrt(1-e);
    real_type t1 = (1+cos(E0_values[0]))/ep;
    real_type t2 = (1+cos(m_theta0))/em;
    real_type t3 = t1/t2;

    grad[0] = t3*grad_theta0[0];
    grad[1] = t3*grad_theta0[1];
    grad[2] = t3*grad_theta0[2];
    grad[3] = t3*grad_theta0[3];
    grad[4] = t3*grad_theta0[4];
    grad[5] = t3*grad_theta0[5];

    real_type f = m_EQ.f;
    real_type g = m_EQ.g;
    if ( e > 0 ) {
      real_type t4 = (t1/(2*e))*(tan(m_theta0/2)/em+tan(E0_values[0]/2)/ep);
      grad[1] -= f*t4;
      grad[2] -= g*t4;
    }
    return E0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::H0_EQ_gradient( real_type grad[6] ) const {
    real_type H_values[4];
    mean_anomaly_to_H( m_M0, H_values, 0 );
    real_type H0 = H_values[0];

    //
    // M = e*sinh(H) - H = sqrt(f^2+g^2)*sinh(H) - H
    //
    // dM0      dMdot                             dH             d sqrt(f^2+g^2)
    // --- + DT ----- = (sqrt(f^2+g^2)*cosh(H)-1) --- -  sinh(H) ---------------
    // d{}       d{}                              d{}                 d{}
    //

    real_type e      = m_K.e;
    real_type tmp    = e*cosh(H0)-1;
    real_type sinhHe = sinh(H0)/e;
    real_type f      = m_EQ.f;
    real_type g      = m_EQ.g;

    real_type grad_M0[6];
    M0_EQ_gradient( grad_M0 );

    grad[0] = grad_M0[0]/tmp;
    grad[1] = ( grad_M0[1] - sinhHe*f )/tmp;
    grad[2] = ( grad_M0[2] - sinhHe*g )/tmp;
    grad[3] = grad_M0[3]/tmp;
    grad[4] = grad_M0[4]/tmp;
    grad[5] = grad_M0[5]/tmp;

    return H0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::E_EQ_gradient( real_type t, real_type grad[6] ) const {
    real_type DT = t-m_t0;
    real_type M  = m_M0 + DT * m_Mdot;
    real_type E_values[4];
    mean_anomaly_to_E( M, E_values, 0 );
    real_type E = E_values[0];

    //
    // M = E-e*sin(E) = E - sqrt(f^2+g^2) * E
    //
    // dM0      dMdot                            dE            d sqrt(f^2+g^2)
    // --- + DT ----- = (1-sqrt(f^2+g^2)*cos(E)) --- -  sin(E) ---------------
    // d{}       d{}                             d{}                 d{}
    //

    real_type e     = m_K.e;
    real_type tmp   = 1-e*cos(E);
    real_type sinEe = sin(E)/e;
    real_type f     = m_EQ.f;
    real_type g     = m_EQ.g;

    real_type grad_M0[6];
    M0_EQ_gradient( grad_M0 );

    grad[0] = ( grad_M0[0] + DT*m_Mdot_p           )/tmp;
    grad[1] = ( grad_M0[1] + DT*m_Mdot_f + sinEe*f )/tmp;
    grad[2] = ( grad_M0[2] + DT*m_Mdot_g + sinEe*g )/tmp;
    grad[3] = grad_M0[3]/tmp;
    grad[4] = grad_M0[4]/tmp;
    grad[5] = grad_M0[5]/tmp;

    return E;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::H_EQ_gradient( real_type t, real_type grad[6] ) const {
    real_type DT = t-m_t0;
    real_type M  = m_M0 + DT * m_Mdot;
    real_type H_values[4];
    mean_anomaly_to_H( M, H_values, 0 );
    real_type H = H_values[0];

    //
    // |M| = e*sinh(H) - H = sqrt(f^2+g^2)e*sinh(H) - H
    //
    // dM0      dMdot                             dH             d sqrt(f^2+g^2)
    // --- + DT ----- = (sqrt(f^2+g^2)*cosh(H)-1) --- -  sinh(H) ---------------
    // d{}       d{}                              d{}                 d{}
    //

    real_type e      = m_K.e;
    real_type tmp    = e*cosh(H)-1;
    real_type sinhHe = sinh(H)/e;
    real_type f      = m_EQ.f;
    real_type g      = m_EQ.g;

    real_type grad_M0[6];
    M0_EQ_gradient( grad_M0 );

    real_type signM = M < 0 ? -1 : 1;

    grad[0] = ( signM*(grad_M0[0] + DT*m_Mdot_p)            )/tmp;
    grad[1] = ( signM*(grad_M0[1] + DT*m_Mdot_f) - sinhHe*f )/tmp;
    grad[2] = ( signM*(grad_M0[2] + DT*m_Mdot_g) - sinhHe*g )/tmp;
    grad[3] = signM*grad_M0[3]/tmp;
    grad[4] = signM*grad_M0[4]/tmp;
    grad[5] = signM*grad_M0[5]/tmp;

    return H;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Astro::L_orbital_EQ_gradient( real_type t, real_type grad[6] ) const {

    //real_type p = m_EQ.p;
    real_type p  = m_EQ.p;
    real_type f  = m_EQ.f;
    real_type g  = m_EQ.g;
    real_type h  = m_EQ.h;
    real_type k  = m_EQ.k;
    real_type L0 = m_L0;
    bool retrograde = m_EQ.retrograde;

    Astro P,M;

    real_type DT = t-m_t0;

    real_type dp = 1e-6*p;
    P.setup_Equinoctial("plus",m_t0,p+dp,f,g,h,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p-dp,f,g,h,k,retrograde,L0,m_muS);
    grad[0] = (P.L_orbital(m_t0,DT)-M.L_orbital(m_t0,DT))/(2*dp);

    real_type df = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f+df,g,h,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f-df,g,h,k,retrograde,L0,m_muS);
    grad[1] = (P.L_orbital(m_t0,DT)-M.L_orbital(m_t0,DT))/(2*df);

    real_type dg = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f,g+dg,h,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g-dg,h,k,retrograde,L0,m_muS);
    grad[2] = (P.L_orbital(m_t0,DT)-M.L_orbital(m_t0,DT))/(2*dg);

    real_type dh = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f,g,h+dh,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g,h-dh,k,retrograde,L0,m_muS);
    grad[3] = (P.L_orbital(m_t0,DT)-M.L_orbital(m_t0,DT))/(2*dh);

    real_type dk = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f,g,h,k+dk,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g,h,k-dk,retrograde,L0,m_muS);
    grad[4] = (P.L_orbital(m_t0,DT)-M.L_orbital(m_t0,DT))/(2*dk);

    real_type dL = 1e-5;
    P.setup_Equinoctial("plus",m_t0,p,f,g,h,k,retrograde,L0+dL,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g,h,k,retrograde,L0-dL,m_muS);
    grad[5] = (P.L_orbital(m_t0,DT)-M.L_orbital(m_t0,DT))/(2*dL);

    /*\

      E(t) - e*sin(E(t))  = M0( theta0, e ) + (t-t0)* Mdot( p, f, g )
      e*sinh(H(t)) - H(t) = M0( theta0, e ) + (t-t0)* Mdot( p, f, g )

      E = 2*arctan( sqrt(1-e)/sqrt(1+e) * tan( theta/2 ) )
      H = 2*arctanh( sqrt(e-1)/sqrt(1+e) * tan( theta/2 ) )

      L  = theta  + arctan(g,f);
      L0 = theta0 + arctan(g,f);

      theta(E,e) = 2*arctan( sqrt(1+e)/sqrt(1-e)*tan(E/2) )
      theta(H,e) = 2*arctan( sqrt(1+e)/sqrt(e-1)*tanh(H/2) )

      theta0(L0,f,g) = L0 - arctan(g,f)
      e(f,g)         = sqrt(f^2+g^2)

      E(p,f,g,L0) - e*sin(E(p,f,g,L0))  = M0( theta0, e ) + (t-t0)* Mdot( p, f, g )
      e*sinh(H(p,f,g,L0)) - H(p,f,g,L0) = M0( theta0, e ) + (t-t0)* Mdot( p, f, g )

      L_p  = theta_E * E_p
      L_f  = theta_E * E_f + theta_e * e_f
      L_g  = theta_E * E_g + theta_e * e_g
      L_h  = theta_E * 0
      L_k  = theta_E * 0
      L_L0 = theta_E * E_L0

      L_p  = theta_H * H_p
      L_f  = theta_H * H_f + theta_e * e_f
      L_g  = theta_H * H_g + theta_e * e_g
      L_h  = theta_H * 0
      L_k  = theta_H * 0
      L_L0 = theta_H * H_L0

    \*/



    /*\

      E0(p,f,g,L0) - e*sin(E0(p,f,g,L0)) = M0( theta0, e )

      L0 = theta0 + arctan(g,f);

      theta(E0(p,f,g,L0),e) = 2*arctan( sqrt(1+e)/sqrt(1-e)*tan(E0/2) )
                            = L0 - arctan(g,f)

      theta_E * E0_p = 0
      theta_E * E0_f + theta_e * e_f

      theta0(L0,f,g) = L0 - arctan(g,f)
      e(f,g)         = sqrt(f^2+g^2)

      E(p,f,g,L0) - e*sin(E(p,f,g,L0))  = M0( theta0, e ) + (t-t0)* Mdot( p, f, g )
      e*sinh(H(p,f,g,L0)) - H(p,f,g,L0) = M0( theta0, e ) + (t-t0)* Mdot( p, f, g )

      L_p  = theta_E * E_p
      L_f  = theta_E * E_f + theta_e * e_f
      L_g  = theta_E * E_g + theta_e * e_g
      L_h  = theta_E * 0
      L_k  = theta_E * 0
      L_L0 = theta_E * E_L0

      L_p  = theta_H * H_p
      L_f  = theta_H * H_f + theta_e * e_f
      L_g  = theta_H * H_g + theta_e * e_g
      L_h  = theta_H * 0
      L_k  = theta_H * 0
      L_L0 = theta_H * H_L0

    \*/

    return L_orbital(m_t0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::radius_EQ_gradient( real_type t, real_type grad[6] ) const {
    real_type p = m_EQ.p;
    real_type f = m_EQ.f;
    real_type g = m_EQ.g;
    real_type L_grad[6];
    real_type L = L_orbital_EQ_gradient( t, L_grad );

    real_type t1 = cos(L);
    real_type t2 = sin(L);
    real_type t3 = 1/(f * t1 + g * t2 + 1);
    real_type t4 = t3*t3;
    real_type t5 = p * t4;

    grad[0] = t3;
    grad[1] = -t5 * t1;
    grad[2] = -t5 * t2;
    grad[3] = 0;
    grad[4] = 0;
    grad[5] = p * (f * t2 - g * t1) * t4;

    real_type tmp = grad[5];
    grad[0] += tmp * L_grad[0];
    grad[1] += tmp * L_grad[1];
    grad[2] += tmp * L_grad[2];
    grad[3] += tmp * L_grad[3];
    grad[4] += tmp * L_grad[4];
    grad[5]  = tmp * L_grad[5];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::absolute_velocity_EQ_gradient( real_type t, real_type grad[6] ) const {
    real_type p = m_EQ.p;
    real_type f = m_EQ.f;
    real_type g = m_EQ.g;
    real_type L_grad[6];
    real_type L = L_orbital_EQ_gradient( t, L_grad );

    real_type t1 = sqrt(m_muS);
    real_type t2 = pow(p, -1.5);
    real_type t3 = sin(L);
    real_type t4 = cos(L);
    real_type t5 = f * f + g * g + 2 * (f * t4 + g * t3) + 1;
    real_type t6 = 1/sqrt(t5);
    real_type t7 = p * t2;

    grad[0] = -t1 * t2 * t5 * t6 / 2;
    grad[1] = t1 * (f + t4) * t7 * t6;
    grad[2] = t1 * (g + t3) * t7 * t6;
    grad[3] = 0;
    grad[4] = 0;
    grad[5] = -t1 * (f * t3 - g * t4) * t7 * t6;

    real_type tmp = grad[5];
    grad[0] += tmp * L_grad[0];
    grad[1] += tmp * L_grad[1];
    grad[2] += tmp * L_grad[2];
    grad[3] += tmp * L_grad[3];
    grad[4] += tmp * L_grad[4];
    grad[5]  = tmp * L_grad[5];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::position0_EQ_jacobian( real_type JP[3][6], real_type L0 ) const {
    real_type p = m_EQ.p;
    real_type f = m_EQ.f;
    real_type g = m_EQ.g;
    real_type h = m_EQ.h;
    real_type k = m_EQ.k;
    real_type I = m_EQ.retrograde ? -1 : 1;

    real_type t1  = h * h;
    real_type t2  = k * k;
    real_type t3  = 1 + t1 - t2;
    real_type t4  = cos(L0);
    real_type t5  = sin(L0);
    real_type t6  = t5 * h;
    real_type t7  = 2 * t6;
    real_type t8  = t7 * I;
    real_type t9  = t3 * t4;
    real_type t10 = t8 * k + t9;
    real_type t11 = f * t4;
    real_type t12 = g * t5;
    real_type t13 = 1 + t12 + t11;
    real_type t14 = 1 + t1 + t2;
    real_type t15 = 1 - t1 + t2;
    real_type t16 = 2 * h;
    real_type t17 = t16 * k;
    real_type t18 = t17 * t4;
    real_type t19 = I * t15;
    real_type t20 = -t19 * t5 - t18;
    real_type t21 = 2 * k;
    real_type t22 = -t3 / 2;
    real_type t23 = f * h * k * I + t22 * g;
    real_type t24 = t5*t5;
              t13 = 1 / t13;
              t14 = 1 / t14;
    real_type t25 = t13*t13;
    real_type t26 = t14*t14;
    real_type t27 = 2 * p;
    real_type t28 = p * t10;
    real_type t29 = t28 * t25 * t14;
    real_type t30 = -t19 * f + t17 * g;
    real_type t31 = t20 * p * t25 * t14;
    real_type t32 = k * t4 * I - t6;
    real_type t33 = t27 * t32 * t25 * t14;

    JP[0][0] = t10 * t13 * t14;
    JP[0][1] = -t29 * t4;
    JP[0][2] = -t29 * t5;
    JP[0][3] = -t21 * t20 * p * t13 * t26;
    JP[0][4] = (2 * t6 * I * t3 - 2 * t21 * t4 * (1 + t1)) * p * t13 * t26;
    JP[0][5] = t27 * (t22 * t5 + t23 * t24 + t4 * (h * k * I + t23 * t4)) * t25 * t14;

    JP[1][0] = -t20 * t13 * t14;
    JP[1][1] = t31 * t4;
    JP[1][2] = t31 * t5;
    JP[1][3] = -t27 * (t8 * (1 + t2) - k * t15 * t4) * t13 * t26;
    JP[1][4] = t28 * t16 * t13 * t26;
    JP[1][5] = p * (t4 * (-t30 * t4 + t19) - t5 * (t30 * t5 + t17)) * t25 * t14;

    JP[2][0] = -2 * t32 * t13 * t14;
    JP[2][1] = t33 * t4;
    JP[2][2] = t33 * t5;
    JP[2][3] = -t27 * (-t15 * t5 - t18 * I) * t13 * t26;
    JP[2][4] = -t27 * (t7 * k + t9 * I) * t13 * t26;
    JP[2][5] = t27 * (h * (f * t24 + t4 * (1 + t11)) + (g * (t4*t4) + t5 * (1 + t12)) * I * k) * t25 * t14;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::position_EQ_jacobian( real_type t, real_type JP[3][6] ) const {

    real_type L_grad[6];
    real_type L = L_orbital_EQ_gradient( t, L_grad );

    position0_EQ_jacobian( JP, L );

    /*
       P(t,p,f,g,h,k,L(t,p,f,g,h,k,L0,t0))

       /      dP      | 0 \   / dP        dL         \
       | ------------ | 0 | + | -- x --------------  |
       \ d{p,f,g,h,k} | 0 /   \ dL   d{p,f,g,h,k,L0} /

     */

    for ( integer i = 0; i < 3; ++i ) {
      real_type tmp = JP[i][5];
      JP[i][0] += tmp * L_grad[0];
      JP[i][1] += tmp * L_grad[1];
      JP[i][2] += tmp * L_grad[2];
      JP[i][3] += tmp * L_grad[3];
      JP[i][4] += tmp * L_grad[4];
      JP[i][5]  = tmp * L_grad[5];
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::position_EQ_jacobian_FD( real_type t, real_type JP[3][6] ) const {

    //real_type p = m_EQ.p;
    real_type p  = m_EQ.p;
    real_type f  = m_EQ.f;
    real_type g  = m_EQ.g;
    real_type h  = m_EQ.h;
    real_type k  = m_EQ.k;
    real_type L0 = m_L0;
    bool retrograde = m_EQ.retrograde;

    Astro P,M;

    real_type xp, yp, zp, xm, ym, zm;

    real_type dp = 1e-6*p;
    P.setup_Equinoctial("plus",m_t0,p+dp,f,g,h,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p-dp,f,g,h,k,retrograde,L0,m_muS);
    P.position(t,xp,yp,zp);
    M.position(t,xm,ym,zm);
    JP[0][0] = (xp-xm)/(2*dp);
    JP[1][0] = (yp-ym)/(2*dp);
    JP[2][0] = (zp-zm)/(2*dp);

    real_type df = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f+df,g,h,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f-df,g,h,k,retrograde,L0,m_muS);
    P.position(t,xp,yp,zp);
    M.position(t,xm,ym,zm);
    JP[0][1] = (xp-xm)/(2*df);
    JP[1][1] = (yp-ym)/(2*df);
    JP[2][1] = (zp-zm)/(2*df);

    real_type dg = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f,g+dg,h,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g-dg,h,k,retrograde,L0,m_muS);
    P.position(t,xp,yp,zp);
    M.position(t,xm,ym,zm);
    JP[0][2] = (xp-xm)/(2*dg);
    JP[1][2] = (yp-ym)/(2*dg);
    JP[2][2] = (zp-zm)/(2*dg);

    real_type dh = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f,g,h+dh,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g,h-dh,k,retrograde,L0,m_muS);
    P.position(t,xp,yp,zp);
    M.position(t,xm,ym,zm);
    JP[0][3] = (xp-xm)/(2*dh);
    JP[1][3] = (yp-ym)/(2*dh);
    JP[2][3] = (zp-zm)/(2*dh);

    real_type dk = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f,g,h,k+dk,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g,h,k-dk,retrograde,L0,m_muS);
    P.position(t,xp,yp,zp);
    M.position(t,xm,ym,zm);
    JP[0][4] = (xp-xm)/(2*dk);
    JP[1][4] = (yp-ym)/(2*dk);
    JP[2][4] = (zp-zm)/(2*dk);

    real_type dL = 1e-5;
    P.setup_Equinoctial("plus",m_t0,p,f,g,h,k,retrograde,L0+dL,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g,h,k,retrograde,L0-dL,m_muS);
    P.position(t,xp,yp,zp);
    M.position(t,xm,ym,zm);
    JP[0][5] = (xp-xm)/(2*dL);
    JP[1][5] = (yp-ym)/(2*dL);
    JP[2][5] = (zp-zm)/(2*dL);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::velocity0_EQ_jacobian( real_type JV[3][6], real_type L0 ) const {
    real_type p = m_EQ.p;
    real_type f = m_EQ.f;
    real_type g = m_EQ.g;
    real_type h = m_EQ.h;
    real_type k = m_EQ.k;
    real_type I = m_EQ.retrograde ? -1 : 1;

    real_type t1  = pow(p, -1.5);
    real_type t2  = sqrt(m_muS);
    real_type t3  = h * h;
    real_type t4  = k * k;
    real_type t5  = 1 + t3 - t4;
    real_type t6  = sin(L0);
    real_type t7  = cos(L0);
    real_type t8  = f + t7;
    real_type t9  = t8 * h;
    real_type t10 = t9 * k;
    real_type t11 = 2;
    real_type t12 = -t10 * t11 * I + g * t5 + t5 * t6;
    real_type t13 = 1 + t3 + t4;
    real_type t14 = p * t1;
    real_type t15 = 1 - t3 + t4;
    real_type t16 = I * t15;
    real_type t17 = g + t6;
    real_type t18 = h * k;
    real_type t19 = t18 * t17;
              t8  = t19 -t16 * t8 / 2;
    real_type t20 = t5 * t7;
    real_type t21 = I * h;
              t18 *= t11;
    real_type t22 = t18 * I;
              t13 = 1/t13;
    real_type t23 = t13*t13;
    real_type t24 = t14 * t2;
    real_type t25 = t24 * t11;
    real_type t26 = k * I;

    JV[0][0] = t1 * t2 * t12 * t13 / 2;
    JV[0][1] = t22 * t14 * t2 * t13;
    JV[0][2] = -t24 * t5 * t13;
    JV[0][3] = -4 * t8 * t14 * t2 * t23 * k;
    JV[0][4] = 4 * t24 * (k * ((t3 + 1) * g + t6 * (t3 + 1)) + t21 * (((1 - k) * (k + 1) + t3) * f + t20) / 2) * t23;
    JV[0][5] = -t24 * (t22 * t6 + t20) * t13;

    JV[1][0] = t8 * t1 * t2 * t13;
    JV[1][1] = t24 * t16 * t13;
    JV[1][2] = -t24 * t18 * t13;
    JV[1][3] = t25 * (-t11 * t21 * ((t4 + 1) * f + t7 * (t4 + 1)) - (g * t4 + t15 * t6 + g * (1 - t3)) * k) * t23;
    JV[1][4] = -t25 * t12 * t23 * h;
    JV[1][5] = t24 * (-t16 * t6 - t18 * t7) * t13;

    JV[2][0] = -t2 * (t26 * t17 + t9) * t1 * t13;
    JV[2][1] = t25 * h * t13;
    JV[2][2] = t26 * t25 * t13;
    JV[2][3] = t24 * (t11 * (f * t15 + t15 * t7) - 4 * t19 * I) * t23;
    JV[2][4] = t25 * (I * t5 * t17 - t10 * t11) * t23;
    JV[2][5] = t25 * (-h * t6 + t26 * t7) * t13;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::velocity_EQ_jacobian( real_type t, real_type JV[3][6] ) const {

    real_type L_grad[6];
    real_type L = L_orbital_EQ_gradient( t, L_grad );

    velocity0_EQ_jacobian( JV, L );

    /*
       P(t,p,f,g,h,k,L(t,p,f,g,h,k,L0,t0))

       /      dP      | 0 \   / dP        dL         \
       | ------------ | 0 | + | -- x --------------  |
       \ d{p,f,g,h,k} | 0 /   \ dL   d{p,f,g,h,k,L0} /

     */

    for ( integer i = 0; i < 3; ++i ) {
      real_type tmp = JV[i][5];
      JV[i][0] += tmp * L_grad[0];
      JV[i][1] += tmp * L_grad[1];
      JV[i][2] += tmp * L_grad[2];
      JV[i][3] += tmp * L_grad[3];
      JV[i][4] += tmp * L_grad[4];
      JV[i][5]  = tmp * L_grad[5];
    }

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Astro::velocity_EQ_jacobian_FD( real_type t, real_type JP[3][6] ) const {

    //real_type p = m_EQ.p;
    real_type p  = m_EQ.p;
    real_type f  = m_EQ.f;
    real_type g  = m_EQ.g;
    real_type h  = m_EQ.h;
    real_type k  = m_EQ.k;
    real_type L0 = m_L0;
    bool retrograde = m_EQ.retrograde;

    Astro P,M;

    real_type xp, yp, zp, xm, ym, zm;

    real_type dp = 1e-6*p;
    P.setup_Equinoctial("plus",m_t0,p+dp,f,g,h,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p-dp,f,g,h,k,retrograde,L0,m_muS);
    P.velocity(t,xp,yp,zp);
    M.velocity(t,xm,ym,zm);
    JP[0][0] = (xp-xm)/(2*dp);
    JP[1][0] = (yp-ym)/(2*dp);
    JP[2][0] = (zp-zm)/(2*dp);

    real_type df = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f+df,g,h,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f-df,g,h,k,retrograde,L0,m_muS);
    P.velocity(t,xp,yp,zp);
    M.velocity(t,xm,ym,zm);
    JP[0][1] = (xp-xm)/(2*df);
    JP[1][1] = (yp-ym)/(2*df);
    JP[2][1] = (zp-zm)/(2*df);

    real_type dg = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f,g+dg,h,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g-dg,h,k,retrograde,L0,m_muS);
    P.velocity(t,xp,yp,zp);
    M.velocity(t,xm,ym,zm);
    JP[0][2] = (xp-xm)/(2*dg);
    JP[1][2] = (yp-ym)/(2*dg);
    JP[2][2] = (zp-zm)/(2*dg);

    real_type dh = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f,g,h+dh,k,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g,h-dh,k,retrograde,L0,m_muS);
    P.velocity(t,xp,yp,zp);
    M.velocity(t,xm,ym,zm);
    JP[0][3] = (xp-xm)/(2*dh);
    JP[1][3] = (yp-ym)/(2*dh);
    JP[2][3] = (zp-zm)/(2*dh);

    real_type dk = 1e-7;
    P.setup_Equinoctial("plus",m_t0,p,f,g,h,k+dk,retrograde,L0,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g,h,k-dk,retrograde,L0,m_muS);
    P.velocity(t,xp,yp,zp);
    M.velocity(t,xm,ym,zm);
    JP[0][4] = (xp-xm)/(2*dk);
    JP[1][4] = (yp-ym)/(2*dk);
    JP[2][4] = (zp-zm)/(2*dk);

    real_type dL = 1e-5;
    P.setup_Equinoctial("plus",m_t0,p,f,g,h,k,retrograde,L0+dL,m_muS);
    M.setup_Equinoctial("minus",m_t0,p,f,g,h,k,retrograde,L0-dL,m_muS);
    P.velocity(t,xp,yp,zp);
    M.velocity(t,xm,ym,zm);
    JP[0][5] = (xp-xm)/(2*dL);
    JP[1][5] = (yp-ym)/(2*dL);
    JP[2][5] = (zp-zm)/(2*dL);
  }

}
