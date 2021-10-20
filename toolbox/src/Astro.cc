/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2007                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                          /|                   /|                        |
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

  Astro::Astro()
  : m_name("noname")
  , m_t0(0)
  , m_M0(0)
  , m_Mdot(0)
  {}

  Astro::Astro( string const & __name )
  : m_name(__name)
  , m_t0(0)
  , m_M0(0)
  , m_Mdot(0)
  {}

  Astro::Astro( Astro const & ast ) {
    *this = ast;
  }

  Astro::~Astro() {
  }

  Astro const &
  Astro::operator = ( Astro const & ast ) {
    this->m_name = ast.m_name;
    m_EQ   = ast.m_EQ;
    m_K    = ast.m_K;
    m_t0   = ast.m_t0; // days
    m_M0   = ast.m_M0; // Angle corresponding to time t0
    m_Mdot = ast.m_Mdot;
    m_muS  = ast.m_muS;
    return *this;
  }

  void
  Astro::check_for_consistency() const {
    //ASSERT( p>0,         "Astro::check p = " << p << " must be positive!");
    //ASSERT( f*f+g*g < 1, "Astro::check f*f+g*g = " << f*f+g*g << " must be less than 1!");
    real_type h  = m_EQ.h;
    real_type k  = m_EQ.k;
    real_type e2 = h*h+k*k;
    UTILS_ASSERT(
      e2 < 1,
      "Astro::check_for_consistency h*h+k*k = {} must be less than 1!\n", e2
    );
    UTILS_ASSERT(
      m_muS > 0,
      "Astro::check_for_consistency muS = {} must be positive!\n", m_muS
    );
  }

  real_type
  Astro::latitude_of_periapsis() const {
    real_type angle = m_K.omega+m_K.Omega;
    angle_in_range(angle);
    return angle;
  }

  real_type
  Astro::latitude_of_apoapsis() const {
    real_type angle = m_K.omega+m_K.Omega+m_pi;
    angle_in_range(angle);
    return angle;
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
      "------------------------------------------\n"
      "t0  = {:<12} Mdot   = {:<12}\n"
      "M0  = {:<12} L(t0)  = {:<12}\n"
      "muS = {:<12} period = {:<12}\n"
      "Equinoctial        Keplerian\n"
      "------------------------------------------\n"
      "p = {:<12}   e     = {}\n"
      "f = {:<12}   a     = {}\n"
      "g = {:<12}   i     = {} [degree]\n"
      "h = {:<12}   Omega = {} [degree]\n"
      "k = {:<12}   omega = {} [degree]\n"
      "------------------------------------------\n",
      m_name, (m_EQ.retrograde?"RETROGRADE":"NORMAL"),
      fmt::format("{:.6}",m_t0),
      fmt::format("{:.6}",m_Mdot),
      fmt::format("{:.6}",m_M0),
      fmt::format("{:.6}",L_orbital(m_t0)),
      fmt::format("{:.6}",m_muS),
      fmt::format("{:.6}",period()),
      fmt::format("{:.6}",m_EQ.p), fmt::format("{:.6}",m_K.e), 
      fmt::format("{:.6}",m_EQ.f), fmt::format("{:.6}",m_K.a),
      fmt::format("{:.6}",m_EQ.g), fmt::format("{:.6}",radiants_to_degrees(m_K.i)),
      fmt::format("{:.6}",m_EQ.h), fmt::format("{:.6}",radiants_to_degrees(m_K.Omega)),
      fmt::format("{:.6}",m_EQ.k), fmt::format("{:.6}",radiants_to_degrees(m_K.omega))
    );
  }

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
    real_type absa = m_K.a > 0 ? m_K.a : -m_K.a;
    m_Mdot = sqrt(m_muS/absa)/absa;

    from_Keplerian_to_Equinoctial( m_K, m_EQ );
    check_for_consistency();
    return *this;
  }

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
    real_type      L,
    real_type      muS
  ) {

    m_EQ.p          = p;
    m_EQ.f          = f;
    m_EQ.g          = g;
    m_EQ.h          = h;
    m_EQ.k          = k;
    m_EQ.retrograde = retrograde;

    from_equinoctial_to_Keplerian( m_EQ, m_K );

    m_name = n;
    m_t0   = t0;
    m_muS  = muS;
    // Angle corresponding to time t0
    real_type theta = L-m_K.omega;

    if ( retrograde ) theta += m_K.Omega;
    else              theta -= m_K.Omega;

    m_M0 = true_anomaly_to_mean_anomaly( theta, m_K.e );

    real_type absa = m_K.a > 0 ? m_K.a : -m_K.a;
    m_Mdot = sqrt(m_muS/absa)/absa;

    check_for_consistency();

    return *this;
  }

  Astro const &
  Astro::setup_Keplerian( string const & n, GenericContainer & vars ) {
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

  Astro const &
  Astro::setup_Equinoctial( string const & n, GenericContainer & vars ) {
    return setup_Equinoctial(
      n,
      vars("t0")         . get_number(),
      vars("p")          . get_number(),
      vars("f")          . get_number(),
      vars("g")          . get_number(),
      vars("h")          . get_number(),
      vars("k")          . get_number(),
      vars("retrograde") . get_bool(),
      vars("L")          . get_number(),
      vars("muS")        . get_number()
    );
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
  // solve e*sinh(F)-F=M
  // solve e*cosh(F)F'-F'=Mdot
  // solve e*cosh(F)F''-F''=-e*sinh(F)(F')^2
  */
  void
  Astro::eval_E(
    real_type t,
    real_type E[],
    integer   nderiv
  ) const {
    real_type M = m_M0 + (t-m_t0) * m_Mdot;
    if ( m_K.e <= 1 ) {
      mean_anomaly_to_eccentric_anomaly_elliptic( M, m_Mdot, m_K.e, E, nderiv );
    } else {
      mean_anomaly_to_eccentric_anomaly_hyperbolic( M, m_Mdot, m_K.e, E, nderiv );
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
    real_type M = m_M0 + (t-m_t0) * m_Mdot;
    if ( m_K.e <= 1 ) { // caso ellittico
      real_type E[4];
      mean_anomaly_to_eccentric_anomaly_elliptic( M, m_Mdot, m_K.e, E, nderiv );

      // da E calcolo theta = 2 * atan( sqrt(1+e)/sqrt(1-e)*tan(E/2))
      real_type theta = eccentric_anomaly_to_true_anomaly(E[0],m_K.e);

      Lvalues[0] = theta + m_K.omega; // L = theta + omega + Omega
      if ( m_EQ.retrograde ) Lvalues[0] -= m_K.Omega; // L = theta + omega + Omega
      else                   Lvalues[0] += m_K.Omega; // L = theta + omega + Omega

      if ( nderiv < 1 ) return;

      real_type cosE      = cos(E[0]);
      real_type tmp       = 1-m_K.e*cosE;
      real_type tmp1      = sqrt(1-m_K.e*m_K.e);
      real_type dtheta_dE = tmp1/tmp;
      Lvalues[1] = dtheta_dE*E[1];
      if ( nderiv < 2 ) return;

      real_type sinE        = sin(E[0]);
      real_type d2theta_d2E = -tmp1*m_K.e*sinE/power2(tmp);
      Lvalues[2] = dtheta_dE * E[2] + d2theta_d2E*power2(E[1]);
      if ( nderiv < 3 ) return;

      real_type d3theta_d3E = m_K.e*tmp1*(cosE*(cosE*m_K.e+1)-2*m_K.e)/power3(tmp);
      Lvalues[3] = d3theta_d3E*E[3]+3*d2theta_d2E*E[1]*power2(E[2])+dtheta_dE*E[3];

    } else {

      real_type F[4];
      mean_anomaly_to_eccentric_anomaly_hyperbolic( M, m_Mdot, m_K.e, F, nderiv );

      // da F calcolo theta = 2 * atan( sqrt(1+e)/sqrt(1-e)*tanh(F/2))
      real_type theta = eccentric_anomaly_to_true_anomaly( F[0], m_K.e );

      Lvalues[0] = theta + m_K.omega; // L = theta + omega + Omega
      if ( m_EQ.retrograde ) Lvalues[0] -= m_K.Omega; // L = theta + omega + Omega
      else                   Lvalues[0] += m_K.Omega; // L = theta + omega + Omega

      if ( nderiv < 1 ) return;

      real_type cosh_F    = cosh(F[0]);
      real_type tmp       = m_K.e*cosh_F-1;
      real_type tmp1      = sqrt(m_K.e*m_K.e-1);
      real_type dtheta_dF = tmp1/tmp;
      Lvalues[1] = dtheta_dF*F[1];
      if ( nderiv < 2 ) return;

      real_type tmp2        = tmp1*m_K.e/power2(tmp);
      real_type sinh_F      = sinh(F[0]);
      real_type d2theta_d2F = -tmp2*sinh_F;
      Lvalues[2] = dtheta_dF * F[2] + d2theta_d2F*power2(F[1]);
      if ( nderiv < 3 ) return;

      real_type d3theta_d3F = ((cosh_F*m_K.e+1)*cosh_F-2*m_K.e)*tmp2/tmp;
      Lvalues[3] = d3theta_d3F*F[3]+3*d2theta_d2F*F[1]*power2(F[2])+dtheta_dF*F[3];

    }
  }

  real_type
  Astro::L_orbital( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    return L[0];
  }

  real_type
  Astro::L_orbital_D( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 1 );
    return L[1];
  }

  real_type
  Astro::L_orbital_DD( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 2 );
    return L[2];
  }

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
      return this->L_orbital( _t0 + dt );
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

  void
  Astro::position_by_L_jacobian_EQ( real_type L, real_type JP[3][6] ) const {
    real_type p = m_EQ.p;
    real_type f = m_EQ.f;
    real_type g = m_EQ.g;
    real_type h = m_EQ.h;
    real_type k = m_EQ.k;
    real_type I = m_EQ.retrograde ? -1 : 1;

    real_type t1 = h * h;
    real_type t2 = k * k;
    real_type t3 = 1 + t1 - t2;
    real_type t4 = cos(L);
    real_type t5 = sin(L);
    real_type t6 = t5 * h;
    real_type t7 = 2 * t6;
    real_type t8 = t7 * I;
    real_type t9 = t3 * t4;
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

  void
  Astro::ray_by_L_gradient( real_type L, real_type grad[6] ) const {
    real_type p = m_EQ.p;
    real_type f = m_EQ.f;
    real_type g = m_EQ.g;

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

  void
  Astro::velocity_by_L_jacobian_EQ( real_type L, real_type JV[3][6] ) const {
    real_type p = m_EQ.p;
    real_type f = m_EQ.f;
    real_type g = m_EQ.g;
    real_type h = m_EQ.h;
    real_type k = m_EQ.k;
    real_type I = m_EQ.retrograde ? -1 : 1;

    real_type t1 = pow(p, -1.5);
    real_type t2 = sqrt(m_muS);
    real_type t3 = h * h;
    real_type t4 = k * k;
    real_type t5 = 1 + t3 - t4;
    real_type t6 = sin(L);
    real_type t7 = cos(L);
    real_type t8 = f + t7;
    real_type t9 = t8 * h;
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

  void
  Astro::absolute_velocity_by_angle_gradient( real_type L, real_type grad[6] ) const {
    real_type p = m_EQ.p;
    real_type f = m_EQ.f;
    real_type g = m_EQ.g;

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

  /////////////////////////////////////////////////////////////////////////////

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

  //////////////////////////////////////////////////////////

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

  //////////////////////////////////////////////////////////

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

  //////////////////////////////////////////////////////////

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

  //////////////////////////////////////////////////////////

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

  //////////////////////////////////////////////////////////

  real_type
  Astro::ray_by_L( real_type L ) const {
    real_type p = m_EQ.p;
    real_type f = m_EQ.f;
    real_type g = m_EQ.g;
    return p/(1+f*cos(L)+g*sin(L));
  }

  real_type
  Astro::ray_by_L_D( real_type L ) const {
    real_type p    = m_EQ.p;
    real_type f    = m_EQ.f;
    real_type g    = m_EQ.g;
    real_type sinL = sin(L);
    real_type cosL = cos(L);
    return p*(f*sinL-g*cosL)/power2(1+f*cosL+g*sinL);
  }

  real_type
  Astro::ray_by_L_DD( real_type L ) const {
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

  //////////////////////////////////////////////////////////

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
}
