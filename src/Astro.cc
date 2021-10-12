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

  Astro::Astro()
  : m_name("noname")
  , t0(0)
  , M0(0)
  , Mdot(0)
  {}

  Astro::Astro( string const & __name )
  : m_name(__name)
  , t0(0)
  , M0(0)
  , Mdot(0)
  {}

  Astro::Astro( Astro const & ast ) {
    *this = ast;
  }

  Astro::~Astro() {
  }

  Astro const &
  Astro::operator = ( Astro const & ast ) {
    this->m_name = ast.m_name;
    this->EQ     = ast.EQ;
    this->K      = ast.K;
    this->t0     = ast.t0; // days
    this->M0     = ast.M0; // Angle corresponding to time t0
    this->Mdot   = ast.Mdot;
    this->muS    = ast.muS;
    return *this;
  }

  void
  Astro::check_for_consistency() const {
    //ASSERT( p>0,         "Astro::check p = " << p << " must be positive!");
    //ASSERT( f*f+g*g < 1, "Astro::check f*f+g*g = " << f*f+g*g << " must be less than 1!");
    real_type h  = EQ.h;
    real_type k  = EQ.k;
    real_type e2 = h*h+k*k;
    UTILS_ASSERT(
      e2 < 1,
      "Astro::check_for_consistency h*h+k*k = {} must be less than 1!\n", e2
    );
    UTILS_ASSERT(
      muS > 0,
      "Astro::check_for_consistency muS = {} must be positive!\n", muS
    );
  }

  real_type
  Astro::latitude_of_periapsis() const {
    real_type angle = K.omega+K.Omega;
    angle_in_range(angle);
    return angle;
  }

  real_type
  Astro::latitude_of_apoapsis() const {
    real_type angle = K.omega+K.Omega+m_pi;
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
      "name: {}\n"
      "mode: {}\n"
      "t0    = {}\n"
      "M0    = {}\n"
      "Mdot  = {}\n"
      "muS   = {}\n"
      "Equinoctial\n"
      "p     = {}\n"
      "f     = {}\n"
      "g     = {}\n"
      "h     = {}\n"
      "k     = {}\n"
      "Keplerian\n"
      "e     = {}\n"
      "a     = {}\n"
      "i     = {} [degree]\n"
      "Omega = {} [degree]\n"
      "omega = {} [degree]\n"
      "Other infos"
      "period = {}",
      m_name, (EQ.retrograde?"RETROGRADE":"NORMAL"),
      t0, M0, Mdot, muS,
      EQ.p, EQ.f, EQ.g, EQ.h, EQ.k,
      K.e, K.a,
      radiants_to_degrees(K.i),
      radiants_to_degrees(K.Omega),
      radiants_to_degrees(K.omega),
      period()
    );
  }

  Astro const &
  Astro::setup_Keplerian(
    string const & n,
    real_type _t0,
    real_type _a,
    real_type _e,
    real_type _Omega, // sempre prima Omega "Grande"!!!
    real_type _omega,
    real_type _i,
    real_type _M0,
    real_type _muS
  ) {
    m_name = n;

    K.a     = _a;
    K.e     = _e;
    K.i     = _i;
    K.Omega = _Omega;
    K.omega = _omega;
    t0      = _t0;
    M0      = _M0;
    muS     = _muS;
    real_type absa = K.a > 0 ? K.a : -K.a;
    this->Mdot = sqrt(muS/absa)/absa;

    from_Keplerian_to_Equinoctial( K, EQ );
    check_for_consistency();
    return *this;
  }

  Astro const &
  Astro::setup_Equinoctial(
    string const & n,
    real_type      _t0,
    real_type      _p,
    real_type      _f,
    real_type      _g,
    real_type      _h,
    real_type      _k,
    bool           _retrograde,
    real_type      _M0,
    real_type      _muS
  ) {

    EQ.p          = _p;
    EQ.f          = _f;
    EQ.g          = _g;
    EQ.h          = _h;
    EQ.k          = _k;
    EQ.retrograde = _retrograde;

    from_equinoctial_to_Keplerian( EQ, K );

    m_name = n;
    t0     = _t0;
    M0     = _M0; // Angle corresponding to time t0
    muS    = _muS;
    real_type absa = K.a > 0 ? K.a : -K.a;
    Mdot = sqrt(muS/absa)/absa;

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
      vars("M0")         . get_number(),
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
  Astro::orbit_energy() const { return -muS/(2*K.a); }

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
    real_type M = M0 + (t-t0) * Mdot;
    if ( K.e <= 1 ) {
      mean_anomaly_to_eccentric_anomaly_elliptic( M, Mdot, K.e, E, nderiv );
    } else {
      mean_anomaly_to_eccentric_anomaly_hyperbolic( M, Mdot, K.e, E, nderiv );
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
    real_type M = M0 + (t-t0) * Mdot;
    if ( K.e <= 1 ) { // caso ellittico
      real_type E[4];
      mean_anomaly_to_eccentric_anomaly_elliptic( M, Mdot, K.e, E, nderiv );

      // da E calcolo theta = 2 * atan( sqrt(1+e)/sqrt(1-e)*tan(E/2))
      real_type theta = eccentric_anomaly_to_true_anomaly(E[0],K.e);

      Lvalues[0] = theta + K.omega; // L = theta + omega + Omega
      if ( EQ.retrograde ) Lvalues[0] -= K.Omega; // L = theta + omega + Omega
      else                 Lvalues[0] += K.Omega; // L = theta + omega + Omega

      if ( nderiv < 1 ) return;

      real_type cosE      = cos(E[0]);
      real_type tmp       = 1-K.e*cosE;
      real_type tmp1      = sqrt(1-K.e*K.e);
      real_type dtheta_dE = tmp1/tmp;
      Lvalues[1] = dtheta_dE*E[1];
      if ( nderiv < 2 ) return;

      real_type sinE        = sin(E[0]);
      real_type d2theta_d2E = -tmp1*K.e*sinE/power2(tmp);
      Lvalues[2] = dtheta_dE * E[2] + d2theta_d2E*power2(E[1]);
      if ( nderiv < 3 ) return;

      real_type d3theta_d3E = K.e*tmp1*(cosE*(cosE*K.e+1)-2*K.e)/power3(tmp);
      Lvalues[3] = d3theta_d3E*E[3]+3*d2theta_d2E*E[1]*power2(E[2])+dtheta_dE*E[3];

    } else {

      real_type F[4];
      mean_anomaly_to_eccentric_anomaly_hyperbolic( M, Mdot, K.e, F, nderiv );

      // da F calcolo theta = 2 * atan( sqrt(1+e)/sqrt(1-e)*tanh(F/2))
      real_type theta = eccentric_anomaly_to_true_anomaly( F[0], K.e );

      Lvalues[0] = theta + K.omega; // L = theta + omega + Omega
      if ( EQ.retrograde ) Lvalues[0] -= K.Omega; // L = theta + omega + Omega
      else                 Lvalues[0] += K.Omega; // L = theta + omega + Omega

      if ( nderiv < 1 ) return;

      real_type cosh_F    = cosh(F[0]);
      real_type tmp       = K.e*cosh_F-1;
      real_type tmp1      = sqrt(K.e*K.e-1);
      real_type dtheta_dF = tmp1/tmp;
      Lvalues[1] = dtheta_dF*F[1];
      if ( nderiv < 2 ) return;

      real_type tmp2        = tmp1*K.e/power2(tmp);
      real_type sinh_F      = sinh(F[0]);
      real_type d2theta_d2F = -tmp2*sinh_F;
      Lvalues[2] = dtheta_dF * F[2] + d2theta_d2F*power2(F[1]);
      if ( nderiv < 3 ) return;

      real_type d3theta_d3F = ((cosh_F*K.e+1)*cosh_F-2*K.e)*tmp2/tmp;
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
    if ( K.e < 1 ) {
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
      K.e < 1,
      "Astro::time_from_L_angle(tbase,L), can be used only for elliptic orbits\n"
    );
    real_type theta = L - K.omega;
    if ( EQ.retrograde ) theta += K.Omega;
    else                 theta -= K.Omega;

    angle_in_range(theta);
    //real_type E     = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
    //real_type E     = 2*atan2(sqrt(1-e)*sin(theta/2),sqrt(1+e)*cos(theta/2));
    //real_type M     = E-e*sin(E);

    //real_type sinE  = sqrt(1-e*e)*(sin(theta)/(1+e*cos(theta)));
    //if ( sinE > 1 ) sinE = 1; else if ( sinE < -1 ) sinE = -1;
    //real_type M     = asin(sinE)-e*sinE;
    real_type E    = 2*atan2(sqrt(1-K.e)*sin(theta/2),sqrt(1+K.e)*cos(theta/2));
    real_type M    = E-K.e*sin(E);
    real_type tres = t0 + (M-M0)/Mdot;
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
    return eccentric_anomaly_to_true_anomaly(E,K.e);
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
    real_type p    = EQ.p;
    real_type f    = EQ.f;
    real_type g    = EQ.g;
    real_type h    = EQ.h;
    real_type k    = EQ.k;
    real_type cosL = cos(L);
    real_type sinL = sin(L);
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type hk   = h*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;

    real_type I = EQ.retrograde ? -1 : 1;

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
    real_type p     = EQ.p;
    real_type f     = EQ.f;
    real_type g     = EQ.g;
    real_type h     = EQ.h;
    real_type k     = EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type hk    = h*k;
    real_type bf    = sqrt(muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cos(L)+f);
    real_type sinLg = bf * (sin(L)+g);

    real_type I = EQ.retrograde ? -1 : 1;

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
    real_type p      = EQ.p;
    real_type h      = EQ.h;
    real_type k      = EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(muS/p)/(1+h2+k2);
    real_type cosL_D = -bf*sin(L[0])*L[1];
    real_type sinL_D =  bf*cos(L[0])*L[1];

    real_type I = EQ.retrograde ? -1 : 1;

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
    real_type p      = EQ.p;
    real_type h      = EQ.h;
    real_type k      = EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(muS/p)/(1+h2+k2);
    real_type cosL_2 = bf * ( -sin(L[0])*L[2] - cos(L[0]) * power2(L[1]) );
    real_type sinL_2 = bf * (  cos(L[0])*L[2] - sin(L[0]) * power2(L[1]) );

    real_type I = EQ.retrograde ? -1 : 1;

    jx = I*2*hk*cosL_2 - (1+h2-k2)*sinL_2;
    jy = I*(1-h2+k2)*cosL_2 - 2*hk*sinL_2;
    jz = 2 * ( h*cosL_2 + I*k*sinL_2 );

  }

  //////////////////////////////////////////////////////////

  real_type
  Astro::x_position( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type p    = EQ.p;
    real_type f    = EQ.f;
    real_type g    = EQ.g;
    real_type h    = EQ.h;
    real_type k    = EQ.k;
    real_type cosL = cos(L[0]);
    real_type sinL = sin(L[0]);
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type hk   = h*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;

    real_type I = EQ.retrograde ? -1 : 1;

    return (1+h2-k2)*X+2*I*hk*Y;
  }

  real_type
  Astro::y_position( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type cosL = cos(L[0]);
    real_type sinL = sin(L[0]);
    real_type p    = EQ.p;
    real_type f    = EQ.f;
    real_type g    = EQ.g;
    real_type h    = EQ.h;
    real_type k    = EQ.k;
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type hk   = h*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;

    real_type I = EQ.retrograde ? -1 : 1;

    return I*(1-h2+k2)*Y+2*hk*X;
  }

  real_type
  Astro::z_position( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type cosL = cos(L[0]);
    real_type sinL = sin(L[0]);
    real_type p    = EQ.p;
    real_type f    = EQ.f;
    real_type g    = EQ.g;
    real_type h    = EQ.h;
    real_type k    = EQ.k;
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;

    real_type I = EQ.retrograde ? -1 : 1;

    return 2*(h*Y-I*k*X);
  }

  //////////////////////////////////////////////////////////

  real_type
  Astro::x_velocity( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type p     = EQ.p;
    real_type f     = EQ.f;
    real_type g     = EQ.g;
    real_type h     = EQ.h;
    real_type k     = EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type hk    = h*k;
    real_type bf    = sqrt(muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cos(L[0])+f);
    real_type sinLg = bf * (sin(L[0])+g);

    real_type I = EQ.retrograde ? -1 : 1;

    return I*2*hk*cosLf - (1+h2-k2)*sinLg;
  }

  real_type
  Astro::y_velocity( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type p     = EQ.p;
    real_type f     = EQ.f;
    real_type g     = EQ.g;
    real_type h     = EQ.h;
    real_type k     = EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type hk    = h*k;
    real_type bf    = sqrt(muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cos(L[0])+f);
    real_type sinLg = bf * (sin(L[0])+g);

    real_type I = EQ.retrograde ? -1 : 1;

    return I*(1-h2+k2)*cosLf - 2*hk*sinLg;
  }

  real_type
  Astro::z_velocity( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 0 );
    real_type p     = EQ.p;
    real_type f     = EQ.f;
    real_type g     = EQ.g;
    real_type h     = EQ.h;
    real_type k     = EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type bf    = sqrt(muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cos(L[0])+f);
    real_type sinLg = bf * (sin(L[0])+g);

    real_type I = EQ.retrograde ? -1 : 1;

    return 2 * ( h*cosLf + I*k*sinLg );
  }

  //////////////////////////////////////////////////////////

  real_type
  Astro::x_acceleration( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 1 );
    real_type p      = EQ.p;
    real_type h      = EQ.h;
    real_type k      = EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(muS/p)/(1+h2+k2);
    real_type cosL_D = -bf*sin(L[0])*L[1];
    real_type sinL_D =  bf*cos(L[0])*L[1];

    real_type I = EQ.retrograde ? -1 : 1;

    return I*2*hk*cosL_D - (1+h2-k2)*sinL_D;
  }

  real_type
  Astro::y_acceleration( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 1 );
    real_type p      = EQ.p;
    real_type h      = EQ.h;
    real_type k      = EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(muS/p)/(1+h2+k2);
    real_type cosL_D = -bf*sin(L[0])*L[1];
    real_type sinL_D =  bf*cos(L[0])*L[1];

    real_type I = EQ.retrograde ? -1 : 1;

    return I*(1-h2+k2)*cosL_D - 2*hk*sinL_D;
  }

  real_type
  Astro::z_acceleration( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 1 );
    real_type p      = EQ.p;
    real_type h      = EQ.h;
    real_type k      = EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type bf     = sqrt(muS/p)/(1+h2+k2);
    real_type cosL_D = -bf*sin(L[0])*L[1];
    real_type sinL_D =  bf*cos(L[0])*L[1];

    real_type I = EQ.retrograde ? -1 : 1;

    return 2 * ( h*cosL_D + I*k*sinL_D );
  }

  //////////////////////////////////////////////////////////

  real_type
  Astro::x_jerk( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 2 );
    real_type p      = EQ.p;
    real_type h      = EQ.h;
    real_type k      = EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(muS/p)/(1+h2+k2);
    real_type cosL   = cos(L[0]);
    real_type sinL   = sin(L[0]);
    real_type cosL_2 = bf * ( -sinL*L[2] - cosL * power2(L[1]) );
    real_type sinL_2 = bf * (  cosL*L[2] - sinL * power2(L[1]) );

    real_type I = EQ.retrograde ? -1 : 1;

    return I*2*hk*cosL_2 - (1+h2-k2)*sinL_2;
  }

  real_type
  Astro::y_jerk( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 2 );
    real_type p      = EQ.p;
    real_type h      = EQ.h;
    real_type k      = EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type hk     = h*k;
    real_type bf     = sqrt(muS/p)/(1+h2+k2);
    real_type cosL   = cos(L[0]);
    real_type sinL   = sin(L[0]);
    real_type cosL_2 = bf * ( -sinL*L[2] - cosL * power2(L[1]) );
    real_type sinL_2 = bf * (  cosL*L[2] - sinL * power2(L[1]) );

    real_type I = EQ.retrograde ? -1 : 1;

    return I*(1-h2+k2)*cosL_2 - 2*hk*sinL_2;
  }

  real_type
  Astro::z_jerk( real_type t ) const {
    real_type L[4];
    eval_L( t, L, 2 );
    real_type p      = EQ.p;
    real_type h      = EQ.h;
    real_type k      = EQ.k;
    real_type h2     = h*h;
    real_type k2     = k*k;
    real_type bf     = sqrt(muS/p)/(1+h2+k2);
    real_type cosL   = cos(L[0]);
    real_type sinL   = sin(L[0]);
    real_type cosL_2 = bf * ( -sinL*L[2] - cosL * power2(L[1]) );
    real_type sinL_2 = bf * (  cosL*L[2] - sinL * power2(L[1]) );

    real_type I = EQ.retrograde ? -1 : 1;

    return 2 * ( EQ.h*cosL_2 + I*EQ.k*sinL_2 );
  }

  //////////////////////////////////////////////////////////

  real_type
  Astro::ray_by_L( real_type L ) const {
    real_type p = EQ.p;
    real_type f = EQ.f;
    real_type g = EQ.g;
    return p/(1+f*cos(L)+g*sin(L));
  }

  real_type
  Astro::ray_by_L_D( real_type L ) const {
    real_type p    = EQ.p;
    real_type f    = EQ.f;
    real_type g    = EQ.g;
    real_type sinL = sin(L);
    real_type cosL = cos(L);
    return p*(f*sinL-g*cosL)/power2(1+f*cosL+g*sinL);
  }

  real_type
  Astro::ray_by_L_DD( real_type L ) const {
    real_type p     = EQ.p;
    real_type f     = EQ.f;
    real_type g     = EQ.g;
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
    real_type h  = EQ.h;
    real_type k  = EQ.k;
    real_type k2 = k*k;
    real_type h2 = h*h;
    real_type tmp = 1/(k2+h2+1);
    N[0] =  2*tmp*k;
    N[1] = -2*tmp*h;
    N[2] =  tmp*(1-k2-h2);
    if ( EQ.retrograde ) N[2] = -N[2];
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
    real_type p    = EQ.p;
    real_type f    = EQ.f;
    real_type g    = EQ.g;
    real_type h    = EQ.h;
    real_type k    = EQ.k;
    real_type cosL = cos(L);
    real_type sinL = sin(L);
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type hk   = h*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;

    real_type I = EQ.retrograde ? -1 : 1;

    real_type x = (1+h2-k2)*X+2*I*hk*Y;
    real_type y = I*(1-h2+k2)*Y+2*hk*X;
    real_type z = 2*(EQ.h*Y-I*EQ.k*X);

    real_type bf1   = sqrt(muS/p)/(1+h2+k2);
    real_type cosLf = bf1 * (cos(L)+f);
    real_type sinLg = bf1 * (sin(L)+g);

    real_type vx = I*2*hk*cosLf - (1+h2-k2)*sinLg;
    real_type vy = I*(1-h2+k2)*cosLf - 2*hk*sinLg;
    real_type vz = 2 * ( EQ.h*cosLf + I*EQ.k*sinLg );

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
    real_type h = EQ.h;
    real_type k = EQ.k;

    real_type h2 = h*h;
    real_type k2 = k*k;
    real_type hk = h*k;
    real_type bf = 1+h2+k2;

    real_type I = EQ.retrograde ? -1 : 1;

    M[0][0] = ( 1-k2+h2 ) / bf;
    M[0][1] = ( 2*hk ) / bf;
    M[0][2] = ( -2*I*EQ.k ) / bf;

    M[1][0] = ( 2*I*hk ) / bf;
    M[1][1] = ( (1+k2-h2)*I ) / bf;
    M[1][2] = ( 2*EQ.h ) / bf;

    M[2][0] = ( 2*EQ.k ) / bf;
    M[2][1] = ( -2*EQ.h ) / bf;
    M[2][2] = ( I*(1-h2-k2) ) / bf;
  }

  void
  Astro::make_retrograde() {
    if ( EQ.retrograde ) return;
    EQ.f = K.e*cos(K.omega-K.Omega);
    EQ.g = K.e*sin(K.omega-K.Omega);
    real_type tg = tan(K.i/2);
    EQ.h = cos(K.Omega)/tg;
    EQ.k = sin(K.Omega)/tg;
    EQ.retrograde = true;
  }

  void
  Astro::make_not_retrograde() {
    if ( !EQ.retrograde ) return;
    EQ.f = K.e*cos(K.omega+K.Omega);
    EQ.g = K.e*sin(K.omega+K.Omega);
    real_type tg = tan(K.i/2);
    EQ.h = tg*cos(K.Omega);
    EQ.k = tg*sin(K.Omega);
    EQ.retrograde = false;
  }
}
