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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Astro.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

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

  /*
  //   _                      _             _
  //  (_)_ ____   ____ _ _ __(_) __ _ _ __ | |_ ___
  //  | | '_ \ \ / / _` | '__| |/ _` | '_ \| __/ __|
  //  | | | | \ V / (_| | |  | | (_| | | | | |_\__ \
  //  |_|_| |_|\_/ \__,_|_|  |_|\__,_|_| |_|\__|___/
  */

  real_type
  equinoctial_to_invariant_A( Equinoctial const & EQ ) {
    return EQ.p/(1-EQ.f*EQ.f-EQ.g*EQ.g);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_invariant_E( Equinoctial const & EQ ) {
    return hypot(EQ.f,EQ.g);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_invariant_I( Equinoctial const & EQ ) {
    real_type ang = 2*atan(hypot(EQ.k,EQ.h));
    if ( EQ.retrograde ) ang = m_pi - ang;
    angle_in_range_symm(ang);
    return ang;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_invariant_lowercase_omega( Equinoctial const & EQ ) {
    real_type Omega = atan2(EQ.k,EQ.h);
    real_type ang   = atan2(EQ.g,EQ.f);
    if ( EQ.retrograde ) ang += Omega;
    else                 ang -= Omega;
    angle_in_range_symm(ang);
    return ang;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_invariant_uppercase_Omega( Equinoctial const & EQ ) {
    real_type ang = atan2(EQ.k,EQ.h);
    angle_in_range_symm(ang);
    return ang;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  equinoctial_to_invariant_Lvec(
    Equinoctial const & EQ,
    real_type           muS,
    real_type           L[3]
  ) {
    real_type p = EQ.p;
    real_type h = EQ.h;
    real_type k = EQ.k;

    real_type h2 = h*h;
    real_type k2 = k*k;
    real_type bf = sqrt(p*muS)/(1+h2+k2);
    L[0] = 2*k*bf;
    L[1] = -2*h*bf;
    L[2] = (1-(h2+k2))*bf;
    if ( EQ.retrograde ) L[2] = -L[2];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  equinoctial_to_invariant_Avec(
    Equinoctial const & EQ,
    real_type           muS,
    real_type           A[3]
  ) {
    // DA CONTROLLARE @@@@@@@@@@@@@@@@@@@@@@@@@
    real_type f  = EQ.f;
    real_type g  = EQ.g;
    real_type h  = EQ.h;
    real_type k  = EQ.k;
    real_type hk = h*k;
    real_type h2 = h*h;
    real_type k2 = k*k;
    real_type bf = muS/(1+h2+k2);
    A[0] = bf*( (1+h2-k2)*f+2*hk*g );
    A[1] = bf*( (1-h2+k2)*g+2*hk*f );
    A[2] = bf*( 2*(h*g-k*f) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_orbit_energy(
    Equinoctial const & EQ,
    real_type           muS
  ) {
    return muS*(EQ.g*EQ.g+EQ.f*EQ.f-1)/EQ.p/2;
  }

  //////////////////////////////////////////////////////////

#if 0
  void
  invariantLAEtoEquinoctial(
    real_type const L[3],
    real_type const A[3],
    real_type       E,
    real_type       muS,
    Equinoctial &   EQ
  ) {
    real_type normL2 = L[0]*L[0]+L[1]*L[1]+L[2]*L[2];
    EQ.p = normL2/muS;
    normL2 = sqrt(normL2);

    real_type e1 = (normL2+L[2])/(normL2-L[2]);
    real_type k1 = L[0]*(1+e1)/(2*normL2);
    real_type h1 = L[1]*(1+e1)/(2*normL2);
    real_type f1 = (A[0]-k1*A[2])/muS;
    real_type g1 = (A[1]+h1*A[2])/muS;
    real_type E1 = muS*(g1*g1+f1*f1-1)/EQ.p/2;

    real_type e2 = 1/e1;
    real_type k2 = L[0]*(1+e2)/(2*normL2);
    real_type h2 = L[1]*(1+e2)/(2*normL2);
    real_type f2 = (A[0]-k2*A[2])/muS;
    real_type g2 = (A[1]+h2*A[2])/muS;
    real_type E2 = muS*(g2*g2+f2*f2-1)/EQ.p/2;

    if ( abs(E-E1) < abs(E-E2) ) {
      EQ.k = k1; EQ.h = EQ.h1; EQ.f = f1; EQ.g = g1;
    } else {
      EQ.k = k2; EQ.h = EQ.h2; EQ.f = f2; EQ.g = g2;
    }
  }
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_x( Equinoctial const & EQ, real_type L ) {
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

    return (1+h2-k2)*X+2*I*hk*Y;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_y( Equinoctial const & EQ, real_type L ) {
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

    return I*(1-h2+k2)*Y+2*hk*X;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_z( Equinoctial const & EQ, real_type L ) {
    real_type p    = EQ.p;
    real_type f    = EQ.f;
    real_type g    = EQ.g;
    real_type h    = EQ.h;
    real_type k    = EQ.k;
    real_type cosL = cos(L);
    real_type sinL = sin(L);
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;

    real_type I = EQ.retrograde ? -1 : 1;

    return 2*(h*Y-I*k*X);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_vx( Equinoctial const & EQ, real_type L, real_type muS ) {
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

    return I*2*hk*cosLf - (1+h2-k2)*sinLg;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_vy( Equinoctial const & EQ, real_type L, real_type muS ) {
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

    return I*(1-h2+k2)*cosLf - 2*hk*sinLg;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_vz( Equinoctial const & EQ, real_type L, real_type muS ) {
    real_type p     = EQ.p;
    real_type f     = EQ.f;
    real_type g     = EQ.g;
    real_type h     = EQ.h;
    real_type k     = EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type bf    = sqrt(muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cos(L)+f);
    real_type sinLg = bf * (sin(L)+g);

    real_type I = EQ.retrograde ? -1 : 1;

    return 2 * ( h*cosLf + I*k*sinLg );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_radius( Equinoctial const & EQ, real_type L ) {
    return EQ.p/(1+EQ.f*cos(L)+EQ.g*sin(L));
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  equinoctial_to_velocity( Equinoctial const & EQ, real_type L, real_type muS ) {
    real_type p   = EQ.p;
    real_type f   = EQ.f;
    real_type g   = EQ.g;
    real_type tmp = muS*(f*(2*cos(L)+f)+g*(2*sin(L)+g)+1)/p;
    return sqrt(tmp);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*
  void
  angle_in_range( real_type & angle ) {
    if ( angle > 0 ) {
      integer n = integer( (angle+m_pi)/m_2pi );
      angle -= n * m_2pi;
    } else {
      integer n = integer( (-angle+m_pi)/m_2pi );
      angle += n * m_2pi;
    }
  }
  */

  real_type
  mean_anomaly_to_E( real_type M, real_type e ) {

    angle_in_range(M);

    real_type dE=0, E = M;
    for ( integer k = 0; k < 100; ++k ) { // risolvo con Newton
      dE = (E-e*sin(E)-M) / (1-e*cos(E));
      if      ( dE >  0.1 ) dE = 0.1;
      else if ( dE < -0.1 ) dE = -0.1; // scalo se passi troppo grandi
      E -= dE; // per la convergenza quando e ~1
      if ( abs( dE ) < 1E-12 ) break;
    }

    UTILS_ASSERT(
      abs( dE ) < 1E-10,
      "mean_anomaly_to_E, do not converge:"
      "\nE  = {}"
      "\ndE = {}"
      "\nM  = {}"
      "\ne  = {}\n",
      E, dE, M, e
    );
    return E;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  mean_anomaly_to_H( real_type M, real_type e ) {
    // M non va messo "in range"!!!

    real_type absM = M > 0 ? M : -M;
    real_type dH=0, H = 5*e-2.5 > absM ? pow( 6*absM/e, 1./3.) : log( 2*absM/e);
    for ( integer k = 0; k < 100; ++k ) { // risolvo con Newton
      dH = (e*sinh(H)-H-absM) / (e*cosh(H)-1);
      H -= dH;
      if ( abs( dH ) < 1E-12 ) break;
    }
    UTILS_ASSERT(
      abs( dH ) < 1E-10,
      "mean_anomaly_to_H, do not converge:"
      "\nH  = {}"
      "\ndE = {}"
      "\nM  = {}"
      "\ne  = {}\n",
      H, dH, M, e
    );
    if ( M < 0 ) H = -H;
    return H;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  eccentric_anomaly_to_true_anomaly( real_type E, real_type e ) {
    // conto quanti giri (positivi o negativi ha fatto)
    if ( e <= 1 ) {
      return 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2)); // theta
    } else {
      return 2*atan(sqrt( (e+1)/(e-1) )*tanh(E/2) ); // theta
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  E_to_true_anomaly( real_type E, real_type e ) {
    return 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2)); // theta
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  H_to_true_anomaly( real_type H, real_type e ) {
    // CONTROLLARE @@@@@@@@@@@@@@@@@@
    // conto quanti giri (positivi o negativi ha fatto)
    return 2*atan(sqrt( (e+1)/(e-1) )*tanh(H/2) ); // theta
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  true_anomaly_to_mean_anomaly( real_type theta, real_type e ) {
    real_type th2 = theta/2;
    if ( e <= 1 ) {
      real_type y = sin(th2) * sqrt( 1-e );
      real_type x = cos(th2) * sqrt( 1+e );
      real_type E = 2*atan2( y, x );
      return E-e*sin(E);
    } else {
      // 2*arctan( sqrt( (e+1)/(e-1) ) ) = theta ;
      real_type tmp = sqrt( (e-1)/(e+1) ) * tan(th2);
      UTILS_ASSERT(
        tmp > -1 && tmp < 1,
        "true_anomaly_to_mean_anomaly( theta = {}, e = {} )\n"
        "sqrt( (e-1)/(e+1) ) * tan(theta/2) = {}\n"
        "bad angle range!\n",
        theta, e, tmp
      );
      // sinh( 2* atanh( X ) ) = 2*x/(1-x^2)
      real_type H = 2*atanh( tmp );
      //real_type H = 2*atan( sqrt( (e-1)/(e+1) ) * tanh(th2) );
      return e*sinh(H)-H;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  from_Keplerian_to_Equinoctial(
    Keplerian const & K,
    Equinoctial &     EQ
  ) {
    real_type ii = K.i;
    angle_in_range_symm(ii);
    EQ.retrograde = ii > m_pi_2 || ii < -m_pi_2;

    real_type oo = K.omega;
    if ( EQ.retrograde ) oo -= K.Omega;
    else                 oo += K.Omega;
    EQ.p = K.a*(1-K.e*K.e);
    EQ.f = K.e*cos(oo);
    EQ.g = K.e*sin(oo);

    real_type ci = cos(ii/2);
    real_type si = sin(ii/2);
    real_type cO = cos(K.Omega);
    real_type sO = sin(K.Omega);
    if ( EQ.retrograde ) {
      real_type ctg = ci/si;
      EQ.h = ctg*cO;
      EQ.k = ctg*sO;
    } else {
      real_type tg = si/ci;
      EQ.h = tg*cO;
      EQ.k = tg*sO;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  from_Keplerian_to_Equinoctial(
    Keplerian const & K,
    real_type         theta,
    Equinoctial &     EQ,
    real_type &       L
  ) {
    from_Keplerian_to_Equinoctial( K, EQ );
    L = K.omega+theta;
    if ( EQ.retrograde ) L -= K.Omega;
    else                 L += K.Omega;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  from_equinoctial_to_Keplerian(
    Equinoctial const & EQ,
    Keplerian &         K
  ) {
    K.e = hypot(EQ.f,EQ.g);
    K.a = EQ.p/(1-K.e*K.e);
    real_type hk = hypot(EQ.k,EQ.h);
    K.i = 2*atan(hk);
    if ( EQ.retrograde ) K.i = m_pi - K.i;
    K.Omega = atan2(EQ.k,EQ.h);
    K.omega = atan2(EQ.g,EQ.f);
    if ( EQ.retrograde ) K.omega += K.Omega;
    else                 K.omega -= K.Omega;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  from_equinoctial_to_Keplerian(
    Equinoctial const & EQ, real_type   L,
    Keplerian &         K,  real_type & theta
  ) {
    from_equinoctial_to_Keplerian( EQ, K );
    theta = L-K.omega;
    if ( EQ.retrograde ) theta += K.Omega;
    else                 theta -= K.Omega;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  equinoctial_to_reference(
    Equinoctial const & EQ,
    real_type           f[3],
    real_type           g[3],
    real_type           w[3]
  ) {
    real_type h = EQ.h;
    real_type k = EQ.k;
    // h = g
    // k = f
    // p = h
    // q = k
    // a = p/(1-e^2)
    real_type h2 = h*h;
    real_type k2 = k*k;
    real_type hk = h*k;

    real_type I = EQ.retrograde ? -1 : 1;

    f[0] = 1-h2+k2;
    f[1] = 2*hk;
    f[2] = -2*I*h;

    g[0] = 2*I*hk;
    g[1] = I*(1+h2-k2);
    g[2] = 2*k;

    w[0] = 2*h;
    w[1] = -2*k;
    w[2] = I*(1-h2-k2);

    real_type bf = 1+h2+k2;
    f[0] /= bf; f[1] /= bf; f[2] /= bf;
    g[0] /= bf; g[1] /= bf; g[2] /= bf;
    w[0] /= bf; w[1] /= bf; w[2] /= bf;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Keplerian_to_reference(
    Keplerian const & K,
    real_type         X[3],
    real_type         Y[3],
    real_type         Z[3]
  ) {
    real_type SO = sin(K.Omega);
    real_type CO = cos(K.Omega);
    real_type so = sin(K.omega);
    real_type co = cos(K.omega);
    real_type Si = sin(K.i);
    real_type Ci = cos(K.i);

    X[0] = CO*co - SO*so*Ci;
    X[1] = SO*co + CO*so*Ci;
    X[2] = Si*so;

    Y[0] = -CO*so - SO*co*Ci;
    Y[1] = -SO*so + CO*co*Ci;
    Y[2] = Si*co;

    Z[0] = SO * Si;
    Z[1] = -CO * Si;
    Z[2] = Ci;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  point_and_velocity_to_Equinoctial_and_Keplerian(
    dvec3_t const & P,
    dvec3_t const & V,
    real_type       muS,
    Equinoctial &   EQ,
    real_type &     L,
    Keplerian &     K,
    real_type &     theta,
    real_type &     M0
  ) {

    dvec3_t n = P.cross(V);
    real_type n2 = n.squaredNorm();
    UTILS_ASSERT( n2 > 0, "n2 = {}\n", n2 );

    dvec3_t ecc = V.cross(n) / muS;

    EQ.p = n2/muS;
    n.normalize();

    EQ.retrograde = n[2] < 0;
    if ( EQ.retrograde ) {
      EQ.h = -n[1]/(1-n[2]);
      EQ.k =  n[0]/(1-n[2]);
    } else {
      EQ.h = -n[1]/(1+n[2]);
      EQ.k =  n[0]/(1+n[2]);
    }

    real_type r = P.norm();
    dvec3_t d = P/r;

    ecc -= d;

    K.i = 2*atan(hypot(EQ.k,EQ.h));
    if ( EQ.retrograde ) K.i = m_pi - K.i;
    K.Omega = atan2(EQ.k,EQ.h);
    K.e     = ecc.norm(); // NON NOTI


    real_type tmp = (EQ.h * ecc.coeff(0)+ EQ.k * ecc.coeff(1))/(K.e*hypot(EQ.h,EQ.k));
    if      ( tmp >=  1 ) K.omega = 0;
    else if ( tmp <= -1 ) K.omega = m_pi;
    else                  K.omega = acos(tmp);
    if ( ecc.coeff(2) < 0 ) K.omega = m_2pi - K.omega;

    //if ( EQ.retrograde ) K.omega += K.Omega;
    //else                 K.omega -= K.Omega;

    dvec3_t X, Y;

    real_type SO = sin(K.Omega);
    real_type CO = cos(K.Omega);
    real_type so = sin(K.omega);
    real_type co = cos(K.omega);
    real_type Si = sin(K.i);
    real_type Ci = cos(K.i);

    X.coeffRef(0) = CO*co - SO*so*Ci;
    X.coeffRef(1) = SO*co + CO*so*Ci;
    X.coeffRef(2) = Si*so;

    Y.coeffRef(0) = -CO*so - SO*co*Ci;
    Y.coeffRef(1) = -SO*so + CO*co*Ci;
    Y.coeffRef(2) = Si*co;

    // compute true anomaly
    real_type PcosTheta = P.dot(X);
    real_type PsinTheta = P.dot(Y);
    theta = atan2(PsinTheta,PcosTheta);

    L = theta + K.omega;

    if ( EQ.retrograde ) L -= K.Omega;
    else                 L += K.Omega;

    K.a = EQ.p/(1-K.e*K.e);

    real_type A = d.dot(V)*sqrt(EQ.p/muS); // equal to -cos(L)*g+f*sin(L)
    real_type B = EQ.p/r-1;                // equal to f*cos(L)+g*sin(L)
    real_type cosL = cos(L);
    real_type sinL = sin(L);

    EQ.f = A*sinL+B*cosL;
    EQ.g = B*sinL-A*cosL;

    M0 = true_anomaly_to_mean_anomaly( theta, K.e );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  equinoctial_to_point(
    Equinoctial const & EQ,
    real_type           L,
    real_type           P[3]
  ) {
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

    P[0] = (1+h2-k2)*X+2*I*hk*Y;
    P[1] = I*(1-h2+k2)*Y+2*hk*X;
    P[2] = 2*(h*Y-I*k*X);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  equinoctial_to_velocity(
    Equinoctial const & EQ,
    real_type           L,
    real_type           muS,
    real_type           V[3]
  ) {
    real_type p  = EQ.p;
    real_type f  = EQ.f;
    real_type g  = EQ.g;
    real_type h  = EQ.h;
    real_type k  = EQ.k;
    real_type h2 = h*h;
    real_type k2 = k*k;
    real_type hk = h*k;

    real_type I = EQ.retrograde ? -1 : 1;

    real_type bf1   = sqrt(muS/p)/(1+h2+k2);
    real_type cosLf = bf1 * (cos(L)+f);
    real_type sinLg = bf1 * (sin(L)+g);

    V[0] = I*2*hk*cosLf - (1+h2-k2)*sinLg;
    V[1] = I*(1-h2+k2)*cosLf - 2*hk*sinLg;
    V[2] = 2 * ( h*cosLf + I*k*sinLg );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  equinoctial_to_point_and_velocity(
    Equinoctial const & EQ,
    real_type           L,
    real_type           muS,
    real_type           P[3],
    real_type           V[3]
  ) {
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

    P[0] = (1+h2-k2)*X+2*I*hk*Y;
    P[1] = I*(1-h2+k2)*Y+2*hk*X;
    P[2] = 2*(h*Y-I*k*X);

    real_type bf1   = sqrt(muS/p)/(1+h2+k2);
    real_type cosLf = bf1 * (cosL+f);
    real_type sinLg = bf1 * (sinL+g);

    V[0] = I*2*hk*cosLf - (1+h2-k2)*sinLg;
    V[1] = I*(1-h2+k2)*cosLf - 2*hk*sinLg;
    V[2] = 2 * ( h*cosLf + I*k*sinLg );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  point_and_velocity_to_Frenet_RTN(
    dvec3_t const & P,
    dvec3_t const & V,
    dvec3_t       & Dr,
    dvec3_t       & Dt,
    dvec3_t       & Dn
  ) {
    Dr.noalias() = P;
    Dr.normalize();
    Dn = Dr.cross(V);
    Dn.normalize();
    Dt = Dn.cross( Dr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  equinoctial_to_Frenet_RTN(
    Equinoctial const & EQ,
    real_type           L,
    dvec3_t           & Dr,
    dvec3_t           & Dt,
    dvec3_t           & Dn
  ) {
    real_type h    = EQ.h;
    real_type k    = EQ.k;
    real_type cosL = cos(L);
    real_type sinL = sin(L);
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type hk   = h*k;
    real_type bf   = 1+h2+k2;
    real_type I    = EQ.retrograde ? -1 : 1;

    Dr.coeffRef(0) = ((h2-k2+1)*cosL+2*I*hk*sinL)/bf;
    Dr.coeffRef(1) = (2*hk*cosL-I*(h2-k2-1)*sinL)/bf;
    Dr.coeffRef(2) = (2*(h*sinL-cosL*k*I))/bf;

    Dn.coeffRef(0) = 2*k/bf;
    Dn.coeffRef(1) = -2*h/bf;
    Dn.coeffRef(2) = I*(1-h2-k2)/bf;

    Dt.coeffRef(0) = ((k2-h2-1)*sinL+2*hk*I*cosL)/bf;
    Dt.coeffRef(1) = (I*(1+k2-h2)*cosL-2*hk*sinL)/bf;
    Dt.coeffRef(2) = (2*k*I*sinL+2*h*cosL)/bf;
  }

  /*
  //                   _ _____
  //    _____   ____ _| |_   _|
  //   / _ \ \ / / _` | | | |
  //  |  __/\ V / (_| | | | |
  //   \___| \_/ \__,_|_| |_|
  */

  void
  equinoctial_Trtn_to_Txyz(
    Equinoctial const & EQ,
    real_type           L,
    dvec3_t     const & Trtn,
    dvec3_t           & Txyz
  ) {
    dvec3_t Dr, Dt, Dn;
    equinoctial_to_Frenet_RTN( EQ, L, Dr, Dt, Dn );
    Txyz.noalias() = Trtn.coeff(0)*Dr +
                     Trtn.coeff(1)*Dt + Trtn.coeff(2)*Dn;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  point_and_velocity_Trtn_to_Txyz(
    dvec3_t const & P,
    dvec3_t const & V,
    dvec3_t const & Trtn,
    dvec3_t       & Txyz
  ) {
    dvec3_t Dr, Dt, Dn;
    point_and_velocity_to_Frenet_RTN( P, V, Dr, Dt, Dn );
    Txyz.noalias() = Trtn.coeff(0)*Dr + Trtn[1]*Dt + Trtn[2]*Dn;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  equinoctial_Txyz_to_Trtn(
    Equinoctial const & EQ,
    real_type           L,
    dvec3_t     const & Txyz,
    dvec3_t           & Trtn
  ) {
    dvec3_t Dr, Dt, Dn;
    equinoctial_to_Frenet_RTN( EQ, L, Dr, Dt, Dn );
    Trtn.coeffRef(0) = Txyz.dot(Dr);
    Trtn.coeffRef(1) = Txyz.dot(Dt);
    Trtn.coeffRef(2) = Txyz.dot(Dn);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  point_and_velocity_Txyz_to_Trtn(
    dvec3_t const & P,
    dvec3_t const & V,
    dvec3_t const & Txyz,
    dvec3_t       & Trtn
  ) {
    dvec3_t Dr, Dt, Dn;
    point_and_velocity_to_Frenet_RTN( P, V, Dr, Dt, Dn );
    Trtn.coeffRef(0) = Txyz.dot(Dr);
    Trtn.coeffRef(1) = Txyz.dot(Dt);
    Trtn.coeffRef(2) = Txyz.dot(Dn);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // matrice
  real_type // value of b^T
  equinoctial_matrix(
    Equinoctial const & EQ,
    real_type           L,
    real_type           muS,
    real_type           A[6][3]
  ) {
    real_type p    = EQ.p;
    real_type f    = EQ.f;
    real_type g    = EQ.g;
    real_type h    = EQ.h;
    real_type k    = EQ.k;
    real_type sinL = sin(L);
    real_type cosL = cos(L);

    real_type w   = 1+(f*cosL+g*sinL);
    real_type s2  = 1+(h*h+k*k);
    real_type bf  = sqrt(p/muS);
    real_type bf1 = bf/w;
    real_type bf2 = (h*sinL-k*cosL)*bf1;

    A[0][0] = 0;        A[0][1] = bf1*2*p;            A[0][2] = 0;
    A[1][0] =  bf*sinL; A[1][1] = bf1*((w+1)*cosL+f); A[1][2] = -bf2*g;
    A[2][0] = -bf*cosL; A[2][1] = bf1*((w+1)*sinL+g); A[2][2] =  bf2*f;
    A[3][0] = 0;        A[3][1] = 0;                  A[3][2] =  bf1*s2*cosL/2;
    A[4][0] = 0;        A[4][1] = 0;                  A[4][2] =  bf1*s2*sinL/2;
    A[5][0] = 0;        A[5][1] = 0;                  A[5][2] =  bf2;

    return sqrt(p*muS)*power2(w/p);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}
