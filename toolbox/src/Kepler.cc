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
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Astro.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

namespace AstroLib {

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

  real_type
  equinoctial_to_invariant_E( Equinoctial const & EQ ) {
    return hypot(EQ.f,EQ.g);
  }

  real_type
  equinoctial_to_invariant_I( Equinoctial const & EQ ) {
    real_type ang = 2*atan(hypot(EQ.k,EQ.h));
    if ( EQ.retrograde ) ang = m_pi - ang;
    angle_in_range_symm(ang);
    return ang;
  }

  real_type
  equinoctial_to_invariant_lowercase_omega( Equinoctial const & EQ ) {
    real_type Omega = atan2(EQ.k,EQ.h);
    real_type ang   = atan2(EQ.g,EQ.f);
    if ( EQ.retrograde ) ang += Omega;
    else                 ang -= Omega;
    angle_in_range_symm(ang);
    return ang;
  }

  real_type
  equinoctial_to_invariant_uppercase_Omega( Equinoctial const & EQ ) {
    real_type ang = atan2(EQ.k,EQ.h);
    angle_in_range_symm(ang);
    return ang;
  }

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

    if ( std::abs(E-E1) < std::abs(E-E2) ) {
      EQ.k = k1; EQ.h = EQ.h1; EQ.f = f1; EQ.g = g1;
    } else {
      EQ.k = k2; EQ.h = EQ.h2; EQ.f = f2; EQ.g = g2;
    }
  }
#endif

  //////////////////////////////////////////////////////////

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

  ///////////////////////////////////////////////////////////////
  real_type
  equinoctial_to_ray( Equinoctial const & EQ, real_type L ) {
    return EQ.p/(1+EQ.f*cos(L)+EQ.g*sin(L));
  }

  real_type
  equinoctial_to_velocity( Equinoctial const & EQ, real_type L, real_type muS ) {
    real_type p   = EQ.p;
    real_type f   = EQ.f;
    real_type g   = EQ.g;
    real_type tmp = muS*(f*(2*cos(L)+f)+g*(2*sin(L)+g)+1)/p;
    return sqrt(tmp);
  }

  //////////////////////////////////////////////////////////
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

  //////////////////////////////////////////////////////////
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
  mean_anomaly_to_eccentric_anomaly_elliptic(
    real_type M,
    real_type Mdot, // derivata anomalia media
    real_type e,    // eccentricità orbita
    real_type Evalues[],
    integer   nderiv
  ) {

    angle_in_range(M);

    real_type dE=0, E = M;
    for ( integer k = 0; k < 100; ++k ) { // risolvo con Newton
      dE = (E-e*sin(E)-M) / (1-e*cos(E));
      for ( integer kk = 0; abs(dE) > 0.1 && kk < 10; ++kk ) dE /= 2; // scalo se passi troppo grandi
      E -= dE; // per la convergenza quando e ~1
      if ( std::abs( dE ) < 1E-12 ) break;
    }

    UTILS_ASSERT(
      std::abs( dE ) < 1E-10,
      "mean_anomaly_to_eccentric_anomaly_elliptic, do not converge:"
      "\nE  = {}"
      "\ndE = {}"
      "\nM  = {}"
      "\ne  = {}\n",
      E, dE, M, e
    );

    Evalues[0] = E;
    if ( nderiv < 1 ) return;

    real_type cos_E = cos(E);
    real_type bf = 1 - e*cos_E;
    Evalues[1] = Mdot/bf;
    if ( nderiv < 2 ) return;

    real_type sin_E = sin(E);
    Evalues[2] = -e*sin_E*power2(Evalues[1])/bf;
    if ( nderiv < 3 ) return;

    Evalues[3] = e*Evalues[1]*(cos_E*power2(Evalues[1])+3*sin_E*Evalues[2])/bf;
  }


  //////////////////////////////////////////////////////////
  // F = Hyperbolic Eccentric Anomaly
  // M = Mean Anomaly
  //
  // solve e*sinh(F)-F=M
  // solve e*cosh(F)F'-F'=Mdot
  // solve e*cosh(F)F''-F''=-e*sinh(F)(F')^2
  // solve e*cosh(F)F'''-F'''=-3*e*sinh(F) F' F'' -e*cosh(F)(F')^3
  //
  void
  mean_anomaly_to_eccentric_anomaly_hyperbolic(
    real_type M,    // non costante
    real_type Mdot, // derivata anomalia media
    real_type e,    // eccentricità orbita
    real_type Fvalues[],
    integer   nderiv
  ) {
    // M non va messo "in range"!!!

    real_type absM = M > 0 ? M : -M;
    real_type dF=0, F = 5*e-2.5 > absM ? pow( 6*absM/e, 1./3.) : log( 2*absM/e);
    for ( integer k = 0; k < 100; ++k ) { // risolvo con Newton
      dF = (e*sinh(F)-F-absM) / (e*cosh(F)-1);
      F -= dF;
      if ( std::abs( dF ) < 1E-12 ) break;
    }
    UTILS_ASSERT(
      std::abs( dF ) < 1E-10,
      "mean_anomaly_to_eccentric_anomaly_hyperbolic, do not converge dF {}\n", dF
    );
    if ( M < 0 ) F = -F;

    Fvalues[0] = F;
    if ( nderiv < 1 ) return;

    real_type cosh_F = cosh(F);
    real_type bf = e*cosh_F-1;
    Fvalues[1] = Mdot/bf;
    if ( nderiv < 2 ) return;

    real_type sinh_F = sinh(F);
    Fvalues[2] = -e*sinh_F*power2(Fvalues[1])/bf;
    if ( nderiv < 3 ) return;

    Fvalues[3] = -e*Fvalues[1]*(cosh_F*power2(Fvalues[1])+3*sinh_F*Fvalues[2])/bf;
  }

  //////////////////////////////////////////////////////////

  real_type
  eccentric_anomaly_to_true_anomaly( real_type E, real_type e ) {
    // conto quanti giri (positivi o negativi ha fatto)
    if ( e <= 1 ) {
      return 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2)); // theta
    } else {
      return 2*atan(sqrt( (e+1)/(e-1) )*tanh(E/2) ); // theta
    }
  }

  //////////////////////////////////////////////////////////

  real_type
  true_anomaly_to_mean_anomaly( real_type theta, real_type e ) {
    real_type th2 = theta/2;
    if ( e <= 1 ) {
      real_type y = sin(th2) * sqrt( 1-e );
      real_type x = cos(th2) * sqrt( 1+e );
      real_type E = 2*atan2( y, x );
      return E-e*sin(E);
    } else {
      // sinh( 2* atanh( X ) ) = 2*x/(1-x^2)
      //real_type F = 2*atanh( sqrt( (e-1)/(e+1) ) * tan(th2) );
      real_type F = 2*atan( sqrt( (e-1)/(e+1) ) * tanh(th2) );
      return e*sinh(F)-F;
    }
  }

  //////////////////////////////////////////////////////////

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

  void
  point_and_velocity_to_Equinoctial_and_Keplerian(
    real_type const P[3],
    real_type const V[3],
    real_type       muS,
    Equinoctial &   EQ,
    real_type &     L,
    Keplerian &     K,
    real_type &     theta,
    real_type &     M0
  ) {

    real_type n[3], d[3], ecc[3];
    cross(P,V,n);
    real_type n2 = dot3(n,n);
    UTILS_ASSERT( n2 > 0, "n2 = {}\n", n2 );

    cross(V,n,ecc);
    ecc[0] /= muS;
    ecc[1] /= muS;
    ecc[2] /= muS;

    EQ.p = n2/muS;

    real_type nlen = sqrt(n2);
    n[0] /= nlen;
    n[1] /= nlen;
    n[2] /= nlen;

    EQ.retrograde = n[2] < 0;
    if ( EQ.retrograde ) {
      EQ.h = -n[1]/(1-n[2]);
      EQ.k =  n[0]/(1-n[2]);
    } else {
      EQ.h = -n[1]/(1+n[2]);
      EQ.k =  n[0]/(1+n[2]);
    }

    real_type r = norm3(P);
    d[0] = P[0]/r;
    d[1] = P[1]/r;
    d[2] = P[2]/r;

    ecc[0] -= d[0];
    ecc[1] -= d[1];
    ecc[2] -= d[2];

    K.i = 2*atan(hypot(EQ.k,EQ.h));
    if ( EQ.retrograde ) K.i = m_pi - K.i;
    K.Omega = atan2(EQ.k,EQ.h);

    real_type tmp = (EQ.h * ecc[0] + EQ.k * ecc[1])/(norm3(ecc)*hypot(EQ.h,EQ.k));
    if      ( tmp >=  1 ) K.omega = 0;
    else if ( tmp <= -1 ) K.omega = m_pi;
    else                  K.omega = acos(tmp);
    if ( ecc[2] < 0 ) K.omega = m_2pi - K.omega;

    //if ( EQ.retrograde ) K.omega += K.Omega;
    //else                 K.omega -= K.Omega;

    real_type X[3], Y[3];

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

    // compute true anomaly
    real_type PcosTheta = dot3(X,P);
    real_type PsinTheta = dot3(Y,P);
    theta = atan2(PsinTheta,PcosTheta);

    L = theta + K.omega;

    if ( EQ.retrograde ) L -= K.Omega;
    else                 L += K.Omega;

    K.e = norm3(ecc); // NON NOTI

    K.a = EQ.p/(1-K.e*K.e);

    real_type A = dot3(d,V)*sqrt(EQ.p/muS); // equal to -cos(L)*g+f*sin(L)
    real_type B = EQ.p/r-1; // equal to f*cos(L)+g*sin(L)
    real_type cosL = cos(L);
    real_type sinL = sin(L);

    EQ.f = A*sinL+B*cosL;
    EQ.g = B*sinL-A*cosL;

    M0 = true_anomaly_to_mean_anomaly( theta, K.e );
  }

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

  void
  point_and_velocity_to_Frenet_RTN(
    real_type const P[3],
    real_type const V[3],
    real_type       Dr[3],
    real_type       Dt[3],
    real_type       Dn[3]
  ) {
    real_type r = norm3(P);
    Dr[0] = P[0]/r;
    Dr[1] = P[0]/r;
    Dr[2] = P[0]/r;
    cross( Dr, V, Dn );
    real_type len = norm3(Dn);
    Dn[0] /= len;
    Dn[1] /= len;
    Dn[2] /= len;
    cross( Dn, Dr, Dt );
  }

  void
  equinoctial_to_Frenet_RTN(
    Equinoctial const & EQ,
    real_type           L,
    real_type           Dr[3],
    real_type           Dt[3],
    real_type           Dn[3]
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

    Dr[0] = ((h2-k2+1)*cosL+2*I*hk*sinL)/bf;
    Dr[1] = (2*hk*cosL-I*(h2-k2-1)*sinL)/bf;
    Dr[2] = (2*(h*sinL-cosL*k*I))/bf;

    Dn[0] = 2*k/bf;
    Dn[1] = -2*h/bf;
    Dn[2] = I*(1-h2-k2)/bf;

    Dt[0] = ((k2-h2-1)*sinL+2*hk*I*cosL)/bf;
    Dt[1] = (I*(1+k2-h2)*cosL-2*hk*sinL)/bf;
    Dt[2] = (2*k*I*sinL+2*h*cosL)/bf;
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
    real_type   const   Trtn[3],
    real_type           Txyz[3]
  ) {
    real_type Dr[3], Dt[3], Dn[3];
    equinoctial_to_Frenet_RTN( EQ, L, Dr, Dt, Dn );
    Txyz[0] = Trtn[0]*Dr[0] + Trtn[1]*Dt[0] + Trtn[2]*Dn[0];
    Txyz[1] = Trtn[0]*Dr[1] + Trtn[1]*Dt[1] + Trtn[2]*Dn[1];
    Txyz[2] = Trtn[0]*Dr[2] + Trtn[1]*Dt[2] + Trtn[2]*Dn[2];
  }

  void
  point_and_velocity_Trtn_to_Txyz(
    real_type const P[3],
    real_type const V[3],
    real_type const Trtn[3],
    real_type       Txyz[3]
  ) {
    real_type Dr[3], Dt[3], Dn[3];
    point_and_velocity_to_Frenet_RTN( P, V, Dr, Dt, Dn );
    Txyz[0] = Trtn[0]*Dr[0] + Trtn[1]*Dt[0] + Trtn[2]*Dn[0];
    Txyz[1] = Trtn[0]*Dr[1] + Trtn[1]*Dt[1] + Trtn[2]*Dn[1];
    Txyz[2] = Trtn[0]*Dr[2] + Trtn[1]*Dt[2] + Trtn[2]*Dn[2];
  }

  void
  equinoctial_Txyz_to_Trtn(
    Equinoctial const & EQ,
    real_type           L,
    real_type   const   Txyz[3],
    real_type           Trtn[3]
  ) {
    real_type Dr[3], Dt[3], Dn[3];
    equinoctial_to_Frenet_RTN( EQ, L, Dr, Dt, Dn );
    Trtn[0] = Txyz[0]*Dr[0] + Txyz[1]*Dr[1] + Txyz[2]*Dr[2];
    Trtn[1] = Txyz[0]*Dt[0] + Txyz[1]*Dt[1] + Txyz[2]*Dt[2];
    Trtn[2] = Txyz[0]*Dn[0] + Txyz[1]*Dn[1] + Txyz[2]*Dn[2];
  }

  void
  point_and_velocity_Txyz_to_Trtn(
    real_type const P[3],
    real_type const V[3],
    real_type const Txyz[3],
    real_type       Trtn[3]
  ) {
    real_type Dr[3], Dt[3], Dn[3];
    point_and_velocity_to_Frenet_RTN( P, V, Dr, Dt, Dn );
    Trtn[0] = Txyz[0]*Dr[0] + Txyz[1]*Dr[1] + Txyz[2]*Dr[2];
    Trtn[1] = Txyz[0]*Dt[0] + Txyz[1]*Dt[1] + Txyz[2]*Dt[2];
    Trtn[2] = Txyz[0]*Dn[0] + Txyz[1]*Dn[1] + Txyz[2]*Dn[2];
  }

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

}