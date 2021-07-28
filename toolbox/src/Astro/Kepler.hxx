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

// raccolta funzioni per il calcolo di traiettorie Kepleriane

#pragma once

#ifndef ASTRO_LOCAL_KEPLER_DOT_HH
#define ASTRO_LOCAL_KEPLER_DOT_HH

namespace Astro {

  struct Equinoctial {
    real_type p;
    real_type f;
    real_type g;
    real_type h;
    real_type k;
    bool      retrograde;
    void
    info( std::ostream & stream ) const {
      fmt::print( stream,
        "p = {}\n"
        "f = {}\n"
        "g = {}\n"
        "h = {}\n"
        "k = {}\n"
        "{}\n",
        p, f, g, h, k, (retrograde?"RETROGRADE":"NORMAL")\
      );
    }
  };

  typedef struct Equinoctial Equinoctial;

  struct Keplerian {
    real_type a;
    real_type e;
    real_type i;
    real_type Omega;
    real_type omega;
    void
    info( std::ostream & stream ) const {
      fmt::print( stream,
        "a     = {}\n"
        "e     = {}\n"
        "i     = {:<10} [rad]   {} [deg]\n"
        "Omega = {:<10} [rad]   {} [deg]\n"
        "omega = {:<10} [rad]   {} [deg]\n",
        a, e,
        i, radiants_to_degrees(i),
        Omega, radiants_to_degrees(Omega),
        omega, radiants_to_degrees(omega)
      );
    }
  };

  typedef struct Keplerian Keplerian;

  /*
  //   _                      _             _
  //  (_)_ ____   ____ _ _ __(_) __ _ _ __ | |_ ___
  //  | | '_ \ \ / / _` | '__| |/ _` | '_ \| __/ __|
  //  | | | | \ V / (_| | |  | | (_| | | | | |_\__ \
  //  |_|_| |_|\_/ \__,_|_|  |_|\__,_|_| |_|\__|___/
  */

  real_type equinoctial_to_invariant_A( Equinoctial const & EQ );
  real_type equinoctial_to_invariant_E( Equinoctial const & EQ );
  real_type equinoctial_to_invariant_I( Equinoctial const & EQ );
  real_type equinoctial_to_invariant_lowercase_omega( Equinoctial const & EQ );
  real_type equinoctial_to_invariant_uppercase_Omega( Equinoctial const & EQ );

  //////////////////////////////////////////////////////////////

  void
  equinoctial_to_invariant_Lvec(
    Equinoctial const & EQ,
    real_type           muS,
    real_type           L[3]
  );

  void
  equinoctial_to_invariant_Avec(
    Equinoctial const & EQ,
    real_type           muS,
    real_type           A[3]
  );

  real_type
  equinoctial_to_orbit_energy(
    Equinoctial const & EQ,
    real_type           muS
  );

#if 0
  // da rivedere
  void
  invariantLAEtoEquinoctial(
    real_type const L[3],
    real_type const A[3],
    real_type       E,
    real_type       muS,
    Equinoctial &   EQ
  );
#endif

  //////////////////////////////////////////////////////////

  static
  inline
  real_type
  equinoctial_to_apoapsis( Equinoctial const & EQ ) {
    real_type e = hypot(EQ.f,EQ.g);
    return EQ.p/(1-e);
  }

  static
  inline
  real_type
  equinoctial_to_periapsis( Equinoctial const & EQ ) {
    real_type e = hypot(EQ.f,EQ.g);
    return EQ.p/(1+e);
  }

  //////////////////////////////////////////////////////////

  real_type equinoctial_to_x( Equinoctial const & EQ, real_type L );
  real_type equinoctial_to_y( Equinoctial const & EQ, real_type L );
  real_type equinoctial_to_z( Equinoctial const & EQ, real_type L );
  real_type equinoctial_to_vx( Equinoctial const & EQ, real_type L, real_type muS );
  real_type equinoctial_to_vy( Equinoctial const & EQ, real_type L, real_type muS );
  real_type equinoctial_to_vz( Equinoctial const & EQ, real_type L, real_type muS );

  ///////////////////////////////////////////////////////////////
  real_type equinoctial_to_ray( Equinoctial const & EQ, real_type L );
  real_type equinoctial_to_velocity( Equinoctial const & EQ, real_type L, real_type muS );

  ///////////////////////////////////////////////////////////////

  void
  mean_anomaly_to_eccentric_anomaly_elliptic(
    real_type M,    // non costante
    real_type Mdot, // derivata anomalia media
    real_type e,    // eccentricità orbita
    real_type Evalues[],
    integer   nderiv
  );

  void
  mean_anomaly_to_eccentric_anomaly_hyperbolic(
    real_type M,    // non costante
    real_type Mdot, // derivata anomalia media
    real_type e,    // eccentricità orbita
    real_type Fvalues[],
    integer   nderiv
  );

  ///////////////////////////////////////////////////////////////

  real_type
  eccentric_anomaly_to_true_anomaly( real_type E, real_type e );

  real_type
  true_anomaly_to_mean_anomaly( real_type theta, real_type e );

  //////////////////////////////////////////////////////////

  void
  from_Keplerian_to_Equinoctial(
    Keplerian const & K,
    Equinoctial &     EQ
  );

  void
  from_Keplerian_to_Equinoctial(
    Keplerian const & K,  real_type   theta,
    Equinoctial &     EQ, real_type & L
  );

  void
  from_equinoctial_to_Keplerian(
    Equinoctial const & EQ,
    Keplerian &         K
  );

  void
  from_equinoctial_to_Keplerian(
    Equinoctial const & EQ, real_type   L,
    Keplerian & K,          real_type & theta
  );

  void
  equinoctial_to_reference(
    Equinoctial const & EQ,
    real_type           f[3],
    real_type           g[3],
    real_type           w[3]
  );

  void
  Keplerian_to_reference(
    Keplerian const & K,
    real_type         X[3],
    real_type         Y[3],
    real_type         Z[3]
  );

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
  );

  void
  equinoctial_to_point(
    Equinoctial const & EQ,
    real_type           L,
    real_type           P[3]
  );

  void
  equinoctial_to_velocity(
    Equinoctial const & EQ,
    real_type           L,
    real_type           muS,
    real_type           V[3]
  );

  void
  equinoctial_to_point_and_velocity(
    Equinoctial const & EQ,
    real_type           L,
    real_type           muS,
    real_type           P[3],
    real_type           V[3]
  );

  void
  point_and_velocity_to_Frenet_RTN(
    real_type const P[3],
    real_type const V[3],
    real_type       Dr[3],
    real_type       Dt[3],
    real_type       Dn[3]
  );

  void
  equinoctial_to_Frenet_RTN(
    Equinoctial const & EQ,
    real_type           L,
    real_type           Dr[3],
    real_type           Dt[3],
    real_type           Dn[3]
  );

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
  );

  void
  point_and_velocity_Trtn_to_Txyz(
    real_type const P[3],
    real_type const V[3],
    real_type const Trtn[3],
    real_type       Txyz[3]
  );

  void
  equinoctial_Txyz_to_Trtn(
    Equinoctial const & EQ,
    real_type           L,
    real_type   const   Txyz[3],
    real_type           Trtn[3]
  );

  void
  point_and_velocity_Txyz_to_Trtn(
    real_type const P[3],
    real_type const V[3],
    real_type const Txyz[3],
    real_type       Trtn[3]
  );

  // matrice
  real_type // value of b^T
  equinoctial_matrix(
    Equinoctial const & EQ,
    real_type           L,
    real_type           muS,
    real_type           A[6][3]
  );

}

#endif
