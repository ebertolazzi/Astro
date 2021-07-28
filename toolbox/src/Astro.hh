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

#pragma once

#ifndef ASTRO_dot_HH
#define ASTRO_dot_HH

#if !defined(GENERIC_CONTAINER_HH) && !defined(GENERIC_CONTAINER_dot_HH)
  #include "GenericContainer.hh"
#endif

#ifndef UTILS_dot_HH
  #include "Utils.hh"
#endif

namespace Astro {
  using GC_namespace::GenericContainer;
  using GC_namespace::ostream_type;
  using GC_namespace::real_type;
  typedef GC_namespace::int_type integer;
}

#include "Astro/Utils.hxx"
#include "Astro/Units.hxx"
#include "Astro/Kepler.hxx"
#include "Astro/Rocket.hxx"

#include <string>

namespace Astro {

  using std::string;
  using Utils::m_pi;
  using Utils::m_2pi;
  using Utils::m_pi_2;

  integer
  Lambert(
    real_type const R1[3],
    real_type const R2[3],
    real_type       tf_in, // tempo di volo
    integer         m_in,  // numero di rivoluzioni
    real_type       muC,
    real_type       V1[3],
    real_type       V2[3]
  );

  integer
  Lambert_Lancaster_Blanchard(
    real_type const r1vec[3],
    real_type const r2vec[3],
    real_type       tf_in,
    integer         m_in,
    real_type       muC,
    real_type       V1[3],
    real_type       V2[3]
  );

  void
  Lambert_minmax_distances(
    real_type const r1vec[3],
    real_type       r1,
    real_type const r2vec[3],
    real_type       r2,
    real_type       dth,
    real_type       a,
    real_type const V1[3],
    real_type const V2[3],
    integer         m,
    real_type       muC,
    real_type &     minimum_distance,
    real_type &     maximum_distance
  );

  class Astro {

    string m_name;

  protected:

    // mantengo entrambe le rappresentazioni
    Equinoctial EQ;
    Keplerian   K;

    real_type   t0;  //!< days
    real_type   M0;  //!< Angle corresponding to time t0
    real_type   Mdot;
    real_type   muS; // attorno a quale astro gira!.

  public:

    Astro( );
    Astro( string const & );
    Astro( Astro const & );

    Astro const & operator = ( Astro const & ast );

    ~Astro();

    string const & name() const { return m_name; }

    void check_for_consistency() const;

    string info() const;

    Astro const &
    info( ostream_type & stream ) const {
      stream << this->info() << '\n';
      return *this;
    }

    Astro const &
    setup_Keplerian(
      string const & n,
      real_type      t0,
      real_type      a,
      real_type      e,
      real_type      Omega,
      real_type      omega,
      real_type      i,
      real_type      M0,
      real_type      muS
    );

    Astro const &
    setup_Equinoctial(
      string const & n,
      real_type      t0,
      real_type      p,
      real_type      f,
      real_type      g,
      real_type      h,
      real_type      k,
      bool           retrograde,
      real_type      M0,
      real_type      muS
    );

    Astro const &
    setup(
      string const &    n,
      real_type         _t0,
      Keplerian const & _K,
      real_type         _M0,
      real_type         _muS
    ) {
      return setup_Keplerian(
        n, _t0, _K.a, _K.e, _K.Omega, _K.omega, _K.i, _M0, _muS
      );
    }

    Astro const &
    setup(
      string const &      n,
      real_type           _t0,
      Equinoctial const & _EQ,
      real_type           _M0,
      real_type           _muS
    ) {
      return setup_Equinoctial(
        n, _t0, _EQ.p, _EQ.f, _EQ.g, _EQ.h, _EQ.k, _EQ.retrograde, _M0, _muS
      );
    }

    bool
    setup_using_point_and_velocity(
      string    const & n,
      real_type const   P[3],
      real_type const   V[3],
      real_type         _muS,
      real_type         _t0
    ) {
      m_name = n;
      return setup_using_point_and_velocity( P, V, _muS, _t0 );
    }

    bool
    setup_using_point_and_velocity(
      real_type const P[3],
      real_type const V[3],
      real_type       _muS,
      real_type       _t0
    ) {
      real_type L0, theta;
      t0  = _t0;
      muS = _muS;
      point_and_velocity_to_Equinoctial_and_Keplerian( P, V, muS, EQ, L0, K, theta, this->M0 );
      real_type absa = std::abs(K.a);
      this->Mdot = sqrt(muS/absa)/absa;
      return true;
    }

    Astro const &
    setup( GenericContainer & vars ) {
      if ( vars.exists("i") ) setup_Keplerian( vars("name").get_string(), vars );
      else                    setup_Equinoctial( vars("name").get_string(), vars );
      return *this;
    }

    Astro const & setup_Keplerian   ( string const & n, GenericContainer & vars );
    Astro const & setup_Equinoctial ( string const & n, GenericContainer & vars );

    Equinoctial const & get_EQ() const { return this->EQ; }
    Keplerian   const & get_K()  const { return this->K;  }

    //real_type fromL0toM0( real_type L0 ) const;

    real_type mean_anomaly ( real_type t ) const { return M0 + (t-t0) * Mdot; }
    real_type true_anomaly ( real_type t ) const; // true anormality
    real_type L_orbital    ( real_type t ) const;
    real_type L_orbital_D  ( real_type t ) const;
    real_type L_orbital_DD ( real_type t ) const;
    real_type L_orbital    ( real_type t0, real_type dt ) const;

    real_type latitude_of_periapsis() const;
    real_type latitude_of_apoapsis() const;

    real_type p_orbital()  const { return EQ.p; }
    real_type f_orbital()  const { return EQ.f; }
    real_type g_orbital()  const { return EQ.g; }
    real_type h_orbital()  const { return EQ.h; }
    real_type k_orbital()  const { return EQ.k; }
    bool      retrograde() const { return EQ.retrograde; }

    real_type a_orbital()     const { return K.a; }
    real_type e_orbital()     const { return K.e; }
    real_type i_orbital()     const { return K.i; }
    real_type Omega_orbital() const { return K.Omega; }
    real_type omega_orbital() const { return K.omega; }

    real_type t0_orbital() const { return t0; }
    real_type M0_orbital() const { return M0; } // mean anomaly al time t0
    //integer numRevolution( real_type t ) const { return std::floor((t-t0)/period()-M0/m_2pi); }
    integer number_of_revolution( real_type t ) const { return integer(std::floor(mean_anomaly(t)/m_2pi)); }

    real_type time_from_L_angle( real_type t, real_type L ) const;

    real_type absolute_position( real_type t ) const { return equinoctial_to_ray(EQ,L_orbital(t)); }
    real_type absolute_velocity( real_type t ) const { return equinoctial_to_velocity(EQ,L_orbital(t),muS); }

    void eval_E( real_type t, real_type E[], integer nderiv ) const;
    void eval_L( real_type t, real_type L[], integer nderiv ) const;

    void position_by_L( real_type L, real_type & x,  real_type & y,  real_type & z  ) const;
    void velocity_by_L( real_type L, real_type & vx, real_type & vy, real_type & vz ) const;

    void position     ( real_type t, real_type & x,  real_type & y,  real_type & z  ) const;
    void velocity     ( real_type t, real_type & vx, real_type & vy, real_type & vz ) const;
    void acceleration ( real_type t, real_type & ax, real_type & ay, real_type & az ) const;
    void jerk         ( real_type t, real_type & jx, real_type & jy, real_type & jz ) const;

    real_type x_position( real_type t ) const;
    real_type y_position( real_type t ) const;
    real_type z_position( real_type t ) const;

    real_type x_velocity( real_type t ) const;
    real_type y_velocity( real_type t ) const;
    real_type z_velocity( real_type t ) const;

    real_type x_acceleration( real_type t ) const;
    real_type y_acceleration( real_type t ) const;
    real_type z_acceleration( real_type t ) const;

    real_type x_jerk( real_type t ) const;
    real_type y_jerk( real_type t ) const;
    real_type z_jerk( real_type t ) const;

    void
    position_by_L( real_type L, real_type P[3] ) const
    { this->position_by_L( L, P[0], P[1], P[2] ); }

    void
    velocity_by_L( real_type L, real_type V[3] ) const
    { this->velocity_by_L( L, V[0], V[1], V[2] ); }

    void
    position( real_type t, real_type P[3] ) const
    { this->position( t, P[0], P[1], P[2] ); }

    void
    velocity( real_type t, real_type V[3] ) const
    { this->velocity( t, V[0], V[1], V[2] ); }

    void
    acceleration( real_type t, real_type A[3] ) const
    { this->acceleration( t, A[0], A[1], A[2] ); }

    void
    jerk( real_type t, real_type J[3] ) const
    { this->jerk( t, J[0], J[1], J[2] ); }

    void
    get_Keplerian_orbital(
      real_type & _e,
      real_type & _a,
      real_type & _i,
      real_type & _Omega,
      real_type & _omega
    ) const {
      _e     = K.e;
      _a     = K.a;
      _i     = K.i;
      _Omega = K.Omega;
      _omega = K.omega;
    }

    void
    get_Equinoctial_orbital(
      real_type & _p,
      real_type & _f,
      real_type & _g,
      real_type & _h,
      real_type & _k
    ) const {
      _p = EQ.p;
      _f = EQ.f;
      _g = EQ.g;
      _h = EQ.h;
      _k = EQ.k;
    }

    real_type get_muS()   const { return muS; }
    real_type period()    const { return m_2pi/Mdot; }
    real_type apoapsis()  const { return K.a*(1+K.e); }
    real_type periapsis() const { return K.a*(1-K.e); }

    real_type
    absolute_velocity_by_angle( real_type L ) const
    { return equinoctial_to_velocity(EQ,L,muS); }

    real_type
    L_from_true_anomaly( real_type const nu ) const {
      if ( EQ.retrograde ) return nu - K.Omega + K.omega;
      else                 return nu + K.Omega + K.omega;
    }

    real_type ray_by_L   ( real_type L ) const;
    real_type ray_by_L_D ( real_type L ) const;
    real_type ray_by_L_DD( real_type L ) const;

    real_type orbit_energy() const;

    // trasformazioni varie
    void normal( real_type N[3] ) const; // normale al piano dell'ellisse
    void local_frame_by_L( real_type L, real_type M[3][3] ) const; // terna solidare al satellite
    void local_frame( real_type t, real_type M[3][3] ) const; // terna solidare al satellite
    void ellipse_frame( real_type M[3][3] ) const; // terna del piano dell'ellisse

    void make_retrograde();
    void make_not_retrograde();

  };

}

#endif
