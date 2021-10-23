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
  #include "GenericContainer/GenericContainer.hh"
#endif

#ifndef UTILS_dot_HH
  #include "Utils.hh"
#endif

namespace AstroLib {
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

namespace AstroLib {

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
    Equinoctial m_EQ;
    Keplerian   m_K;

    real_type   m_t0; //!< days
    real_type   m_M0; //!< Angle corresponding to time t0
    real_type   m_M0_theta0;
    real_type   m_M0_e;
    real_type   m_theta0;
    real_type   m_L0;
    real_type   m_Mdot;
    real_type   m_Mdot_p;
    real_type   m_Mdot_f;
    real_type   m_Mdot_g;
    real_type   m_muS; // attorno a quale astro gira!.

    void M0_Mdot_grad_eval();

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
      real_type      L0,
      real_type      muS
    );

    Astro const &
    setup(
      string const &    n,
      real_type         t0,
      Keplerian const & K,
      real_type         M0,
      real_type         muS
    ) {
      return setup_Keplerian(
        n, t0, K.a, K.e, K.Omega, K.omega, K.i, M0, muS
      );
    }

    Astro const &
    setup(
      string const &      n,
      real_type           t0,
      Equinoctial const & EQ,
      real_type           L,
      real_type           muS
    ) {
      return setup_Equinoctial(
        n, t0, EQ.p, EQ.f, EQ.g, EQ.h, EQ.k, EQ.retrograde, L, muS
      );
    }

    bool
    setup_using_point_and_velocity(
      string    const & n,
      real_type const   P[3],
      real_type const   V[3],
      real_type         muS,
      real_type         t0
    ) {
      m_name = n;
      return setup_using_point_and_velocity( P, V, muS, t0 );
    }

    bool
    setup_using_point_and_velocity(
      real_type const P[3],
      real_type const V[3],
      real_type       muS,
      real_type       t0
    );

    Astro const &
    setup( GenericContainer & vars ) {
      if ( vars.exists("i") ) setup_Keplerian( vars("name").get_string(), vars );
      else                    setup_Equinoctial( vars("name").get_string(), vars );
      return *this;
    }

    Astro const & setup_Keplerian   ( string const & n, GenericContainer & vars );
    Astro const & setup_Equinoctial ( string const & n, GenericContainer & vars );

    Equinoctial const & get_EQ() const { return m_EQ; }
    Keplerian   const & get_K()  const { return m_K;  }

    //real_type fromL0toM0( real_type L0 ) const;

    real_type mean_anomaly ( real_type t ) const { return m_M0 + (t-m_t0) * m_Mdot; }
    real_type true_anomaly ( real_type t ) const; // true anormality
    real_type L_orbital    ( real_type t ) const;
    real_type L_orbital_D  ( real_type t ) const;
    real_type L_orbital_DD ( real_type t ) const;
    real_type L_orbital    ( real_type t0, real_type dt ) const;

    real_type latitude_of_periapsis() const;
    real_type latitude_of_apoapsis() const;

    real_type p_orbital()  const { return m_EQ.p; }
    real_type f_orbital()  const { return m_EQ.f; }
    real_type g_orbital()  const { return m_EQ.g; }
    real_type h_orbital()  const { return m_EQ.h; }
    real_type k_orbital()  const { return m_EQ.k; }
    bool      retrograde() const { return m_EQ.retrograde; }

    real_type a_orbital()     const { return m_K.a; }
    real_type e_orbital()     const { return m_K.e; }
    real_type i_orbital()     const { return m_K.i; }
    real_type Omega_orbital() const { return m_K.Omega; }
    real_type omega_orbital() const { return m_K.omega; }

    real_type t0_orbital()     const { return m_t0; }
    real_type M0_orbital()     const { return m_M0; } // mean anomaly al time t0
    real_type L0_orbital()     const { return m_L0; } // mean anomaly al time t0
    real_type theta0_orbital() const { return m_theta0; } // mean anomaly al time t0
    //integer numRevolution( real_type t ) const { return std::floor((t-t0)/period()-M0/m_2pi); }
    integer number_of_revolution( real_type t ) const { return integer(std::floor(mean_anomaly(t)/m_2pi)); }

    real_type time_from_L_angle( real_type t, real_type L ) const;

    real_type absolute_position( real_type t ) const { return equinoctial_to_ray(m_EQ,L_orbital(t)); }
    real_type absolute_velocity( real_type t ) const { return equinoctial_to_velocity(m_EQ,L_orbital(t),m_muS); }

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
      real_type & e,
      real_type & a,
      real_type & i,
      real_type & Omega,
      real_type & omega
    ) const {
      e     = m_K.e;
      a     = m_K.a;
      i     = m_K.i;
      Omega = m_K.Omega;
      omega = m_K.omega;
    }

    void
    get_Equinoctial_orbital(
      real_type & p,
      real_type & f,
      real_type & g,
      real_type & h,
      real_type & k
    ) const {
      p = m_EQ.p;
      f = m_EQ.f;
      g = m_EQ.g;
      h = m_EQ.h;
      k = m_EQ.k;
    }

    real_type get_muS()   const { return m_muS; }
    real_type period()    const { return m_2pi/m_Mdot; }
    real_type apoapsis()  const { return m_K.a*(1+m_K.e); }
    real_type periapsis() const { return m_K.a*(1-m_K.e); }

    real_type
    absolute_velocity_by_angle( real_type L ) const
    { return equinoctial_to_velocity(m_EQ,L,m_muS); }

    real_type
    L_from_true_anomaly( real_type const nu ) const {
      if ( m_EQ.retrograde ) return nu - m_K.Omega + m_K.omega;
      else                   return nu + m_K.Omega + m_K.omega;
    }

    real_type radius_by_L   ( real_type L ) const;
    real_type radius_by_L_D ( real_type L ) const;
    real_type radius_by_L_DD( real_type L ) const;

    real_type orbit_energy() const;

    // trasformazioni varie
    void normal( real_type N[3] ) const; // normale al piano dell'ellisse
    void local_frame_by_L( real_type L, real_type M[3][3] ) const; // terna solidare al satellite
    void local_frame( real_type t, real_type M[3][3] ) const; // terna solidare al satellite
    void ellipse_frame( real_type M[3][3] ) const; // terna del piano dell'ellisse

    void make_retrograde();
    void make_not_retrograde();

    // gradient and jacobians
    void theta0_EQ_gradient( real_type grad[6] ) const;
    void M0_EQ_gradient( real_type grad[6] ) const;

    real_type E0_angle() const { return mean_anomaly_to_E( m_M0, m_K.e ); }
    real_type E0_EQ_gradient( real_type grad[6] ) const;

    real_type H0_angle() const { return mean_anomaly_to_H( m_M0, m_K.e ); }
    real_type H0_EQ_gradient( real_type grad[6] ) const;

    real_type E_angle( real_type t ) const { return mean_anomaly_to_E( mean_anomaly( t ), m_K.e ); }
    real_type E_EQ_gradient( real_type t, real_type grad[6] ) const;

    real_type H_angle( real_type t ) const { return mean_anomaly_to_H( mean_anomaly( t ), m_K.e ); }
    real_type H_EQ_gradient( real_type t, real_type grad[6] ) const;

    /*!
     * Compute the Jacobian of L(t;p,f,g,h,k,L0) as a function of
     * equinoctial coordinates.
     */
    real_type L_orbital_EQ_gradient( real_type t, real_type grad[6] ) const;

    // gradient respect to p, f, g, h, k, L
    void radius_EQ_gradient( real_type t, real_type grad[6] ) const;

    // gradient respect to p, f, g, h, k, l
    void absolute_velocity_EQ_gradient( real_type t, real_type grad[6] ) const;

    /*!
     * Compute the Jacobian of P(p,f,g,h,k,L(t0)) as a function of
     * equinoctial coordinates.
     */

    void position0_EQ_jacobian( real_type JP[3][6], real_type L0 ) const;
    void velocity0_EQ_jacobian( real_type JV[3][6], real_type L0 ) const;

    // must be checked
    void position_EQ_jacobian( real_type t, real_type JP[3][6] ) const;
    void velocity_EQ_jacobian( real_type t, real_type JV[3][6] ) const;

    void position_EQ_jacobian_FD( real_type t, real_type JP[3][6] ) const;
    void velocity_EQ_jacobian_FD( real_type t, real_type JV[3][6] ) const;

  };

}

#endif
