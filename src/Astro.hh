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

#if defined(__clang__)
#pragma clang diagnostic ignored "-Wsign-compare"
#elif defined(__llvm__) || defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif

#ifdef NO_SYSTEM_UTILS
  #include "Utils.hh"
  #include "Utils_eigen.hh"
  #include "Utils_Poly.hh"
  #include "Utils_autodiff.hh"
#else
  #include <Utils.hh>
  #include <Utils_eigen.hh>
  #include <Utils_Poly.hh>
  #include <Utils_autodiff.hh>
#endif

namespace AstroLib {
  using GC_namespace::GenericContainer;
  using GC_namespace::ostream_type;
  using GC_namespace::real_type;

  using integer = GC_namespace::int_type;

  using std::string;
  using std::string_view;

  using Utils::m_pi;
  using Utils::m_2pi;
  using Utils::m_pi_2;
  using Utils::power2;
  using Utils::power3;

  using dvec3_t = Eigen::Matrix<real_type,3,1>;
  using dvec_t  = Eigen::Matrix<real_type,Eigen::Dynamic,1>;
  using dmat_t  = Eigen::Matrix<real_type,Eigen::Dynamic,Eigen::Dynamic>;
 
  /*
  //               _            _ _  __  __
  //    __ _ _   _| |_ ___   __| (_)/ _|/ _|
  //   / _` | | | | __/ _ \ / _` | | |_| |_
  //  | (_| | |_| | || (_) | (_| | |  _|  _|
  //   \__,_|\__,_|\__\___/ \__,_|_|_| |_|
  */

  using autodiff::detail::DualOrder;
  template <size_t N>      using DualN = autodiff::detail::HigherOrderDual<N,real_type>;
  template <typename... T> using DualT = DualN<DualOrder<T...>::value>;

  template <size_t N> constexpr void static_check_N() { static_assert(N <= 2, "DualOrder exceeded allowed limit"); }

  #define PINS_MAP_AUTODIFF_TO_REAL1(PINS_FUN) \
  template <typename T> inline \
  DualT<T> PINS_FUN( T const & x ) { \
    constexpr size_t N{ DualOrder<T>::value }; \
    static_check_N<N>(); \
    return PINS_FUN<N>( DualT<T>(x) ); \
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //!
  //! Convert angle in degrees to radiants.
  //!
  template <typename T>
  inline
  auto 
  degrees_to_radiants( T const & deg )
  { return 0.01745329251994329576923690768488612713443*deg; } // Pi/180*deg

  //!
  //! Convert angle in radiants to degrees.
  //!
  template <typename T>
  inline
  auto
  radiants_to_degrees( T const & rad )
  { return 57.29577951308232087679815481410517033240*rad; } // (180/Pi)*rad

  /*  _             __  _              __     _
  // |_)  /\  |\ | /__ |_ __ /\  |\ | /__ |  |_
  // | \ /--\ | \| \_| |_   /--\ | \| \_| |_ |_
  */
  //! Add or remove multiple of \f$ 2\pi \f$ to an angle in order to put it in the range \f$ [0,2\pi] \f$.
  template <typename T>
  inline
  void
  angle_in_range( T & ang ) {
    using Utils::m_2pi;
    ang -= m_2pi * floor(ang/m_2pi);
    while ( ang < 0     ) ang += m_2pi;
    while ( ang > m_2pi ) ang -= m_2pi;
  }

  //! Add or remove multiple of \f$ 2\pi \f$ to an angle  in order to put it in the range \f$ [-\pi,\pi]\f$.
  template <typename T>
  inline
  void
  angle_in_range_symm( T & ang ) {
    using Utils::m_2pi;
    using Utils::m_pi;
    ang -= m_2pi * floor(ang/m_2pi);
    while ( ang < -m_pi ) ang += m_2pi;
    while ( ang >  m_pi ) ang -= m_2pi;
  }

}

#include "Astro/Units.hxx"
#include "Astro/Kepler.hxx"
#include "Astro/Rocket.hxx"

#include <string>

namespace AstroLib {

  using std::vector;

  class Astro {

    string m_name;

  protected:

    // mantengo entrambe le rappresentazioni
    Equinoctial m_EQ;
    Keplerian   m_K;

    real_type   m_t0{0}; //!< days
    real_type   m_M0{0}; //!< Angle corresponding to time t0
    real_type   m_M0_theta0{0};
    real_type   m_M0_e{0};
    real_type   m_theta0{0};
    real_type   m_L0{0};
    real_type   m_Mdot{0};
    real_type   m_Mdot_p{0};
    real_type   m_Mdot_f{0};
    real_type   m_Mdot_g{0};
    real_type   m_muS{0}; // attorno a quale astro gira!.

    void M0_Mdot_grad_eval();

  public:

    Astro( );
    Astro( string_view );
    Astro( Astro const & );

    Astro & operator = ( Astro const & ast );

    ~Astro();

    string_view name() const { return m_name; }

    bool check_for_consistency() const;

    string info() const;

    Astro const &
    info( ostream_type & stream ) const {
      stream << this->info() << '\n';
      return *this;
    }

    Astro const &
    setup_Keplerian(
      string_view n,
      real_type   t0,
      real_type   a,
      real_type   e,
      real_type   Omega,
      real_type   omega,
      real_type   i,
      real_type   M0,
      real_type   muS
    );

    Astro const &
    setup_Equinoctial(
      string_view n,
      real_type   t0,
      real_type   p,
      real_type   f,
      real_type   g,
      real_type   h,
      real_type   k,
      bool        retrograde,
      real_type   L0,
      real_type   muS
    );

    Astro const &
    setup(
      string_view       n,
      real_type         t0,
      Keplerian const & K,
      real_type         M0,
      real_type         muS
    ) {
      return setup_Keplerian( n, t0, K.a, K.e, K.Omega, K.omega, K.i, M0, muS );
    }

    Astro const &
    setup(
      string_view         n,
      real_type           t0,
      Equinoctial const & EQ,
      real_type           L,
      real_type           muS
    ) {
      return setup_Equinoctial( n, t0, EQ.p, EQ.f, EQ.g, EQ.h, EQ.k, EQ.retrograde, L, muS );
    }

    bool
    setup_using_point_and_velocity(
      string_view     n,
      dvec3_t const & P,
      dvec3_t const & V,
      real_type       muS,
      real_type       t0
    ) {
      m_name = n;
      return setup_using_point_and_velocity( P, V, muS, t0 );
    }

    bool
    setup_using_point_and_velocity(
      dvec3_t const & P,
      dvec3_t const & V,
      real_type       muS,
      real_type       t0
    );

    bool
    check_parameters( Equinoctial const & EQ ) const {
      real_type p = EQ.p;
      real_type h = EQ.h;
      real_type k = EQ.k;
      return p > 0 && h*h+k*k <= 1;
    }

    Astro const &
    setup( GenericContainer & vars ) {
      if ( vars.exists("i") ) setup_Keplerian( vars("name").get_string(), vars );
      else                    setup_Equinoctial( vars("name").get_string(), vars );
      return *this;
    }

    Astro const & setup_Keplerian   ( string_view n, GenericContainer const & vars );
    Astro const & setup_Equinoctial ( string_view n, GenericContainer const & vars );

    Equinoctial const & get_EQ() const { return m_EQ; }
    Keplerian   const & get_K()  const { return m_K;  }

    //real_type fromL0toM0( real_type L0 ) const;

    real_type mean_anomaly ( real_type t ) const { return m_M0 + (t-m_t0) * m_Mdot; }
    real_type true_anomaly ( real_type t ) const; // true anormality
    real_type L_orbital    ( real_type t ) const;
    real_type L_orbital_D  ( real_type t ) const;
    real_type L_orbital_DD ( real_type t ) const;
    real_type L_orbital    ( real_type t0, real_type dt ) const;

    void mean_anomaly_to_E( real_type M, real_type Evalues[], integer nderiv ) const;
    void mean_anomaly_to_H( real_type M, real_type Evalues[], integer nderiv ) const;

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

    real_type radius( real_type t ) const { return equinoctial_to_radius(m_EQ,L_orbital(t)); }
    real_type absolute_velocity( real_type t ) const { return equinoctial_to_velocity(m_EQ,L_orbital(t),m_muS); }

    void eval_E( real_type t, real_type E[], integer nderiv ) const;
    void eval_L( real_type t, real_type L[], integer nderiv ) const;

    void position_by_L( real_type L, real_type & x,  real_type & y,  real_type & z  ) const;
    void velocity_by_L( real_type L, real_type & vx, real_type & vy, real_type & vz ) const;

    void position     ( real_type t, real_type & x,  real_type & y,  real_type & z  ) const;
    void velocity     ( real_type t, real_type & vx, real_type & vy, real_type & vz ) const;
    void acceleration ( real_type t, real_type & ax, real_type & ay, real_type & az ) const;
    void jerk         ( real_type t, real_type & jx, real_type & jy, real_type & jz ) const;

    template <size_t N> DualN<N> x_position( DualN<N> const & t ) const;
    template <size_t N> DualN<N> y_position( DualN<N> const & t ) const;
    template <size_t N> DualN<N> z_position( DualN<N> const & t ) const;

    template <size_t N> DualN<N> x_velocity( DualN<N> const & t ) const;
    template <size_t N> DualN<N> y_velocity( DualN<N> const & t ) const;
    template <size_t N> DualN<N> z_velocity( DualN<N> const & t ) const;

    //real_type x_position( real_type t ) const;
    //real_type y_position( real_type t ) const;
    //real_type z_position( real_type t ) const;

    //real_type x_velocity( real_type t ) const;
    //real_type y_velocity( real_type t ) const;
    //real_type z_velocity( real_type t ) const;

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

    // Eigen 3 vec
    void
    position_by_L( real_type L, dvec3_t & P ) const
    { position_by_L( L, P.data() ); }

    void
    velocity_by_L( real_type L, dvec3_t & V ) const
    { velocity_by_L( L, V.data() ); }

    void
    position( real_type t, dvec3_t & P ) const
    { position( t, P.data() ); }

    void
    velocity( real_type t, dvec3_t & V ) const
    { velocity( t, V.data() ); }

    void
    acceleration( real_type t, dvec3_t & A ) const
    { acceleration( t, A.data() ); }

    void
    jerk( real_type t, dvec3_t & J ) const
    { jerk( t, J.data() ); }

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

    template <size_t N> DualN<N> radius_by_L( DualN<N> const & L ) const;

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

    real_type
    E0_angle() const {
      real_type E_values[4];
      mean_anomaly_to_E( m_M0, E_values, 0 );
      return E_values[0];
    }
    real_type E0_EQ_gradient( real_type grad[6] ) const;

    real_type
    H0_angle() const {
      real_type H_values[4];
      mean_anomaly_to_H( m_M0, H_values, 0 );
      return H_values[0];
    }
    real_type H0_EQ_gradient( real_type grad[6] ) const;

    real_type
    E_angle( real_type t ) const {
      real_type E_values[4];
      mean_anomaly_to_E( mean_anomaly( t ), E_values, 0 );
      return E_values[0];
    }
    real_type E_EQ_gradient( real_type t, real_type grad[6] ) const;

    real_type
    H_angle( real_type t ) const {
      real_type H_values[4];
      mean_anomaly_to_H( mean_anomaly( t ), H_values, 0 );
      return H_values[0];
    }
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

  integer
  Lambert(
    dvec3_t const & R1,
    dvec3_t const & R2,
    real_type       tf_in, // tempo di volo
    integer         m_in,  // numero di rivoluzioni
    real_type       muC,
    dvec3_t &       V1,
    dvec3_t &       V2
  );

  integer
  Lambert_Lancaster_Blanchard(
    dvec3_t const & r1vec,
    dvec3_t const & r2vec,
    real_type       tf_in,
    integer         m_in,
    real_type       muC,
    dvec3_t &       V1,
    dvec3_t &       V2
  );

  void
  Lambert_minmax_distances(
    dvec3_t const & r1vec,
    real_type       r1,
    dvec3_t const & r2vec,
    real_type       r2,
    real_type       dth,
    real_type       a,
    dvec3_t const & V1,
    dvec3_t const & V2,
    integer         m,
    real_type       muC,
    real_type &     minimum_distance,
    real_type &     maximum_distance
  );

  bool
  check_EQ_for_consistency(
    real_type p,
    real_type f,
    real_type g,
    real_type h,
    real_type k,
    real_type L
  );

  inline
  bool
  check_EQ_for_consistency( Equinoctial const & EQ, real_type L0 ) {
    return check_EQ_for_consistency( EQ.p, EQ.f, EQ.g, EQ.h, EQ.k, L0 );
  }

  typedef struct {
    bool      long_path;   // true se il minimo è sulla "long path"
    bool      left_branch; // true se il minimo è sul "left branch"
    real_type optimal_travel_time;
    real_type period;
    dvec3_t   W1;
    dvec3_t   W2;
  } minimum_DeltaV_extra;

  // return min |DV1|+|DV2|
  real_type
  minimum_DeltaV(
    real_type       mu,
    dvec3_t const & R1,
    dvec3_t const & V1,
    dvec3_t const & R2,
    dvec3_t const & V2,
    real_type     & DeltaV1,
    real_type     & DeltaV2,
    minimum_DeltaV_extra * extra = nullptr
  );

  // return  min |DV1|^2+|DV2|^2
  real_type
  minimum_DeltaV2(
    real_type       mu,
    dvec3_t const & R1,
    dvec3_t const & V1,
    dvec3_t const & R2,
    dvec3_t const & V2,
    real_type     & DeltaV1,
    real_type     & DeltaV2
  );

  class minimum_DeltaV_trip {
  public:
    real_type t_begin;
    real_type t_end;
    dvec3_t   P1, V1, W1;
    dvec3_t   P2, V2, W2;
    real_type DeltaV1;
    real_type DeltaV2;

    minimum_DeltaV_trip( ) { }

    explicit
    minimum_DeltaV_trip(
      real_type _t1,
      dvec3_t   _P1,
      dvec3_t   _V1,
      dvec3_t   _W1,
      real_type _t2,
      dvec3_t   _P2,
      dvec3_t   _V2,
      dvec3_t   _W2
    )
    : t_begin(_t1)
    , t_end(_t2)
    , P1(_P1), V1(_V1), W1(_W1)
    , P2(_P2), V2(_V2), W2(_W2)
    {
      DeltaV1 = (V1-W1).norm();
      DeltaV2 = (V2-W2).norm();
    }

    minimum_DeltaV_trip const &
    operator = ( minimum_DeltaV_trip const & rhs ) {
      this->t_begin = rhs.t_begin;
      this->t_end   = rhs.t_end;
      this->P1      = rhs.P1;
      this->V1      = rhs.V1;
      this->W1      = rhs.W1;
      this->P2      = rhs.P2;
      this->V2      = rhs.V2;
      this->W2      = rhs.W2;
      this->DeltaV1 = rhs.DeltaV1;
      this->DeltaV2 = rhs.DeltaV2;
      return *this;
    }

    minimum_DeltaV_trip( minimum_DeltaV_trip const & rhs ) {
      *this = rhs;
    }

  };

  inline
  bool
  operator == ( minimum_DeltaV_trip const & A,  minimum_DeltaV_trip const & B ) {
    return abs( A.t_begin-B.t_begin ) < 0.1 && abs( A.t_end-B.t_end ) < 0.1;
  }

  void
  minimum_DeltaV(
    integer                       who,
    real_type                     t_begin,
    real_type                     t_end,
    real_type                     delta_t,
    Astro const                 & a_from,
    Astro const                 & a_to,
    vector<minimum_DeltaV_trip> & trips,
    real_type                     max_accepted_DV,
    real_type                     day_tolerance,
    real_type                     day_equal,
    integer                       max_subiter
  );

  class DV_collect {
  public:
    real_type DV{0};
    real_type L_from{0};
    real_type L_to{0};

    explicit
    DV_collect( real_type _DV, real_type _L_from, real_type _L_to )
    : DV(_DV), L_from(_L_from), L_to(_L_to) {}

    DV_collect( ) { }

    DV_collect const &
    operator = ( DV_collect const & rhs ) {
      this->DV     = rhs.DV;
      this->L_from = rhs.L_from;
      this->L_to   = rhs.L_to;
      return *this;
    }

    DV_collect( DV_collect const & rhs ) {
      *this = rhs;
    }

  };

  void
  global_minimum_DeltaV(
    Astro const        & a_from,
    Astro const        & a_to,
    vector<DV_collect> & v_DV,
    integer              TABLE_SIZE,
    integer              HJ_max_iter,
    Utils::Console     * console
  );

  void
  global_minimum_DeltaV(
    Astro const                 & a_from,
    Astro const                 & a_to,
    real_type                     t_begin,
    real_type                     t_end,
    real_type                     t_tolerance,
    real_type                     max_accepted_DV,
    integer                       TABLE_SIZE,
    integer                       HJ_max_iter,
    vector<minimum_DeltaV_trip> & trips,
    Utils::Console              * console
  );

  void
  global_minimum_DeltaV2(
    Astro const        & a_from,
    Astro const        & a_to,
    vector<DV_collect> & v_DV,
    integer              TABLE_SIZE,
    integer              HJ_max_iter,
    Utils::Console     * console
  );

  void
  global_minimum_DeltaV2(
    Astro const                 & a_from,
    Astro const                 & a_to,
    real_type                     t_begin,
    real_type                     t_end,
    real_type                     t_tolerance,
    real_type                     max_accepted_DV,
    integer                       TABLE_SIZE,
    integer                       HJ_max_iter,
    vector<minimum_DeltaV_trip> & trips,
    Utils::Console              * console
  );

  real_type
  D_SH_similarity( Astro const & a, Astro const & b );

  real_type astro_x_position__xo       ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_1   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_1_1 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_1_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_1_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_1_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_1_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_1_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_1_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_2   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_2_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_2_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_2_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_2_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_2_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_2_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_3   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_3_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_3_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_3_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_3_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_3_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_4   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_4_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_4_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_4_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_4_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_5   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_5_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_5_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_5_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_6   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_6_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_6_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_7   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_x_position__xo_D_7_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );

  real_type astro_y_position__xo       ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_1   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_1_1 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_1_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_1_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_1_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_1_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_1_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_1_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_2   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_2_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_2_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_2_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_2_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_2_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_2_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_3   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_3_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_3_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_3_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_3_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_3_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_4   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_4_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_4_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_4_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_4_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_5   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_5_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_5_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_5_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_6   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_6_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_6_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_7   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_y_position__xo_D_7_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );

  real_type astro_z_position__xo       ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_1   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_1_1 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_1_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_1_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_1_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_1_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_1_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_1_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_2   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_2_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_2_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_2_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_2_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_2_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_2_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_3   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_3_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_3_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_3_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_3_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_3_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_4   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_4_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_4_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_4_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_4_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_5   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_5_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_5_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_5_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_6   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_6_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_6_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_7   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );
  real_type astro_z_position__xo_D_7_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde );

  real_type astro_x_velocity__xo       ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_1   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_1_1 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_1_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_1_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_1_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_1_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_1_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_1_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_1_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_2   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_2_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_2_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_2_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_2_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_2_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_2_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_2_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_3   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_3_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_3_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_3_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_3_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_3_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_3_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_4   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_4_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_4_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_4_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_4_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_4_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_5   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_5_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_5_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_5_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_5_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_6   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_6_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_6_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_6_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_7   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_7_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_7_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_8   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_x_velocity__xo_D_8_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );

  real_type astro_y_velocity__xo       ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_1   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_1_1 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_1_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_1_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_1_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_1_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_1_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_1_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_1_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_2   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_2_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_2_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_2_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_2_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_2_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_2_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_2_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_3   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_3_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_3_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_3_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_3_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_3_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_3_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_4   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_4_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_4_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_4_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_4_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_4_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_5   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_5_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_5_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_5_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_5_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_6   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_6_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_6_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_6_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_7   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_7_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_7_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_8   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_y_velocity__xo_D_8_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );

  real_type astro_z_velocity__xo       ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_1   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_1_1 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_1_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_1_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_1_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_1_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_1_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_1_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_1_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_2   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_2_2 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_2_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_2_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_2_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_2_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_2_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_2_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_3   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_3_3 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_3_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_3_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_3_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_3_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_3_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_4   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_4_4 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_4_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_4_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_4_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_4_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_5   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_5_5 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_5_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_5_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_5_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_6   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_6_6 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_6_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_6_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_7   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_7_7 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_7_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_8   ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );
  real_type astro_z_velocity__xo_D_8_8 ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde );

  real_type astro_ray__xo              ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_1          ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_2          ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_3          ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_4          ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_1_1        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_1_2        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_1_3        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_1_4        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_2_2        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_2_3        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_2_4        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_3_3        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_3_4        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );
  real_type astro_ray__xo_D_4_4        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L );

  real_type astro_vel__xo              ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_1          ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_2          ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_3          ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_4          ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_5          ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_1_1        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_1_2        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_1_3        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_1_4        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_1_5        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_2_2        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_2_3        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_2_4        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_2_5        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_3_3        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_3_4        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_3_5        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_4_4        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_4_5        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );
  real_type astro_vel__xo_D_5_5        ( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS );

}

#endif
