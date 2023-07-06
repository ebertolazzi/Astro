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
//#include "Utils_NelderMead.hh"
#include "Utils_HJPatternSearch.hh"
#include "PolynomialRoots.hh"

namespace AstroLib {

  using std::acos;
  using std::abs;
  using std::min;
  using std::max;
  using std::atan2;
  using std::floor;
  using std::ceil;

  /*\
   |
   |  Optimal two-impulse rendezvous using constrained multiple-revolution Lambert solutions
   |  Gang Zhang · Di Zhou · Daniele Mortari
   |  2011
   |
   |  Celest Mech Dyn Astr (2011) 110:305-317
   |  DOI 10.1007/s10569-011-9349-z
   |
  \*/

  real_type
  minimum_DeltaV(
    real_type       mu,
    dvec3_t const & R1,
    dvec3_t const & V1,
    dvec3_t const & R2,
    dvec3_t const & V2,
    real_type     & DeltaV1,
    real_type     & DeltaV2,
    minimum_DeltaV_extra * extra
  ) {

    using std::abs;

    static real_type tolerance = 1e-9;

    dvec3_t   ch = R2-R1;
    real_type c  = ch.norm();
    ch          /= c;

    real_type r1  = R1.norm();
    real_type r2  = R2.norm();
    dvec3_t   rr1 = R1/r1;
    dvec3_t   rr2 = R2/r2;

    // coefficienti del polinomio le cui radici sono candidati minimi

    real_type V__1   = V1.norm();
    real_type V__1r  = V1.dot(rr1);
    real_type V__1c  = V1.dot(ch);

    real_type V__2   = V2.norm();
    real_type V__2r  = V2.dot(rr2);
    real_type V__2c  = V2.dot(ch);

    real_type K           = mu*c/(r1*r2+R1.dot(R2)); // equation (7)
    real_type cos__phi__1 = ch.dot(rr1);
    real_type cos__phi__2 = ch.dot(rr2);

    real_type t1 = K*K;
    real_type t2 = t1*t1;
    real_type t3 = t2*K;
    real_type t8 = V__1*V__1;
    real_type t10 = V__1r*V__1r;
    real_type t12 = V__2*V__2;
    real_type t14 = V__2r*V__2r;

    real_type poly0 = t10*t2+t12*t2-t14*t2-t8*t2-2.0*cos__phi__1*t3-2.0*cos__phi__2*t3;

    real_type t22 = t1*K;
    real_type t26 = t10*t22;
    real_type t29 = V__1r*t22;

    real_type poly1 = 4.0*cos__phi__2*V__1r*t2-4.0*cos__phi__1*V__2r*t2-2.0*V__2r*t8*t22+4.0*V__1c*t2-4.0*V__2c*t2+2.0*V__2r*t26-2.0*t12*t29+2.0*t14*t29;

    real_type t43 = t8*t1;
    real_type t45 = t10*t1;
    real_type t47 = V__1c*t22;
    real_type t54 = V__2c*t22;

    real_type poly2 = -2.0*cos__phi__1*t14*t22-2.0*V__1r*t47+8.0*V__2c*t29+8.0*V__2r*t47-2.0*V__2r*t54+t12*t45-t14*t43-2.0*cos__phi__2*t26;

    real_type t63 = V__1c*t1;
    real_type t73 = V__1r*t1;
    real_type t74 = V__2c*V__2r;

    real_type poly3 = -4.0*V__1r*V__2r*t63-2.0*V__2c*t43-2.0*V__2c*t45-4.0*V__2r*t22+2.0*t12*t63+2.0*t14*t63-4.0*cos__phi__2*t47-4.0*cos__phi__1*t54+4.0*t74*t73-4.0*t29;

    real_type t83 = V__2c*t1;
    real_type t87 = t8*K;
    real_type t90 = K*V__1c;
    real_type t99 = V__1c*V__1c;
    real_type t103 = V__2c*V__2c;

    real_type poly4 = -2.0*t12*V__1r*t90+4.0*V__1r*cos__phi__2*t63-4.0*V__2r*cos__phi__1*t83-t103*t1-2.0*t12*t1-t14*t1+t99*t1+4.0*cos__phi__1*t22+4.0*cos__phi__2*t22-2.0*t74*t87+2.0*t43+t45;

    real_type t113 = t99*K;
    real_type t121 = K*V__1r;

    real_type poly5 = 4.0*V__1r*V__2c*t90+4.0*cos__phi__1*V__2r*t1+2.0*V__2r*t113+2.0*V__2r*t87+2.0*t103*t121+2.0*t12*t121-4.0*cos__phi__2*t73+4.0*t74*t90-4.0*t63+4.0*t83;

    real_type t141 = K*V__2c;

    real_type poly6 = -2.0*cos__phi__1*t103*K-2.0*V__1r*t90-8.0*V__2c*t121-2.0*V__2r*t141-8.0*V__2r*t90-t103*t8-2.0*cos__phi__2*t113+t99*t12;
    real_type poly7 = 4.0*K*V__2r+2.0*t103*V__1c-2.0*t12*V__1c+2.0*V__2c*t8-2.0*V__2c*t99+4.0*cos__phi__1*t141+4.0*cos__phi__2*t90+4.0*t121;
    real_type poly8 = -2.0*K*cos__phi__1-2.0*K*cos__phi__2-t103+t12-t8+t99;

    std::vector<real_type> roots;
    roots.clear();
    roots.reserve(20);

    // cerco soluzione
    Utils::Sturm<real_type> sturm;
    Utils::Poly<real_type>::dvec_t coeffs8;
    coeffs8.resize(9); // ordine 9
    coeffs8 << poly0, poly1, poly2, poly3, poly4, poly5, poly6, poly7, poly8;
    for ( int i=0; i < 9;++i ) if ( !Utils::is_finite(coeffs8.coeff(i)) ) return Utils::Inf<real_type>();
    Utils::Poly<real_type> poly(coeffs8);

    {
      poly.normalize();
      sturm.build( poly );

      integer n_roots = sturm.separate_roots();
      sturm.refine_roots();

      Utils::Poly<real_type>::dvec_t R = sturm.roots();
      for ( integer k = 0; k < n_roots; ++k ) roots.push_back( R.coeff(k) );
    }
    // cerco eventuali radici doppie
    {
      Utils::Poly<real_type>::dvec_t coeffs7;
      coeffs7.resize(8); // ordine 8
      coeffs7 << poly1, 2*poly2, 3*poly3, 4*poly4, 5*poly5, 6*poly6, 7*poly7, 8*poly8;

      Utils::Poly<real_type> dpoly(coeffs7);
      dpoly.normalize();
      sturm.build( dpoly );

      integer n_roots = sturm.separate_roots();
      sturm.refine_roots();

      Utils::Poly<real_type>::dvec_t R = sturm.roots();
      for ( integer k = 0; k < n_roots; ++k ) {
        if ( abs(poly.eval(R.coeff(k))) < tolerance ) roots.push_back( R.coeff(k) );
      }
    }

    //fmt::print( "N.roots = {}\nCheck\n", n_roots );
    //for ( auto & x : roots) {
    //  fmt::print( "P({}) = {}\n", x, poly.eval(x) );
    //}

    real_type minDV = Utils::Inf<real_type>();
    for ( real_type & omega__c : roots ) {

      real_type t1 = V__1*V__1;
      real_type t4 = omega__c*omega__c;
      real_type t5 = K*K;
      real_type t6 = 1/t4;
      real_type t7 = t6*t5;
      real_type t11 = K/omega__c;
      real_type t14 = V__2*V__2;

      real_type DV1_2 = 2*(K*cos__phi__1-V__1c*omega__c-V__1r*t11) + t1+t4+t7;
      real_type DV2_2 = 2*(V__2r*t11-K*cos__phi__2-V__2c*omega__c) + t14+t4+t7;

      real_type t22 = 1/t4/omega__c;
      real_type t23 = t22*t5;
      real_type t24 = t6*K;

      real_type DV1_D = 2*(V__1r*t24-V__1c-t23+omega__c);
      real_type DV2_D = 2*(omega__c-V__2r*t24-V__2c-t23);

      real_type t29 = t4*t4;
      real_type t32 = 6*t5/t29;
      real_type t33 = t22*K;

      real_type DV1_DD = 2 - 4*V__1r*t33 + t32;
      real_type DV2_DD = 2 + 4*V__2r*t33 + t32;

      real_type DV1 = sqrt(DV1_2);
      real_type DV2 = sqrt(DV2_2);

      real_type F  = DV1+DV2;
      real_type DF = (DV1_D/DV1+DV2_D/DV2)/2;

      real_type DDF = (DV1_DD/DV1+DV2_DD/DV2)/2
                    - ( (DV1_D*DV1_D)/(DV1*DV1_2) + (DV2_D*DV2_D)/(DV2*DV2_2) )/4;

      // fmt::print( "F = {:10.4} DF = {:10.4} DDF = {:10.4}\n", F, DF, DDF );

      if ( abs(DF) < tolerance && DDF >= 0 && F < minDV ) { // trovato minimo
        minDV   = F;
        DeltaV1 = DV1;
        DeltaV2 = DV2;
        if ( extra != nullptr ) {
          bool      & long_path   = extra->long_path;
          bool      & left_branch = extra->left_branch;
          real_type & ott         = extra->optimal_travel_time;
          real_type & period      = extra->period;
          dvec3_t   & W1          = extra->W1;
          dvec3_t   & W2          = extra->W2;

          real_type Kc = K/omega__c;
          W1.noalias() = omega__c*ch+Kc*rr1;
          W2.noalias() = omega__c*ch-Kc*rr2;

          // calcolo alpha e beta
          real_type a     = mu*r1/(2*mu-r1*W1.squaredNorm()); // equation (18)
          real_type s     = (r1+r2+c)/2;
          real_type alpha = acos(1-s/a);     // equation  (2)
          real_type beta  = acos(1-(s-c)/a); // equation  (2)
          // Optimal two-impulse rendezvous using constrained multiple-revolution Lambert solutions
          // Gang Zhang - Di Zhou -  Daniele Mortari
          // per determinare long/short part calcolo p (eq. 6)
          dvec3_t r1_x_r2 = R1.cross(R2);
          real_type p = (omega__c/c);
          p = r1_x_r2.squaredNorm()*p*p/mu;

          // calcolo pm da equazione (21)
          real_type pm = 2*(s-r1)*(s-r2)/c;
          long_path   = p <= pm;
          left_branch = omega__c < 0; // da equazione (7) segno(omegac)=sign(sin(theta))
          if ( left_branch ) {
            beta = -beta;
            if ( !long_path ) alpha = m_2pi - alpha;
          } else {
            if ( long_path ) alpha = m_2pi - alpha;
          }
          period = a*sqrt(a/mu); // periodo / (2*pi)
          ott = (alpha-sin(alpha))-(beta-sin(beta)); // angolo della traiettoria
          if ( ott < 0 ) ott += m_2pi;
          ott    *= period;
          period *= m_2pi;
        }
      }
    }
    return minDV;
  }

  //
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //

  /*\
   |
   |  A closed-form solution to the minimum deltaV^2 problem Lambert's problem
   |  Martín Avendaño, Daniele Mortari
   |  2009
   |
   |  Celest Mech Dyn Astr (2010) 106:25-37
   |  DOI 10.1007/s10569-009-9238-x
   |
  \*/

  real_type
  minimum_DeltaV2(
    real_type       mu,
    dvec3_t const & R1,
    dvec3_t const & V1,
    dvec3_t const & R2,
    dvec3_t const & V2,
    real_type     & DeltaV1,
    real_type     & DeltaV2
  ) {

    static real_type tolerance = 1e-9;

    dvec3_t   ch = R2-R1;
    real_type c  = ch.norm();
    ch          /= c;

    real_type r1  = R1.norm();
    real_type r2  = R2.norm();
    dvec3_t   rr1 = R1/r1;
    dvec3_t   rr2 = R2/r2;

    // coefficienti del polinomio le cui radici sono candidati minimi

    real_type V__1  = V1.norm();
    real_type V__1r = V1.dot(rr1);
    real_type V__1c = V1.dot(ch);

    real_type V__2  = V2.norm();
    real_type V__2r = V2.dot(rr2);
    real_type V__2c = V2.dot(ch);

    real_type K           = mu*c/(r1*r2+R1.dot(R2)); // equation (7)
    real_type cos__phi__1 = ch.dot(rr1);
    real_type cos__phi__2 = ch.dot(rr2);

    real_type poly0 = -4*K*K;
    real_type poly1 = 2*K*(V__1r-V__2r);
    real_type poly2 = 0;
    real_type poly3 = -2*(V__1c+V__2c);
    real_type poly4 = 4;

    if ( ! ( Utils::is_finite(poly0) &&
             Utils::is_finite(poly1) &&
             Utils::is_finite(poly3) ) ) return Utils::Inf<real_type>();

    PolynomialRoots::Quartic quartic( poly4, poly3, poly2, poly1, poly0 );

    real_type roots[4];
    integer n_roots = quartic.get_real_roots( roots );
    real_type minDV2 = Utils::Inf<real_type>();
    for ( int k = 0; k < n_roots; ++k ) {
      real_type omega__c = roots[k];

      real_type t1    = V__1*V__1;
      real_type t4    = omega__c*omega__c;
      real_type t5    = K*K;
      real_type t6    = 1/t4;
      real_type t7    = t6*t5;
      real_type t11   = K/omega__c;
      real_type t15   = V__2*V__2;
      real_type t25   = 1/t4/omega__c;
      real_type t29   = t6*K;
      real_type t36   = t4*t4;
      real_type t40   = t25*K;
      real_type DV1   = 2*(K*cos__phi__1-V__1c*omega__c-V__1r*t11)+t1+t4+t7;
      real_type DV2   = 2*(V__2r*t11-K*cos__phi__2-V__2c*omega__c)+t15+t4+t7;
      real_type DV    = DV1+DV2;
      real_type DV_D  = 2*((V__1r-V__2r)*t29-V__1c-V__2c)+4*(omega__c-t25*t5);
      real_type DV_DD = 4+12*t5/t36+4*(V__2r-V__1r)*t40;

      //fmt::print( "DV = {} DV_D = {} DV_DD = {}\n", DV, DV_D, DV_DD );
      if ( abs(DV_D) < tolerance && DV_DD >= 0 && DV < minDV2 ) {
        minDV2  = DV; // trovato minimo
        DeltaV1 = sqrt(DV1);
        DeltaV2 = sqrt(DV2);
      }
    }
    return minDV2;
  }

  // ---------------------------------------------------------------------------

  #define TABLE_SIZE 16

  //!
  //!  Given two plates/asteroids find the best possible theoretical DV.
  //!
  //!  \param[IN]  a_from        starting planet/asteroid
  //!  \param[IN]  a_to          arrival planet/asteroid
  //!  \param[IN]  day_tolerance tolerance for pattern search
  //!  \param[IN]  console       pointer to console for printing
  //!  \return  the minimal DV
  //!
  void
  global_minimum_DeltaV(
    Astro const        & a_from,
    Astro const        & a_to,
    vector<DV_collect> & v_DV,
    Utils::Console     * console
  ) {

    real_type L_tolerance = m_pi/360;
    real_type L_equal     = 5*L_tolerance;
    integer   HJ_max_iter = 200;

    v_DV.clear();
    v_DV.reserve(10);

    real_type L_table[TABLE_SIZE];
    real_type DV_table[TABLE_SIZE][TABLE_SIZE];

    real_type muS = a_from.get_muS();

    real_type DeltaV1, DeltaV2;

    dvec3_t R1[TABLE_SIZE], V1[TABLE_SIZE],
            R2[TABLE_SIZE], V2[TABLE_SIZE],
            W1, W2;

    // prima tabella
    real_type delta_L = m_2pi/TABLE_SIZE;
    for ( integer i = 0 ; i < TABLE_SIZE ; ++i ) {
      L_table[i] = i*delta_L;
      a_from.position_by_L( L_table[i], R1[i] );
      a_from.velocity_by_L( L_table[i], V1[i] );
      a_to.position_by_L( L_table[i], R2[i] );
      a_to.velocity_by_L( L_table[i], V2[i] );
    }

    real_type L_from_min, L_from_max, L_to_min, L_to_max;
    auto fun = [ &a_from, &a_to, &L_from_min, &L_from_max, &L_to_min, &L_to_max, muS ] ( real_type const X[] )->real_type {
      real_type La = X[0];
      real_type Lb = X[1];
      if ( La < L_from_min || La > L_from_max || Lb < L_to_min || Lb > L_to_max ) Utils::Inf<real_type>();
      dvec3_t P1, V1, P2, V2;
      real_type DV1, DV2;
      a_from.position_by_L( La, P1 );  a_from.velocity_by_L( La, V1 );
      a_to.position_by_L( Lb, P2 );    a_to.velocity_by_L( Lb, V2 );
      return minimum_DeltaV( muS, P1, V1, P2, V2, DV1, DV2, nullptr );
    };
    std::function<real_type(real_type const[])> F(fun);

    // riempio tabella 8*8
    for ( integer i = 0; i < TABLE_SIZE; ++i ) {
      for ( integer j = 0; j < TABLE_SIZE; ++j ) {
        DV_table[i][j] = minimum_DeltaV( muS, R1[i], V1[i], R2[j], V2[j], DeltaV1, DeltaV2, nullptr );
      }
    }

    // cerco minimi locali
    for ( integer i = 0; i < TABLE_SIZE; ++i ) {

      integer im1{i-1}, ip1{i+1};
      if      ( im1 < 0           ) im1 = TABLE_SIZE-1;
      else if ( ip1 >= TABLE_SIZE ) ip1 = 0;

      integer ii[] = { im1, i, ip1,
                       im1,    ip1,
                       im1, i, ip1 };

      real_type La = L_table[i];
      L_from_min = La-2*delta_L;
      L_from_max = La+2*delta_L;

      for ( integer j = 0 ; j < TABLE_SIZE ; ++j ) {
        real_type DV = DV_table[i][j] ;

        // calcolo primi vicini
        integer jm1{j-1}, jp1{j+1};
        if      ( jm1 < 0           ) jm1 = TABLE_SIZE-1;
        else if ( jp1 >= TABLE_SIZE ) jp1 = 0;

        integer jj[] = { jm1, jm1, jm1,
                         j,        j,
                         jp1, jp1, jp1 };

        // controllo se minimo locale
        bool ok_min = true;
        for ( integer k=0 ; k < 8 && ok_min ; ++k )
          ok_min = DV <= DV_table[ ii[k] ][ jj[k] ];

        if ( !ok_min ) continue ;

        // trovato minimo locale, raffino con Nelder Mead
        real_type Lb = L_table[j];
        L_to_min = Lb-2*delta_L;
        L_to_max = Lb+2*delta_L;

        //Utils::NelderMead<real_type> solver("NMsolver");
        Utils::HJPatternSearch<real_type> solver("HJPatternSearch");
        solver.setup( 2, F, console );
        solver.set_tolerance( L_tolerance );
        solver.set_max_iterations(HJ_max_iter);
        real_type X[2] = {La,Lb};
        solver.run( X, delta_L );
        DV = solver.get_last_solution( X );
        // discard probaby double solutions
        if ( abs(La-X[0]) > 1.5*delta_L || abs(Lb-X[1]) > 1.5*delta_L ) continue;
        v_DV.emplace_back(DV,X[0],X[1]);
      }
    }

    std::sort(
      v_DV.begin(), v_DV.end(),
      [L_equal]( DV_collect const & A, DV_collect const & B ) -> bool {
        real_type dL = B.L_from-A.L_from;
        if ( abs(dL) < L_equal ) return B.L_to-A.L_to;
        else                     return dL > 0;
      }
    );

    // remove duplicates
    v_DV.erase(
      std::unique(
        v_DV.begin(), v_DV.end(),
        [L_equal]( DV_collect const & A, DV_collect const & B ) -> bool {
          real_type dL1 = B.L_from-A.L_from;
          real_type dL2 = B.L_to-A.L_to;
          return abs(dL1) <= L_equal && abs(dL2) <= L_equal;
        }
      ),
      v_DV.end()
    );
  }

  // ==========================================================================
  // ==========================================================================
  // ==========================================================================

  void
  global_minimum_DeltaV(
    Astro const                 & a_from,
    Astro const                 & a_to,
    real_type                     t_begin,
    real_type                     t_end,
    real_type                     t_tolerance,
    real_type                     max_accepted_DV,
    vector<minimum_DeltaV_trip> & trips,
    Utils::Console              * console
  ) {

    using std::min;

    real_type muS = a_from.get_muS();
    integer   HJ_max_iter = 200;

    vector<DV_collect> v_DV;
    global_minimum_DeltaV( a_from, a_to, v_DV, console );
    // loop sui candidati
    real_type from_period = a_from.period();
    real_type to_period   = a_to.period();
    trips.clear();
    trips.reserve( 4*v_DV.size());

    real_type t_from_min, t_from_max, t_to_min, t_to_max;
    auto fun = [ &a_from, &a_to, &t_from_min, &t_from_max, &t_to_min, &t_to_max, muS ] ( real_type const X[] )->real_type {
      real_type ta = X[0];
      real_type tb = X[1];
      if ( ta < t_from_min || ta > t_from_max || tb < t_to_min || tb > t_to_max ) Utils::Inf<real_type>();
      dvec3_t P1, V1, VV1, P2, V2, VV2;
      a_from.position( ta, P1 );  a_from.velocity( ta, V1 );
      a_to.position( tb, P2 );    a_to.velocity( tb, V2 );
      integer okk = Lambert( P1, P2, tb-ta, 0, muS, VV1, VV2 );
      if ( okk != 1 ) return Utils::Inf<real_type>();
      return (V1-VV1).norm()+(V2-VV2).norm();
    };
    std::function<real_type(real_type const[])> F(fun);

    real_type DeltaT = min(from_period,to_period)/16;

    for ( DV_collect const & S : v_DV ) {

      if ( S.DV > max_accepted_DV ) continue; // skip

      real_type a_t0 = a_from.time_from_L_angle( t_begin, S.L_from );
      real_type b_t0 = a_to.time_from_L_angle( t_begin, S.L_to );

      for ( real_type ta = a_t0; ta < t_end; ta += from_period ) {

        t_from_min = ta - DeltaT;
        t_from_max = ta + DeltaT;

        for ( real_type tb = b_t0; tb < t_end; tb += to_period ) {
          if ( tb <= ta ) continue;

          t_to_min = tb - DeltaT;
          t_to_max = tb + DeltaT;

          Utils::HJPatternSearch<real_type> solver("HJPatternSearch");
          solver.setup( 2, F, console );
          solver.set_tolerance( t_tolerance );
          solver.set_max_iterations(HJ_max_iter);
          real_type X[2] = {ta,tb};
          solver.run( X, DeltaT/2 );
          real_type DV = solver.get_last_solution( X );

          if ( DV > max_accepted_DV ) continue;

          real_type tt1 = X[0];
          real_type tt2 = X[1];

          dvec3_t  P1, V1, W1, P2, V2, W2;
          a_from.position( tt1, P1 );
          a_from.velocity( tt1, V1 );
          a_to.position( tt2, P2 );
          a_to.velocity( tt2, V2 );

          integer okk = Lambert( P1, P2, tt2-tt1, 0, muS, W1, W2 );
          trips.emplace_back( tt1, P1, V1, W1, tt2, P2, V2, W2 );
        }
      }
    }
  }

  //!
  //!  Given two plates/asteroids find the best possible theoretical DV.
  //!
  //!  \param[IN]  a_from        starting planet/asteroid
  //!  \param[IN]  a_to          arrival planet/asteroid
  //!  \param[IN]  day_tolerance tolerance for pattern search
  //!  \param[IN]  console       pointer to console for printing
  //!  \return  the minimal DV
  //!
  void
  global_minimum_DeltaV2(
    Astro const        & a_from,
    Astro const        & a_to,
    vector<DV_collect> & v_DV,
    Utils::Console     * console
  ) {

    real_type L_tolerance = m_pi/360;
    real_type L_equal     = 5*L_tolerance;
    integer   HJ_max_iter = 200;

    v_DV.clear();
    v_DV.reserve(10);

    real_type L_table[TABLE_SIZE];
    real_type DV_table[TABLE_SIZE][TABLE_SIZE];

    real_type muS = a_from.get_muS();

    real_type DeltaV1, DeltaV2;

    dvec3_t R1[TABLE_SIZE], V1[TABLE_SIZE],
            R2[TABLE_SIZE], V2[TABLE_SIZE],
            W1, W2;

    // prima tabella
    real_type delta_L = m_2pi/TABLE_SIZE;
    for ( integer i = 0 ; i < TABLE_SIZE ; ++i ) {
      L_table[i] = i*delta_L;
      a_from.position_by_L( L_table[i], R1[i] );
      a_from.velocity_by_L( L_table[i], V1[i] );
      a_to.position_by_L( L_table[i], R2[i] );
      a_to.velocity_by_L( L_table[i], V2[i] );
    }

    // trovato minimo locale, raffino con Nelder Mead
    real_type L_from_min, L_from_max, L_to_min, L_to_max;
    auto fun = [ &a_from, &a_to, &L_from_min, &L_from_max, &L_to_min, &L_to_max, muS ] ( real_type const X[] )->real_type {
      real_type La = X[0];
      real_type Lb = X[1];
      if ( La < L_from_min || La > L_from_max || Lb < L_to_min || Lb > L_to_max ) Utils::Inf<real_type>();
      dvec3_t P1, V1, P2, V2;
      real_type DV1, DV2;
      a_from.position_by_L( La, P1 ); a_from.velocity_by_L( La, V1 );
      a_to.position_by_L( Lb, P2 );   a_to.velocity_by_L( Lb, V2 );
      return minimum_DeltaV2( muS, P1, V1, P2, V2, DV1, DV2 );
    };
    std::function<real_type(real_type const[])> F(fun);

    // riempio tabella 8*8
    for ( integer i = 0; i < TABLE_SIZE; ++i ) {
      for ( integer j = 0; j < TABLE_SIZE; ++j ) {
        DV_table[i][j] = minimum_DeltaV2( muS, R1[i], V1[i], R2[j], V2[j], DeltaV1, DeltaV2 );
      }
    }

    // cerco minimi locali
    for ( integer i = 0; i < TABLE_SIZE; ++i ) {

      integer im1{i-1}, ip1{i+1};
      if      ( im1 < 0           ) im1 = TABLE_SIZE-1;
      else if ( ip1 >= TABLE_SIZE ) ip1 = 0;

      integer ii[] = { im1, i, ip1,
                       im1,    ip1,
                       im1, i, ip1 };

      real_type La = L_table[i];
      L_from_min = La-2*delta_L;
      L_from_max = La+2*delta_L;

      for ( integer j = 0 ; j < TABLE_SIZE ; ++j ) {
        real_type DV = DV_table[i][j] ;

        // calcolo primi vicini
        integer jm1{j-1}, jp1{j+1};
        if      ( jm1 < 0           ) jm1 = TABLE_SIZE-1;
        else if ( jp1 >= TABLE_SIZE ) jp1 = 0;

        integer jj[] = { jm1, jm1, jm1,
                         j,        j,
                         jp1, jp1, jp1 };

        // controllo se minimo locale
        bool ok_min = true;
        for ( integer k=0 ; k < 8 && ok_min ; ++k )
          ok_min = DV <= DV_table[ ii[k] ][ jj[k] ];

        if ( !ok_min ) continue ;

        // trovato minimo locale, raffino con Nelder Mead
        real_type Lb = L_table[j];
        L_to_min = Lb-2*delta_L;
        L_to_max = Lb+2*delta_L;

        //Utils::NelderMead<real_type> solver("NMsolver");
        Utils::HJPatternSearch<real_type> solver("HJPatternSearch");
        solver.setup( 2, F, console );
        solver.set_tolerance( L_tolerance );
        solver.set_max_iterations(HJ_max_iter);
        real_type X[2] = {La,Lb};
        solver.run( X, delta_L );
        DV = solver.get_last_solution( X );
        // discard probaby double solutions
        if ( abs(La-X[0]) > 1.5*delta_L || abs(Lb-X[1]) > 1.5*delta_L ) continue;
        v_DV.emplace_back(DV,X[0],X[1]);
      }
    }

    std::sort(
      v_DV.begin(), v_DV.end(),
      [L_equal]( DV_collect const & A, DV_collect const & B ) -> bool {
        real_type dL = B.L_from-A.L_from;
        if ( abs(dL) < L_equal ) return B.L_to-A.L_to;
        else                     return dL > 0;
      }
    );

    // remove duplicates
    v_DV.erase(
      std::unique(
        v_DV.begin(), v_DV.end(),
        [L_equal]( DV_collect const & A, DV_collect const & B ) -> bool {
          real_type dL1 = B.L_from-A.L_from;
          real_type dL2 = B.L_to-A.L_to;
          return abs(dL1) <= L_equal && abs(dL2) <= L_equal;
        }
      ),
      v_DV.end()
    );
  }

  // ==========================================================================
  // ==========================================================================
  // ==========================================================================

  void
  global_minimum_DeltaV2(
    Astro const                 & a_from,
    Astro const                 & a_to,
    real_type                     t_begin,
    real_type                     t_end,
    real_type                     t_tolerance,
    real_type                     max_accepted_DV,
    vector<minimum_DeltaV_trip> & trips,
    Utils::Console              * console
  ) {

    using std::min;

    real_type muS = a_from.get_muS();
    integer   HJ_max_iter = 200;

    vector<DV_collect> v_DV;
    global_minimum_DeltaV2( a_from, a_to, v_DV, console );
    // loop sui candidati
    real_type from_period = a_from.period();
    real_type to_period   = a_to.period();
    trips.clear();
    trips.reserve( 4*v_DV.size());

    real_type t_from_min, t_from_max, t_to_min, t_to_max;
    auto fun = [ &a_from, &a_to, &t_from_min, &t_from_max, &t_to_min, &t_to_max, muS ] ( real_type const X[] )->real_type {
      real_type ta = X[0];
      real_type tb = X[1];
      if ( ta < t_from_min || ta > t_from_max || tb < t_to_min || tb > t_to_max ) Utils::Inf<real_type>();
      dvec3_t P1, V1, VV1, P2, V2, VV2;
      a_from.position( ta, P1 );  a_from.velocity( ta, V1 );
      a_to.position( tb, P2 );    a_to.velocity( tb, V2 );
      integer okk = Lambert( P1, P2, tb-ta, 0, muS, VV1, VV2 );
      if ( okk != 1 ) return Utils::Inf<real_type>();
      return (V1-VV1).norm()+(V2-VV2).norm();
    };
    std::function<real_type(real_type const[])> F(fun);

    real_type DeltaT = min(from_period,to_period)/16;

    for ( DV_collect const & S : v_DV ) {

      if ( S.DV > max_accepted_DV ) continue; // skip

      real_type a_t0 = a_from.time_from_L_angle( t_begin, S.L_from );
      real_type b_t0 = a_to.time_from_L_angle( t_begin, S.L_to );

      for ( real_type ta = a_t0; ta < t_end; ta += from_period ) {

        t_from_min = ta - DeltaT;
        t_from_max = ta + DeltaT;

        for ( real_type tb = b_t0; tb < t_end; tb += to_period ) {
          if ( tb <= ta ) continue;

          t_to_min = tb - DeltaT;
          t_to_max = tb + DeltaT;

          Utils::HJPatternSearch<real_type> solver("HJPatternSearch");
          solver.setup( 2, F, console );
          solver.set_tolerance( t_tolerance );
          solver.set_max_iterations(HJ_max_iter);
          real_type X[2] = {ta,tb};
          solver.run( X, DeltaT/2 );
          real_type DV = solver.get_last_solution( X );

          if ( DV > max_accepted_DV ) continue;

          real_type tt1 = X[0];
          real_type tt2 = X[1];

          dvec3_t  P1, V1, W1, P2, V2, W2;
          a_from.position( tt1, P1 );
          a_from.velocity( tt1, V1 );
          a_to.position( tt2, P2 );
          a_to.velocity( tt2, V2 );

          integer okk = Lambert( P1, P2, tt2-tt1, 0, muS, W1, W2 );
          trips.emplace_back( tt1, P1, V1, W1, tt2, P2, V2, W2 );
        }
      }
    }
  }

  // ---------------------------------------------------------------------------

  //!
  //!  Given two plates/asteroids find the best trips in term of DV
  //!  in a specific time range
  //!
  //!  \param[IN]  who -1  minimize initial DV1 +1 minimize final DV2 0 minimize DV1+DV2
  //!  \param[IN]  t_begin initial time
  //!  \param[IN]  t_end   final time
  //!  \param[IN]  delta_t initial granularity of the search
  //!  \param[IN]  a_from  starting planet/asteroid
  //!  \param[IN]  a_to    arrival planet/asteroid
  //!  \param[OUT] trips   list of the minimal DV trips
  //!  \param[IN]  max_accepted_DV
  //!  \param[IN]  day_tolerance
  //!  \param[IN]  max_subiter
  //!
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
  ) {

    trips.clear();
    trips.reserve(10);

    integer N_TABLE = ceil((t_end-t_begin)/delta_t);

    dvec_t t_table;
    dmat_t DV_table;

    t_table.resize(N_TABLE+1);
    DV_table.resize(N_TABLE+1,N_TABLE+1);

    std::vector<dvec3_t> R1_vec(N_TABLE+1),
                         V1_vec(N_TABLE+1),
                         R2_vec(N_TABLE+1),
                         V2_vec(N_TABLE+1);

    real_type dt = (t_end-t_begin)/N_TABLE;
    for ( integer i = 0; i <= N_TABLE; ++i ) {
      real_type t = t_begin+i*dt;
      t_table.coeffRef(i) = t;
      a_from.position( t, R1_vec[i] );
      a_from.velocity( t, V1_vec[i] );
      a_to.position( t, R2_vec[i] );
      a_to.velocity( t, V2_vec[i] );
    }

    real_type muS = a_from.get_muS();

    auto fun = [ &a_from, &a_to, &t_begin, &t_end, muS, who ]( real_type const X[] )->real_type {
      real_type DV = Utils::Inf<real_type>();
      real_type t1 = X[0];
      real_type t2 = X[1];
      if ( ! ( t_begin <= t1 && t1 < t2 && t2 <= t_end ) ) return DV;
      dvec3_t P1, V1, P2, V2, VV1, VV2;
      a_from.position( t1, P1 );
      a_from.velocity( t1, V1 );
      a_to.position( t2, P2 );
      a_to.velocity( t2, V2 );
      integer okk = Lambert( P1, P2, t2-t1, 0, muS, VV1, VV2 );
      if ( okk == 1 ) {
        real_type DV1 = (V1-VV1).norm();
        real_type DV2 = (V2-VV2).norm();
        switch ( who ) {
          case -1: DV = DV1;     break;
          case  1: DV = DV2;     break;
          default: DV = DV1+DV2; break;
        }
      }
      return DV;
    };

    std::function<real_type(real_type const[])> F(fun);

    // riempio tabella 8*8
    DV_table.fill( Utils::Inf<real_type>() );
    real_type X[2];
    for ( integer i = 0; i < N_TABLE; ++i ) {
      X[0] = t_table.coeff(i);
      for ( integer j = i+1; j <= N_TABLE; ++j ) {
        X[1] = t_table.coeff(j);
        DV_table(i,j) = fun( X );
        //fmt::print( "DV_table({},{}) = {}, DT={}\n", i, j, DV_table(i,j), X[1]-X[0] );
      }
    }

    // cerco minimi locali
    for ( integer i = 0; i < N_TABLE; ++i ) {

      // trovato minimo locale, raffino con Nelder Mead
      real_type t1 = t_table.coeff(i);

      // calcolo primi vicini
      integer im1{i-1}, ip1{i+1};
      if      ( im1 < 0       ) im1 = N_TABLE;
      else if ( ip1 > N_TABLE ) ip1 = 0;

      integer ii[] = { im1, i, ip1,
                       im1,    ip1,
                       im1, i, ip1 };

      for ( integer j = i+1; j <= N_TABLE; ++j ) {
        real_type DV = DV_table.coeff(i,j);

        // calcolo primi vicini
        integer jm1{j-1}, jp1{j+1};
        if      ( jm1 < 0       ) jm1 = N_TABLE;
        else if ( jp1 > N_TABLE ) jp1 = 0;

        integer jj[] = { jm1, jm1, jm1,
                         j,        j,
                         jp1, jp1, jp1 };

        // controllo se minimo locale
        bool ok_min = true;
        for ( integer k=0; k < 8 && ok_min; ++k )
          ok_min = DV <= DV_table.coeff( ii[k], jj[k] );

        if ( !ok_min ) continue ;

        // trovato minimo locale, raffino con Nelder Mead
        real_type t2 = t_table.coeff(j);

        Utils::Console console(&std::cout,0);
        //Utils::NelderMead<real_type> solver("NMsolver");
        Utils::HJPatternSearch<real_type> solver("HJPatternSearch");

        solver.setup( 2, F, &console );
        solver.set_tolerance( day_tolerance );
        solver.set_max_iterations( max_subiter );
        real_type X[2] = {t1,t2};
        solver.run( X, dt );
        DV = solver.get_last_solution( X );

        if ( DV > max_accepted_DV ) continue;

        real_type tt1 = X[0];
        real_type tt2 = X[1];

        //fmt::print( "X={}, Y={}, DV={}\n", tt1, tt2, DV );

        // confronto con minimo
        dvec3_t  P1, V1, W1, P2, V2, W2;
        //trip.t_begin = tt1;
        //trip.t_end   = tt2;

        a_from.position( tt1, P1 );
        a_from.velocity( tt1, V1 );
        a_to.position( tt2, P2 );
        a_to.velocity( tt2, V2 );

        integer okk = Lambert( P1, P2, tt2-tt1, 0, muS, W1, W2 );
        //trip.DeltaV1 = (V1-trip.W1).norm();
        //trip.DeltaV2 = (V2-trip.W2).norm();

        trips.emplace_back( tt1, P1, V1, W1, tt2, P2, V2, W2 );
      }
    }

    std::sort(
      trips.begin(), trips.end(),
      [day_equal]( minimum_DeltaV_trip const & A, minimum_DeltaV_trip const & B ) -> bool {
        real_type dt = B.t_end-A.t_end;
        if ( abs(dt) < day_equal ) return B.t_begin-A.t_begin;
        else                       return dt > 0;
      }
    );

    // remove duplicates
    trips.erase(
      std::unique(
        trips.begin(), trips.end(),
        [day_equal]( minimum_DeltaV_trip const & A, minimum_DeltaV_trip const & B ) -> bool {
          real_type dt1 = B.t_end-A.t_end;
          real_type dt2 = B.t_begin-A.t_begin;
          return abs(dt1) <= day_equal && abs(dt2) <= day_equal;
        }
      ),
      trips.end()
    );
  }

}
