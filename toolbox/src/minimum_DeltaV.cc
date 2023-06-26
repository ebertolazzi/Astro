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

  using std::acos;
  using std::abs;
  using std::min;
  using std::max;
  using std::atan2;

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

    // cerco soluzione
    Utils::Poly<real_type>::dvec_t coeffs8;
    coeffs8.resize(9); // ordine 9
    coeffs8 << poly0, poly1, poly2, poly3, poly4, poly5, poly6, poly7, poly8;
    for ( int i=0; i < 9;++i ) if ( isnan(coeffs8.coeff(i)) || isinf(coeffs8.coeff(i)) ) return Utils::Inf<real_type>();

    Utils::Poly<real_type>  poly(coeffs8);
    Utils::Sturm<real_type> sturm;

    poly.normalize();
    sturm.build( poly );

    integer n_roots = sturm.separate_roots();
    sturm.refine_roots();

    Utils::Poly<real_type>::dvec_t roots = sturm.roots();

    //fmt::print( "N.roots = {}\nCheck\n", n_roots );
    //for ( auto & x : roots) {
    //  fmt::print( "P({}) = {}\n", x, poly.eval(x) );
    //}

    real_type minDV = Utils::Inf<real_type>();
    for ( int k = 0; k < n_roots; ++k ) {

      real_type omega__c = roots.coeff(k);

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
   |    A closed-form solution to the minimum deltaV^2 problem Lambert's problem
   |    Martín Avendaño, Daniele Mortari
   |    2009
   |
   |    Celest Mech Dyn Astr (2010) 106:25-37
   |    DOI 10.1007/s10569-009-9238-x
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

    real_type V__1   = V1.norm();
    real_type V__1r  = V1.dot(rr1);
    real_type V__1c  = V1.dot(ch);

    real_type V__2   = V2.norm();
    real_type V__2r  = V2.dot(rr2);
    real_type V__2c  = V2.dot(ch);

    real_type K           = mu*c/(r1*r2+R1.dot(R2)); // equation (7)
    real_type cos__phi__1 = ch.dot(rr1);
    real_type cos__phi__2 = ch.dot(rr2);

    real_type poly0 = -4*K*K;
    real_type poly1 = 2*K*(V__1r-V__2r);
    real_type poly2 = 0;
    real_type poly3 = -2*(V__1c+V__2c);
    real_type poly4 = 4;

    // cerco soluzione
    Utils::Poly<real_type>::dvec_t coeffs4;
    coeffs4.resize(5); // ordine 5
    coeffs4 << poly0, poly1, poly2, poly3, poly4;
    for ( int i=0; i <= 4;++i ) if ( isnan(coeffs4.coeff(i)) || isinf(coeffs4.coeff(i)) ) return Utils::Inf<real_type>();

    Utils::Poly<real_type>  poly(coeffs4);
    Utils::Sturm<real_type> sturm;

    poly.normalize();
    sturm.build( poly );

    integer n_roots = sturm.separate_roots();
    sturm.refine_roots();

    Utils::Poly<real_type>::dvec_t roots = sturm.roots();

    //fmt::print( "N.roots = {}\nCheck\n", n_roots );
    //for ( auto & x : roots) {
    //  fmt::print( "P({}) = {}\n", x, poly.eval(x) );
    //}

    real_type minDV2 = Utils::Inf<real_type>();
    for ( int k = 0; k < n_roots; ++k ) {
      real_type omega__c  = roots.coeff(k);

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

  //#define SUB_TABLE_SIZE 64
  //#define TABLE_SIZE     12
  //#define TABLE_SIZE     18

  bool
  minimum_DeltaV(
    real_type     muC,
    real_type     t_begin,
    real_type     t_end,
    real_type     t_delta,
    Astro const & a_from,
    Astro const & a_to,
    minimum_DeltaV_trip & trip
  ) {

    bool ok = false;

    dvec3_t Pa, Va, Pb, Vb, VVa, VVb;
    real_type minDV = Utils::Inf<real_type>();

    for ( real_type ta = t_begin; ta < t_end; ta += t_delta ) {
      a_from.position( ta, Pa );
      a_from.velocity( ta, Va );
      for ( real_type tb = ta+t_delta; tb < t_end; tb += t_delta ) {
        a_to.position( tb, Pb );
        a_to.velocity( tb, Vb );
        // cerco soluzione con Lambert
        // real_type tf = trip.optimalTravelTime;
        // int       m  = 0;
        // if ( trip.long_path   );
        // if ( trip.left_branch ) tf = -tf;
        integer   m    = 0;     // numero rivoluzioni
        real_type tfly = tb-ta; // tempo di volo
        integer   okk  = Lambert( Pa, Pb, tfly, m, muC, VVa, VVb );
        if ( okk == 1 ) {
          real_type DV1 = (Va-VVa).norm();
          real_type DV2 = (Vb-VVb).norm();
          real_type DV  = DV1+DV2;
          if ( DV < minDV ) {
            ok           = true;
            minDV        = DV;
            trip.t_begin = ta;
            trip.t_end   = tb;
            trip.DeltaV1 = DV1;
            trip.DeltaV2 = DV2;
          }
        }
      }
    }
    return ok;
  }

}
