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

  #if 0

  // ---------------------------------------------------------------------------

  #define SUB_TABLE_SIZE 64
  //#define TABLE_SIZE     12
  #define TABLE_SIZE     18

  void
  globalFreeTimeMinimumDeltaV( real_type const   tBegin,
                               int               fromAst,
                               Asteroid  const & from,
                               int               toAst,
                               Asteroid  const & to,
                               vector<LambertFreeTimeTrip> & trips ) {

    real_type E_table[TABLE_SIZE];

    real_type L_from_table[TABLE_SIZE];
    real_type L_to_table[TABLE_SIZE];
    real_type DV_table[TABLE_SIZE][TABLE_SIZE];

    real_type L_from_sub_table[2*SUB_TABLE_SIZE+1];
    real_type L_to_sub_table[2*SUB_TABLE_SIZE+1];
    real_type DV_sub_table[2*SUB_TABLE_SIZE+1][2*SUB_TABLE_SIZE+1];

    bool      long_path, left_branch;
    real_type DeltaV0, DeltaV1, travelTime, period;

    real_type R0[TABLE_SIZE][3], V0[TABLE_SIZE][3],
              R1[TABLE_SIZE][3], V1[TABLE_SIZE][3],
              R0_sub[2*SUB_TABLE_SIZE+1][3], V0_sub[2*SUB_TABLE_SIZE+1][3],
              R1_sub[2*SUB_TABLE_SIZE+1][3], V1_sub[2*SUB_TABLE_SIZE+1][3],
              W0[3], W1[3];

    real_type dE  = (2*m_pi)/TABLE_SIZE;
    real_type ddE = dE/SUB_TABLE_SIZE;
    for ( indexType i = 0; i < TABLE_SIZE; ++i ) {
      E_table[i]      = i*dE;
      L_from_table[i] = from.LfromTrueAnomaly( eccentricAnomalyToTrueAnomaly( E_table[i], from.eOrbital() ) );
      L_to_table[i]   = to.LfromTrueAnomaly( eccentricAnomalyToTrueAnomaly( E_table[i], to.eOrbital() ) );

      from.position_by_L( L_from_table[i], R0[i][0], R0[i][1], R0[i][2] ); // UA
      from.velocity_by_L( L_from_table[i], V0[i][0], V0[i][1], V0[i][2] ); // UA/day

      to.position_by_L( L_to_table[i], R1[i][0], R1[i][1], R1[i][2] ); // UA
      to.velocity_by_L( L_to_table[i], V1[i][0], V1[i][1], V1[i][2] ); // UA/day
    }

    // riempio tabella 8*8
    for ( indexType i = 0; i < TABLE_SIZE; ++i ) {
      for ( indexType j = 0; j < TABLE_SIZE; ++j ) {
        DV_table[i][j] = freeTimeMinimumDeltaV( R0[i], V0[i],
                                                R1[j], V1[j],
                                                false,
                                                long_path,
                                                left_branch,
                                                DeltaV0,
                                                DeltaV1,
                                                travelTime,
                                                period,
                                                W0, W1 );
      }
    }

    // cerco minimi locali
    trips.clear();
    real_type DV_min = 1e10;
    for ( indexType i = 0; i < TABLE_SIZE; ++i ) {
      indexType im1 = (i+TABLE_SIZE-1)%TABLE_SIZE;
      indexType ip1 = (i+1)%TABLE_SIZE;

      for ( indexType j = 0; j < TABLE_SIZE; ++j ) {
        indexType jm1 = (j+TABLE_SIZE-1)%TABLE_SIZE;
        indexType jp1 = (j+1)%TABLE_SIZE;
        real_type DV  = DV_table[i][j];

        if ( DV > DV_table[ip1][jm1] ||
             DV > DV_table[ip1][j]   ||
             DV > DV_table[ip1][jp1] ||
             DV > DV_table[im1][jm1] ||
             DV > DV_table[im1][j]   ||
             DV > DV_table[im1][jp1] ||
             DV > DV_table[i][jm1]   ||
             DV > DV_table[i][jp1] ) continue;

        // trovato minimo locale, raffino
        std::fill( DV_sub_table[0], DV_sub_table[0] + (2*SUB_TABLE_SIZE+1)*(2*SUB_TABLE_SIZE+1), -1.0 );

        DV_sub_table[0*SUB_TABLE_SIZE][0*SUB_TABLE_SIZE] = DV_table[im1][jm1];
        DV_sub_table[1*SUB_TABLE_SIZE][0*SUB_TABLE_SIZE] = DV_table[i  ][jm1];
        DV_sub_table[2*SUB_TABLE_SIZE][0*SUB_TABLE_SIZE] = DV_table[ip1][jm1];

        DV_sub_table[0*SUB_TABLE_SIZE][1*SUB_TABLE_SIZE] = DV_table[im1][j];
        DV_sub_table[1*SUB_TABLE_SIZE][1*SUB_TABLE_SIZE] = DV_table[i  ][j];
        DV_sub_table[2*SUB_TABLE_SIZE][1*SUB_TABLE_SIZE] = DV_table[ip1][j];

        DV_sub_table[0*SUB_TABLE_SIZE][2*SUB_TABLE_SIZE] = DV_table[im1][jp1];
        DV_sub_table[1*SUB_TABLE_SIZE][2*SUB_TABLE_SIZE] = DV_table[i  ][jp1];
        DV_sub_table[2*SUB_TABLE_SIZE][2*SUB_TABLE_SIZE] = DV_table[ip1][jp1];

        for ( indexType ii = -SUB_TABLE_SIZE; ii <= SUB_TABLE_SIZE; ++ii ) {
          indexType idx = ii+SUB_TABLE_SIZE;
          real_type L_from = from.LfromTrueAnomaly( eccentricAnomalyToTrueAnomaly( E_table[i]+ddE*ii, from.eOrbital() ) );
          real_type L_to   = to.LfromTrueAnomaly( eccentricAnomalyToTrueAnomaly( E_table[j]+ddE*ii, to.eOrbital() ) );

          L_from_sub_table[idx] = L_from;
          L_to_sub_table[idx]   = L_to;

          from.position_by_L( L_from, R0_sub[idx][0], R0_sub[idx][1], R0_sub[idx][2] ); // UA
          from.velocity_by_L( L_from, V0_sub[idx][0], V0_sub[idx][1], V0_sub[idx][2] ); // UA/day

          to.position_by_L( L_to, R1_sub[idx][0], R1_sub[idx][1], R1_sub[idx][2] ); // UA
          to.velocity_by_L( L_to, V1_sub[idx][0], V1_sub[idx][1], V1_sub[idx][2] ); // UA/day
        }

        indexType i_pivot    = SUB_TABLE_SIZE;
        indexType j_pivot    = SUB_TABLE_SIZE;
        real_type DV_min_loc = DV_sub_table[i_pivot][j_pivot];
        for ( indexType kstep = SUB_TABLE_SIZE; kstep > 0; kstep /= 2 ) {

          // calcolo primi vicini
          indexType ii[] = { i_pivot-kstep, i_pivot,       i_pivot+kstep,
                             i_pivot-kstep,                i_pivot+kstep,
                             i_pivot-kstep, i_pivot,       i_pivot+kstep };
          indexType jj[] = { j_pivot-kstep, j_pivot-kstep, j_pivot+kstep,
                             j_pivot,                      j_pivot,
                             j_pivot+kstep, j_pivot+kstep, j_pivot+kstep };
          int kk_min = -1;
          DV_min_loc = DV_sub_table[i_pivot][j_pivot];
          // cerco minimo
          for ( int kk = 0; kk < 8; ++kk ) {
            indexType & iii = ii[kk];
            indexType & jjj = jj[kk];
            if ( iii < 0 || iii > 2*SUB_TABLE_SIZE ||
                 jjj < 0 || jjj > 2*SUB_TABLE_SIZE ) continue;
            // controllo se gia calcolato
            real_type & DV = DV_sub_table[iii][jjj];
            if ( DV < 0 ) {
              DV = freeTimeMinimumDeltaV( R0_sub[iii], V0_sub[iii],
                                          R1_sub[jjj], V1_sub[jjj],
                                          false,
                                          long_path,
                                          left_branch,
                                          DeltaV0, DeltaV1,
                                          travelTime, period,
                                          W0, W1 );
            }

            // confronto con minimo
            if ( DV < DV_min_loc ) {
              DV_min_loc = DV;
              kk_min     = kk;
            }
          }
          // kk_min contiene posizione minimo, se kk_min = -1 pozizione centrale
          if ( kk_min >= 0 ) {
            i_pivot = ii[kk_min];
            j_pivot = jj[kk_min];
          }
        }
        // fine raffinamento trovato minimo locale
        // lo aggiungo alla lista
        LambertFreeTimeTrip tmp;

        copy3(R0_sub[i_pivot],tmp.P0);
        copy3(V0_sub[i_pivot],tmp.V0);

        copy3(R1_sub[j_pivot],tmp.P1);
        copy3(V1_sub[j_pivot],tmp.V1);

        tmp.fromAst = fromAst;
        tmp.toAst   = toAst;

        tmp.Lfrom = L_from_sub_table[i_pivot];
        tmp.Lto   = L_to_sub_table[j_pivot];

        freeTimeMinimumDeltaV( tmp.P0, tmp.V0,
                               tmp.P1, tmp.V1,
                               true,
                               tmp.long_path,
                               tmp.left_branch,
                               tmp.DeltaV0,
                               tmp.DeltaV1,
                               tmp.optimalTravelTime,
                               tmp.period,
                               tmp.W0,
                               tmp.W1 );

        // controllo che non sia troppo vicino ad un minimo precedente
        bool found_nearby = false;
        for ( int i = 0; i < trips.size() && !found_nearby; ++i ) {
          LambertFreeTimeTrip & itrip = trips[i];
          found_nearby = abs(tmp.Lfrom-itrip.Lfrom) < 0.05*m_pi &&
                         abs(tmp.Lto-itrip.Lto)     < 0.05*m_pi;
          if ( found_nearby ) {
            // trovati minimi vicini, scelgo il piu piccolo o tengo vecchio
            if ( tmp.DeltaV0+tmp.DeltaV1 < itrip.DeltaV0+itrip.DeltaV1 ) {
              copyto( tmp, itrip );
              DV_min = min( DV_min, tmp.DeltaV0+tmp.DeltaV1 );
            }
          }
        }

        if ( !found_nearby ) { // nessun vicino aggiungo alla lista
          if ( DV_min_loc < DV_min ) {
            DV_min = DV_min_loc;
            trips.insert ( trips.begin(), tmp );
          } else {
            trips.push_back( tmp );
          }
        }
      }
    }

    #if 0
    cout << "Lfrom = " << trip.Lfrom << " apo = " << from.LatApoapsis() << " peri = " << from.LatPeriapsis() << '\n';
    cout << "Lto   = " << trip.Lto   << " apo = " << to.LatApoapsis()   << " peri = " << to.LatPeriapsis() << '\n';
    #endif

    #if 0
    real_type tol = 1e-6;
    bool ok;
    real_type WW0[3], WW1[3], PP[3];
    Asteroid ast;

    cout << '\n';
    cout << "long = " << trip.long_path   << '\n';
    cout << "left = " << trip.left_branch << '\n';

    // controllo salti
    real_type dv = dist3(trip.W0,trip.V0);
    ASSERT( abs( dv-trip.DeltaV0 ) < tol,
            "DeltaV0 errato ERR " << dv-trip.DeltaV0 );

    dv = dist3(trip.W1,trip.V1);
    ASSERT( abs( dv-trip.DeltaV1 ) < tol,
            "DeltaV1 errato ERR " << dv-trip.DeltaV1 );

    ok = ast.setupUsingPointAndVelocity( "pippo",
                                         trip.P0,
                                         trip.W0,
                                         muSun_UA3DAY2,
                                         trip.timeFrom );
    ASSERT( ok, "failed to build asteroid" );

    ast.position( trip.timeFrom+trip.optimalTravelTime, PP[0], PP[1], PP[2] );
    ASSERT( dist3(trip.P1,PP) < tol,
            "P1 trip failed, difference = " << dist3(trip.P1,PP) <<
            " optimalTravelTime = " << trip.optimalTravelTime );

    ok = ast.setupUsingPointAndVelocity( "pippo",
                                          trip.P1,
                                          trip.W1,
                                          muSun_UA3DAY2,
                                          trip.timeFrom+trip.optimalTravelTime );

    ASSERT( ok, "failed to build asteroid" );

    ast.position( trip.timeFrom, PP[0], PP[1], PP[2] );
    ASSERT( dist3(trip.P0,PP) < tol,
            "P0 trip failed, difference = " << dist3(trip.P0,PP)  <<
            " optimalTravelTime = " << trip.optimalTravelTime);


    // controllo soluzione con Lambert
    real_type tf = trip.optimalTravelTime;
    int       m  = 0;
    if ( trip.long_path   );
    if ( trip.left_branch ) tf = -tf;
    int okk = lambert( trip.P0,
                       trip.P1,
                       tf,
                       m,
                       muSun_UA3DAY2,
                       WW0, WW1 );
    ASSERT( okk == 1,
            "lambert failed, ok = " << okk <<
            (trip.long_path?" LONG " : " SHORT ") <<
            (trip.left_branch?" LEFT " : " RIGHT ") );

    // controllo che interpola la lambert
    ok = ast.setupUsingPointAndVelocity( "pippo",
                                         trip.P0,
                                         WW0,
                                         muSun_UA3DAY2,
                                         trip.timeFrom );
    ASSERT( ok, "failed to build asteroid" );

    ast.position( trip.timeFrom+trip.optimalTravelTime, PP[0], PP[1], PP[2] );
    ASSERT( dist3(trip.P1,PP) < tol, "P1 lambert failed, difference = " << dist3(trip.P1,PP) );

    ok = ast.setupUsingPointAndVelocity( "pippo",
                                         trip.P1,
                                         WW1,
                                         muSun_UA3DAY2,
                                         trip.timeFrom+trip.optimalTravelTime );
    ASSERT( ok, "failed to build asteroid" );

    ast.position( trip.timeFrom, PP[0], PP[1], PP[2] );
    ASSERT( dist3(trip.P0,PP) < tol, "P0 lambert failed, difference = " << dist3(trip.P0,PP) );


    ASSERT( dist3(trip.W0,WW0) < tol, "W0 failed, difference = " << dist3(trip.W0,WW0) <<
            (trip.long_path?" LONG " : " SHORT ") <<
            (trip.left_branch?" LEFT " : " RIGHT ") << "  tf = " << tf );
    ASSERT( dist3(trip.W1,WW1) < tol, "W1 failed, difference = " << dist3(trip.W1,WW1) <<
            (trip.long_path?" LONG " : " SHORT ") <<
            (trip.left_branch?" LEFT " : " RIGHT ") << "  tf = " << tf );
    #endif
  }

  #endif

}
