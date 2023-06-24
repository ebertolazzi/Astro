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

  real_type
  minimum_DeltaV(
    real_type       mu,
    dvec3_t const & R0,
    dvec3_t const & V0,
    dvec3_t const & R1,
    dvec3_t const & V1,
    minimum_DeltaV_extra * extra
  ) {

    dvec3_t ch = R1-R0;
    real_type c = ch.norm();
    ch /= c;
    real_type r0   = R0.norm();
    real_type r1   = R1.norm();
    real_type r0r1 = R0.dot(R1);
    real_type K    = mu*c/(r0*r1+r0r1);
    // -----------------------------------
    real_type v0   = V0.norm();
    real_type v1   = V1.norm();
    real_type cr0  = ch.dot(R0)/r0;
    real_type cr1  = ch.dot(R1)/r1;
    real_type cv0  = ch.dot(V0);
    real_type cv1  = ch.dot(V1);
    real_type v0r0 = V0.dot(R0)/r0;
    real_type v1r1 = V1.dot(R1)/r1;

    // coefficienti del polinomio le cui radici sono candidati minimi
    real_type t1  = cr0 + cr1;
    real_type t2  = v1*v1;
    real_type t3  = v0*v0;
    real_type t4  = cv1*cv1;
    real_type t5  = cv0*cv0;
    real_type t6  = 2;
    real_type t7  = t6 * t1 * K;
    real_type t8  = t7 - t5 + t3 + t4 - t2;
    real_type t9  = cr0 * cv1;
    real_type t10 = (cr1 * cv0);
    real_type t11 = cv0 * cv1;
    real_type t12 = t2 * cv0;
    real_type t13 = 4;
    real_type t14 = t13 * (v0r0 + v1r1 + t10 + t9) * K;
    real_type t15 = -t6 * ((t11 - t5 + t3) * cv1 - t12) - t14;
    real_type t16 = t6 * v0r0;
    real_type t17 = 8 * v1r1;
    real_type t18 = t6 * v1r1;
    real_type t19 = 8 * v0r0;
    real_type t20 = (t6 * (cr0 * t4 + cr1 * t5) + (t17 + t16) * cv0 + (t19 + t18) * cv1) * K - t2 * t5 + t3 * t4;
    real_type t21 = t13 * (-cr0 * v1r1 + cr1 * v0r0 + cv0 - cv1) * K;
              t11 = K * (t21 - t6 * ((t4 + t2) * v0r0 + (t5 + t3) * v1r1) - t11 * t13 * (v0r0 + v1r1));
    real_type t22 = v0r0*v0r0;
    real_type t23 = v1r1*v1r1;
    real_type t24 = t3 * cv1;
              t1  = ((-t1 * t13 * K - t13 * (t10 * v0r0 - t9 * v1r1) - t6 * (t3 - t2) - t22 + t23 + t4 - t5) * K
                     + t6 * (t12 * v0r0 + t24 * v1r1)) * K;
              t4   = K*K;
              t5   = t4 * (t14 + t6 * (t24 - t12) + v1r1 * (t13 * v0r0 - t18) * cv0 + v0r0 * (-t13 * v1r1 + t16) * cv1);
              t9   = t4 * ((t6 * (cr0 * t23 + cr1 * t22) + (-t17 + t16) * cv0 + (-t19 + t18) * cv1) * K - t2 * t22 + t23 * t3);
              t6   = K * t4 * (t6 * (t2 * v0r0 + (-v0r0 * v1r1 - t22 + t3) * v1r1) - t21);
              t2   = (t4*t4) * (t7 + t23 + t3 - t22 - t2);

    // cerco soluzione
    Utils::Poly<real_type>::dvec_t coeffs8;
    coeffs8.resize(9); // ordine 9
    coeffs8 << t8, t15, t20, t11, t1, t5, t9, t6, t2;
    for ( int i=0; i < 9;++i ) if ( isnan(coeffs8.coeff(i)) || isinf(coeffs8.coeff(i)) ) return 1e10;

    Utils::Poly<real_type> poly(coeffs8);
    Utils::Sturm<real_type> sturm;

    poly.normalize();
    sturm.build( poly );

    // Cauchy's bounds for roots
    real_type bnd     = 1+coeffs8.cwiseAbs().maxCoeff()/abs(coeffs8.coeff(0));
    integer   n_roots = sturm.separate_roots( -bnd, bnd );
    sturm.refine_roots();
    Utils::Poly<real_type>::dvec_t roots = sturm.roots();

    real_type minDV = 1e100;
    for ( int k = 0; k < n_roots; ++k ) {
      real_type omegac  = roots.coeff(k);
      real_type omegac2 = omegac*omegac;
      real_type Kc      = K/omegac;

      // valore
      real_type delta0 = sqrt(v0*v0 + 2*K*cr0 + omegac*(omegac-2*cv0)+(Kc)*(Kc-2*v0r0));
      real_type delta1 = sqrt(v1*v1 - 2*K*cr1 + omegac*(omegac-2*cv1)+(Kc)*(Kc+2*v1r1));
      real_type DV     = delta0+delta1;
      // derivata prima
      real_type d2_0 = 2*(omegac-cv0+K*(v0r0-Kc)/omegac2);
      real_type d2_1 = 2*(omegac-cv1-K*(v1r1+Kc)/omegac2);
      real_type d    = 0.5*(d2_0/delta0 + d2_1/delta1);

      if ( abs(d) < 1e-9 ) {
        // derivata seconda
        real_type dd2_0 = 2+Kc*(6*Kc-4*v0r0)/omegac2;
        real_type dd2_1 = 2+Kc*(6*Kc+4*v1r1)/omegac2;
        real_type dd    = (dd2_0 - d2_0*d2_0/(2*delta0*delta0))/(2*delta0)
                        + (dd2_1 - d2_1*d2_1/(2*delta1*delta1))/(2*delta1);
        if ( dd > 0 ) {
          // trovato minimo
          if ( DV < minDV ) {
            minDV = DV;
            if ( extra != nullptr ) {
              bool      & long_path   = extra->long_path;
              bool      & left_branch = extra->left_branch;
              real_type & DeltaV0     = extra->DeltaV0;
              real_type & DeltaV1     = extra->DeltaV1;
              real_type & ott         = extra->optimal_travel_time;
              real_type & period      = extra->period;
              dvec3_t   & W0          = extra->W0;
              dvec3_t   & W1          = extra->W1;

              DeltaV0 = delta0;
              DeltaV1 = delta1;
              W0      = omegac*ch+(Kc/r0)*R0;
              W1      = omegac*ch-(Kc/r1)*R1;
              // calcolo alpha e beta
              real_type a     = mu*r0/(2*mu-r0*W0.squaredNorm());
              real_type s     = (r0+r1+c)/2;
              real_type alpha = acos(1-s/a);
              real_type beta  = acos(1-(s-c)/a);
              // Optimal two-impulse rendezvous using constrained multiple-revolution Lambert solutions
              // Gang Zhang - Di Zhou -  Daniele Mortari
              // per determinare long/short part calcolo p (eq. 6)
              dvec3_t r0_x_r1 = R0.cross(R1);
              real_type p = (omegac/c);
              p = r0_x_r1.squaredNorm()*p*p/mu;

              // calcolo pm da equazione (21)
              real_type pm = 2*(s-r0)*(s-r1)/c;
              long_path   = p <= pm;
              left_branch = omegac < 0; // da equazione (7) segno(omegac)=sign(sin(theta))
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
        } else {
          // Gradiente zero ma non minimo
        }
      } else {
        // Gradiente NON zero
      }
    }
    return minDV;
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
