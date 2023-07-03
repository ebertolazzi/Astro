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

using AstroLib::real_type;
using AstroLib::integer;
using Utils::m_pi;

real_type const MJD_begin = 64328; // January 1, 2035, 00:00:00
real_type const MJD_end   = 69807; // January 1, 2050, 00:00:00

int
main() {

  real_type AU_to_km = 1.49597870691E8;
  real_type day_to_s = 86400;
  real_type mu_SUN   = 1.32712440018E11;
  mu_SUN /= (AU_to_km*AU_to_km*AU_to_km)/(day_to_s*day_to_s);

  AstroLib::Astro E, M;

  real_type t0        = 54000;
  real_type a         = 1.00000261;
  real_type e         = 0.01671123;
  real_type i         = -0.00001531*(m_pi/180);
  real_type L         = 100.46457166*(m_pi/180);
  real_type omega_bar = 102.93768193*(m_pi/180);
  real_type Omega     = 0*(m_pi/180);
  real_type omega     = omega_bar - Omega;
  real_type M0        = L - omega_bar;

  E.setup_Keplerian( "Terra", t0, a, e, Omega, omega, i, M0, mu_SUN );

  t0        = 54000;
  a         = 1.52371034;
  e         = 0.09339410;
  i         = 1.84969142*(m_pi/180);
  L         = -4.55343205*(m_pi/180);
  omega_bar = -23.94362959*(m_pi/180);
  Omega     = 49.55953891*(m_pi/180);
  omega     = omega_bar - Omega;
  M0        = L - omega_bar;

  M.setup_Keplerian( "Mars", t0, a, e, Omega, omega, i, M0, mu_SUN );

  //std::cout << E.info() << "\n";
  //std::cout << M.info() << "\n";

  AstroLib::dvec3_t P0, V0, P1, V1;
  real_type         DV0, DV1;

  #if 0
  for ( real_type tt = 0; tt <= 20; tt += 10 ) {

    E.position( t0+tt, P0.data() );
    E.velocity( t0+tt, V0.data() );
    M.position( t0+tt, P1.data() );
    M.velocity( t0+tt, V1.data() );

    real_type DV = AstroLib::minimum_DeltaV( mu_SUN, P0, V0, P1, V1, DV0, DV1, nullptr );
    fmt::print(
      "DV1+DV2     = {:.4}\n"
      "DV1         = {:.4}\n"
      "DV2         = {:.4}\n\n",
      DV, DV0, DV1
    );

    real_type DV2 = AstroLib::minimum_DeltaV2( mu_SUN, P0, V0, P1, V1, DV0, DV1 );
    fmt::print(
      "DV1^2+DV2^2 = {:.4}\n"
      "(UA/Y) DV1  = {:.4}\n"
      "(UA/Y) DV2  = {:.4}\n",
      "(km/s) DV1  = {:.4}\n"
      "(km/s) DV2  = {:.4}\n\n",
      DV2, DV0, DV1, DV0*(AU_to_km/Y_to_s), DV1*(AU_to_km/Y_to_s)
    );
  }
  #endif

  #if 1

  std::vector<AstroLib::minimum_DeltaV_trip> trips;
  real_type maxDV = 6/((AU_to_km/day_to_s));
  AstroLib::minimum_DeltaV( -1, mu_SUN, MJD_begin, MJD_end, 10, E, M, trips, maxDV );
  integer it = 0;
  for ( auto t : trips ) {
    fmt::print("\n\n\n{}\n\n\n",++it);
    fmt::print(
      "(Year)   ta    = {:.4}\n"
      "(Year)   tb    = {:.4}\n"
      "(Day)    tb-ta = {:.4}\n"
      "(UA/day) DV1   = {:.4}\n"
      "(UA/day) DV2   = {:.4}\n"
      "(km/s)   DV1   = {:.4}\n"
      "(km/s)   DV2   = {:.4}\n"
      "(UA/day) DV1+DV2 = {:.4}\n"
      "(km/s)   DV1+DV2 = {:.4}\n\n",
      (t.t_begin-MJD_begin)/365.25, (t.t_end-MJD_begin)/365.25,
      t.t_end-t.t_begin,
      t.DeltaV1, t.DeltaV2,
      t.DeltaV1*(AU_to_km/day_to_s),
      t.DeltaV2*(AU_to_km/day_to_s),
      t.DeltaV1+t.DeltaV2,
      (t.DeltaV1+t.DeltaV2)*(AU_to_km/day_to_s)
    );
  }

  #endif

  #if 1
  Utils::Console console(&std::cout,-1);
  real_type DV = minimum_DeltaV( E, M, 0.01, &console );
  fmt::print(
    "(UA/Y) DV = {:.8}\n"
    "(km/s) DV = {:.8}\n\n",
    DV, DV*(AU_to_km/day_to_s)
  );
  #endif

  std::cout << "All done folks!!\n";
  return 0;
}
