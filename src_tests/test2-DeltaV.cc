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
using Utils::m_pi;

int
main() {

  real_type AU_to_km = 1.49597870691E8;
  real_type Y_to_s   = 86400*365.25;
  real_type mu_SUN   = 1.32712440018E11;
  mu_SUN /= (AU_to_km*AU_to_km*AU_to_km)/(Y_to_s*Y_to_s);

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

  AstroLib::dvec3_t P0, V0, P1, V1;
  real_type         DV0, DV1;

  E.position( t0, P0.data() );
  E.velocity( t0, V0.data() );
  M.position( t0, P1.data() );
  M.velocity( t0, V1.data() );

  real_type DV  = AstroLib::minimum_DeltaV( mu_SUN, P0, V0, P1, V1, DV0, DV1, nullptr );
  fmt::print(
    "DV1+DV2     = {:.4}\n"
    "DV1         = {:.4}\n"
    "DV2         = {:.4}\n\n",
    DV, DV0, DV1
  );

  real_type DV2 = AstroLib::minimum_DeltaV2( mu_SUN, P0, V0, P1, V1, DV0, DV1 );
  fmt::print(
    "DV1^2+DV2^2 = {:.4}\n"
    "DV1         = {:.4}\n"
    "DV2         = {:.4}\n\n",
    DV2, DV0, DV1
  );

  std::cout << "All done folks!!\n";
  return 0;
}
