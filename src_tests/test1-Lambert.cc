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
int
main() {

  AstroLib::dvec3_t P1, P2, V1, V2;

  real_type t1 = 95739; P1 << -0.622167037586749, -0.795634261754565, 0.000048989370791;
  real_type t2 = 96104; P2 << 3.110290382511421, 0.271081501337672, 0.279828671989069;
  real_type muSun_UA3DAY2 = 2.959122082856201e-04;
  int m = 0;
  int ok = AstroLib::Lambert( P1, P2, t2-t1, m, muSun_UA3DAY2, V1, V2 );
  std::cout << "ok = " << ok << "\n";


  AstroLib::Astro A, B;

  real_type AU_to_km = 1.49597870691E8;
  real_type Y_to_s   = 86400*365.25;
  real_type mu_SUN   = 1.32712440018E11;
  mu_SUN /= (AU_to_km*AU_to_km*AU_to_km)/(Y_to_s*Y_to_s);

  std::string name = "Earth";
  real_type t0   = 54000;
  real_type muS  = mu_SUN;
  real_type p    = 0.999723;
  real_type f    = 0;
  real_type g    = 0.016287;
  real_type h    = 0.3;
  real_type k    = 0;
  real_type L0   = 1.1;
  bool retrograde = false;

  A.setup_Equinoctial( name, t0, p, f, g, h, k, retrograde, L0, muS );

  f += 0.001;
  B.setup_Equinoctial( name, t0, p, f, g, h, k, retrograde, L0, muS );

  A.info( std::cout );
  B.info( std::cout );

  fmt::print(
    "A.L_orbital(...) = {}\n"
    "B.L_orbital(...) = {}\n",
    A.L_orbital(t0,0),
    B.L_orbital(t0,0)
  );

  real_type grad[6];
  real_type L = A.L_orbital_EQ_gradient( t0, grad );

  fmt::print(
    "g_p  = {}\n"
    "g_f  = {}\n"
    "g_g  = {}\n"
    "g_h  = {}\n"
    "g_k  = {}\n"
    "g_L0 = {}\n",
    grad[0], grad[1], grad[2],
    grad[3], grad[4], grad[5]
  );

  t0  = 0;
  muS = mu_SUN;
  p   = 2.97152;
  f   = 0.12134639999;
  g   = -0.9926102;
  h   = 0.0402935;
  k   = -0.149672;
  L0  = 11.4597;
  retrograde = false;

  A.setup_Equinoctial( name, t0, p, f, g, h, k, retrograde, L0, muS );
  A.info( std::cout );
  real_type x, y, z;
  A.position( t0, x, y, z );
  real_type res[4];
  real_type M = AstroLib::E_to_true_anomaly( Utils::m_2pi, hypot(f,g) );
  A.mean_anomaly_to_E( 2.5299780298288088e-06, res, 0 );

  fmt::print(
    "x = {}, y = {}, z = {}\n",
    x, y, z
  );

  std::cout << "All done folks!!\n";
  return 0;
}
