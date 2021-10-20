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

  real_type t1   = 95739;
  real_type P1[] = {-0.622167037586749, -0.795634261754565, 0.000048989370791};
  real_type t2   = 96104;
  real_type P2[] = { 3.110290382511421, 0.271081501337672, 0.279828671989069};
  real_type muSun_UA3DAY2 = 2.959122082856201e-04;
  real_type V1[3], V2[3];
  int m = 0;
  int ok = AstroLib::Lambert( P1, P2, t2-t1, m, muSun_UA3DAY2, V1, V2 );
  std::cout << "ok = " << ok << "\n";
  std::cout << "All done folks!!\n";
  return 0;
}
