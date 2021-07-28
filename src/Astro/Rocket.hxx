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

#pragma once

#ifndef MECHATRONIX_ASTRO_LOCAL_ROCKET_dot_HH
#define MECHATRONIX_ASTRO_LOCAL_ROCKET_dot_HH

#ifndef ASTRO_dot_HH
  #include "Astro.hh"
#endif

namespace Astro {

  void
  odeXYZ(
    real_type       t,
    real_type       Thrust_r,
    real_type       Thrust_t,
    real_type       Thrust_n,
    real_type       m,
    real_type       mu,
    real_type const Z[],
    real_type       res[]
  );

  void
  odeEQ(
    real_type       t,
    real_type       Thrust_r,
    real_type       Thrust_t,
    real_type       Thrust_n,
    real_type       m,
    real_type       mu,
    real_type const Z[],
    real_type       res[]
  );
}

#endif
