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

namespace AstroLib {

  /*
  //    ____                              _
  //   / ___|___  _ ____   _____ _ __ ___(_) ___  _ __
  //  | |   / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
  //  | |__| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
  //   \____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
  */

  real_type const one_UA_to_km       = 1.49597870691E8; // KM
  real_type const one_DAY_to_second  = 86400.0; // sec

  real_type const one_km_to_UA       = 1/one_UA_to_km;
  real_type const one_m_to_UA        = one_km_to_UA/1000;
  real_type const one_second_to_DAY  = 1.0/one_DAY_to_second;

  real_type const oneNewton_to_kg_UA_DAY2 = one_m_to_UA/power2(one_second_to_DAY); // Kg * m/s^2 => Kg * UA / day^2

  real_type const gravity_kg_m_s2    = 9.80665;
  real_type const gravity_kg_UA_DAY2 = gravity_kg_m_s2 * oneNewton_to_kg_UA_DAY2;

  // costanti gravitazionali

  real_type const muSun_km3s2        = 1.32712440018E11; // Km^3/s^2
  real_type const muSun_UA3DAY2      = muSun_km3s2*power2(one_km_to_UA/one_second_to_DAY)*one_km_to_UA; // Km^3/s^2

}