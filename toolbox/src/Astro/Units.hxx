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

// raccolta funzioni per il calcolo di traiettorie Kepleriane

#pragma once

#ifndef MECHATRONIX_ASTRO_LOCAL_UNITS_dot_HH
#define MECHATRONIX_ASTRO_LOCAL_UNITS_dot_HH

#ifndef ASTRO_dot_HH
  #include "Astro.hh"
#endif

namespace AstroLib {


  /*
  //    ____                              _
  //   / ___|___  _ ____   _____ _ __ ___(_) ___  _ __
  //  | |   / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
  //  | |__| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
  //   \____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
  */

  constexpr real_type const one_UA_to_km{1.49597870691E8}; // KM
  constexpr real_type const one_DAY_to_second{86400.0}; // sec

  constexpr real_type const one_km_to_UA{1/one_UA_to_km};
  constexpr real_type const one_m_to_UA{one_km_to_UA/1000};
  constexpr real_type const one_second_to_DAY{1/one_DAY_to_second};

  constexpr real_type const oneNewton_to_kg_UA_DAY2{one_m_to_UA/(one_second_to_DAY*one_second_to_DAY)}; // Kg * m/s^2 => Kg * UA / day^2

  constexpr real_type const gravity_kg_m_s2{9.80665};
  constexpr real_type const gravity_kg_UA_DAY2{gravity_kg_m_s2*oneNewton_to_kg_UA_DAY2};

  // costanti gravitazionali

  constexpr real_type const muSun_km3s2{1.32712440018E11}; // Km^3/s^2
  constexpr real_type const muSun_UA3DAY2{muSun_km3s2*(one_km_to_UA*one_km_to_UA*one_km_to_UA)/(one_second_to_DAY*one_second_to_DAY)}; // Km^3/s^2

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  real_type
  UA_to_km( real_type pos )
  { return pos * one_UA_to_km; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  real_type
  UA_by_DAY_to_km_by_s( real_type vel )
  { return vel * (one_UA_to_km/one_DAY_to_second); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  real_type
  UA_by_DAY2_to_km_by_s2( real_type vel )
  { return vel * (one_UA_to_km/(one_DAY_to_second*one_DAY_to_second)); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  real_type
  UA_by_DAY3_to_km_by_s3( real_type vel )
  { return vel * (one_UA_to_km/(one_DAY_to_second*one_DAY_to_second*one_DAY_to_second)); }

}

#endif
