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

namespace Astro {

  /*
  //    ____                              _
  //   / ___|___  _ ____   _____ _ __ ___(_) ___  _ __
  //  | |   / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
  //  | |__| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
  //   \____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
  */

  extern real_type const one_UA_by_km;
  extern real_type const one_DAY_by_second;

  extern real_type const one_km_by_UA;
  extern real_type const one_m_by_UA;
  extern real_type const one_second_by_DAY;

  extern real_type const oneNewton_UA;

  extern real_type const gravity_kg_m_s2;
  extern real_type const gravity_kg_UA_DAY2;

  // sole
  extern real_type const muSun_km3s2;
  extern real_type const muSun_UA3DAY2;

  inline
  real_type
  UA_to_km( real_type pos )
  { return pos * one_UA_by_km; }

  inline
  real_type
  UA_by_DAY_to_km_by_s( real_type vel )
  { return vel * (one_UA_by_km/one_DAY_by_second); }

  inline
  real_type
  UA_by_DAY2_to_km_by_s2( real_type vel )
  { return vel * (one_UA_by_km/(one_DAY_by_second*one_DAY_by_second)); }

  inline
  real_type
  UA_by_DAY3_to_km_by_s3( real_type vel )
  { return vel * (one_UA_by_km/(one_DAY_by_second*one_DAY_by_second*one_DAY_by_second)); }

}

#endif
