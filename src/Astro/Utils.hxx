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

#ifndef ASTRO_LOCAL_ASTRO_UTILS_dot_HH
#define ASTRO_LOCAL_ASTRO_UTILS_dot_HH

#ifndef ASTRO_dot_HH
  #include "Astro.hh"
#endif

namespace Astro {

  inline
  real_type
  power2( real_type a )
  { return a*a; }

  inline
  real_type
  power3( real_type a )
  { return a*a*a; }

  //!
  //! Convert angle in degrees to radiants.
  //!
  static
  inline
  real_type
  degrees_to_radiants( real_type deg )
  { return 0.01745329251994329576923690768488612713443*deg; } // Pi/180*deg

  //!
  //! Convert angle in radiants to degrees.
  //!
  static
  inline
  real_type
  radiants_to_degrees( real_type rad )
  { return 57.29577951308232087679815481410517033240*rad; } // (180/Pi)*rad

  /*  _             __  _              __     _
  // |_)  /\  |\ | /__ |_ __ /\  |\ | /__ |  |_
  // | \ /--\ | \| \_| |_   /--\ | \| \_| |_ |_
  */
  //! Add or remove multiple of \f$ 2\pi \f$ to an angle in order to put it in the range \f$ [0,2\pi] \f$.
  inline
  void
  angle_in_range( real_type & ang ) {
    using Utils::m_2pi;
    ang = fmod( ang, m_2pi );
    while ( ang < 0     ) ang += m_2pi;
    while ( ang > m_2pi ) ang -= m_2pi;
  }

  //! Add or remove multiple of \f$ 2\pi \f$ to an angle  in order to put it in the range \f$ [-\pi,\pi]\f$.
  inline
  void
  angle_in_range_symm( real_type & ang ) {
    using Utils::m_2pi;
    using Utils::m_pi;
    ang = fmod( ang, m_2pi );
    while ( ang < -m_pi ) ang += m_2pi;
    while ( ang >  m_pi ) ang -= m_2pi;
  }

  inline
  void
  zero3( real_type a[3] )
  { a[0] = a[1] = a[2] = 0;  }

  inline
  void
  copy3( real_type const a[3], real_type b[3] )
  { b[0] = a[0]; b[1] = a[1]; b[2] = a[2]; }

  inline
  real_type
  norm3( real_type const v[3] )
  { return hypot(hypot(v[0],v[1]),v[2]); }

  inline
  real_type
  dot3( real_type const a[3], real_type const b[3] )
  { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

  inline
  real_type
  dist3( real_type const a[3], real_type const b[3] )
  { return hypot(hypot(a[0]-b[0],a[1]-b[1]),a[2]-b[2]); }

  inline
  void
  cross( real_type const a[3], real_type const b[3], real_type c[3] ) {
    c[0] = a[1]*b[2]-a[2]*b[1];
    c[1] = a[2]*b[0]-a[0]*b[2];
    c[2] = a[0]*b[1]-a[1]*b[0];
  }
}

#endif
