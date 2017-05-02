/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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

///
/// file: alglin.hh
///

/*\
 |      _        _
 |     / \   ___| |_ _ __ ___
 |    / _ \ / __| __| '__/ _ \
 |   / ___ \\__ \ |_| | | (_) |
 |  /_/   \_\___/\__|_|  \___/
\*/

#ifndef ASTRO_HH
#define ASTRO_HH

#include <stdexcept>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#define DO_ERROR(MSG) { \
  std::ostringstream ost ; ost << MSG ; \
  throw std::runtime_error(ost.str()) ; }

#ifndef ASSERT
  #define ASSERT(COND,MSG) if ( !(COND) ) DO_ERROR(MSG)
#endif

#ifndef ASSERT_DEBUG
  #ifdef DEBUG
    #define ASSERT_DEBUG(COND,MSG) ASSERT(COND,MSG)
  #else
    #define ASSERT_DEBUG(COND,MSG)
  #endif
#endif

namespace astro {

  typedef double real_type ;
  typedef int    int_type ;

  using std::hypot ;
  using std::abs ;

  //============================================================================

  //! `m_pi` the value of \f$ \pi \f$.
  static real_type const m_pi = 3.141592653589793238462643383279502884197 ;

  //! `m_pi_2` the value of \f$ \pi/2 \f$.
  static real_type const m_pi_2 = 1.570796326794896619231321691639751442098 ;

  inline
  void
  zero3( real_type a[3] )
  { a[0] = a[1] = a[2] = 0 ;  }

  inline
  void
  copy3( real_type const a[3], real_type b[3] )
  { b[0] = a[0] ; b[1] = a[1] ; b[2] = a[2] ;  }

  inline
  real_type
  norm3( real_type const v[3] )
  { return hypot(hypot(v[0],v[1]),v[2]) ; }

  inline
  real_type
  dot3( real_type const a[3], real_type const b[3] )
  { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2] ; }

  inline
  real_type
  dist3( real_type const a[3], real_type const b[3] )
  { return hypot(hypot(a[0]-b[0],a[1]-b[1]),a[2]-b[2]) ; }

  inline
  void
  cross( real_type const a[3], real_type const b[3], real_type c[3] ) {
    c[0] = a[1]*b[2]-a[2]*b[1] ;
    c[1] = a[2]*b[0]-a[0]*b[2] ;
    c[2] = a[0]*b[1]-a[1]*b[0] ;
  }

  inline
  real_type
  pos( real_type const a )
  { return a > 0 ? a : 0 ; }

  inline
  real_type
  power2( real_type const a )
  { return a*a; }


  /*\
   % LAMBERT            Lambert-targeter for ballistic flights
   %                    (Izzo, and Lancaster, Blanchard & Gooding)
   %
   % This function solves any Lambert problem *robustly*. It uses two separate
   % solvers; the first one tried is a new and unpublished algorithm developed
   % by Dr. D. Izzo from the European Space Agency [1]. This version is extremely
   % fast, but especially for larger [m] it still fails quite frequently. In such
   % cases, a MUCH more robust algorithm is started (the one by Lancaster &
   % Blancard [2], with modifications, initial values and other improvements by
   % R.Gooding [3]), which is a lot slower partly because of its robustness.
   %
   % INPUT ARGUMENTS:
   % ======================================================================
   %    name        units    description
   % ======================================================================
   %   r1, r2        [m]     position vectors of the two terminal points.
   %     tf          [s]     time of flight to solve for
   %      m          [-]     specifies the number of complete orbits to complete
   %                         (should be an integer)
   %
   %  GM_central   [m^3/s^2]  std. grav. parameter (GM = mu) of the central body
   %
   %  units must be consistent, if r1 and r2 are in km then velocities must be
   %  in km/s and muC on [km^3/s^2]. The same happen is unit time is for example
   %  in days and not seconds the velocity must be in m/day and so on.
   %
   % OUTPUT ARGUMENTS:
   % ======================================================================
   %   name             units   description
   % ======================================================================
   %  V1, V2              [m/s]  terminal velocities at the end-points
   %  extremal_distances    [m]  minimum(1) and maximum(2) distance of the
   %                             spacecraft to the central body.
   %  exitflag             [-]   Integer containing information on why the
   %                             routine terminated. A value of +1 indicates
   %                             success; a normal exit. A value of -1
   %                             indicates that the given problem has no
   %                             solution and cannot be solved. A value of -2
   %                             indicates that both algorithms failed to find
   %                             a solution. This should never occur since
   %                             these problems are well-defined, and at the
   %                             very least it can be determined that the
   %                             problem has no solution. Nevertheless, it
   %                             still occurs sometimes for accidental
   %                             erroneous input, so it provides a basic
   %                             mechanism to check any application using this
   %                             algorithm.
   %
   %
   % References:
   %
   % [1] Izzo, D. ESA Advanced Concepts team. Code used available in MGA.M, on
   %     http://www.esa.int/gsp/ACT/inf/op/globopt.htm. Last retreived Nov, 2009.
   % [2] Lancaster, E.R. and Blanchard, R.C. "A unified form of Lambert's theorem."
   %     NASA technical note TN D-5368,1969.
   % [3] Gooding, R.H. "A procedure for the solution of Lambert's orbital boundary-value
   %     problem. Celestial Mechanics and Dynamical Astronomy, 48:145�165,1990.
   %
  \*/

  /*!
  //  Compute Lambert maneuvre
  //
  //  \param R1_m      initial position [m]
  //  \param R2_m      final position [m]
  //  \param tf_s      travel time [s]
  //  \param m         number of complete orbits [/]
  //  \param mu_m3_s2  gravitational constant [m^3/s^2]
  //  \param V1_ms     initial computed velocity [m/s]
  //  \param V2_ms     final computed velocity [m/s]
  //
  //  \return  +1 computation done correctly
  //           -1 No Solution found
  //           -2 Computation Failed
  */
  int_type
  lambert( real_type const R1_m[3],
           real_type const R2_m[3],
           real_type const tf_s,
           int       const m,
           real_type const mu_m3_s2,
           real_type       V1_ms[3],
           real_type       V2_ms[3] ) ;

}

#endif
