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
 |      Via Mesiano 77, I-38050 Trento, Italy                               |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#if defined(__clang__)
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wunused-macros"
#elif defined(__llvm__) || defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wunused-macros"
#elif defined(_MSC_VER)
#pragma warning( disable : 4100 )
#pragma warning( disable : 4101 )
#endif

#include "Astro.hh"

namespace AstroLib {

  static bool const m_debug = true;

  real_type
  astro_x_position__xo( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type result__ = 1.0 / (t1 + t2 + 1) / (t4 * xo__f + t7 * xo__g + 1) * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type result__ = 1.0 / (t1 + t2 + 1) / (t4 * xo__f + t7 * xo__g + 1) * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_1( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_1_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_1_1( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_1_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type t15  = pow(t4 * xo__f + t7 * xo__g + 1, 2);
    real_type result__ = -t4 / (t1 + t2 + 1) / t15 * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_1_2( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_1_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type t15  = pow(t4 * xo__f + t7 * xo__g + 1, 2);
    real_type result__ = -t7 / (t1 + t2 + 1) / t15 * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_1_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_1_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t11  = 1.0 / (t1 * xo__f + t4 * xo__g + 1);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t26  = t15 * t15;
    real_type result__ = 1.0 / t15 * t11 * (2 * t4 * xo__k * xo__retrograde + 2 * t1 * xo__h) - 2 * xo__h / t26 * t11 * (t1 * (t13 - t14 + 1) + 2 * t4 * xo__k * xo__retrograde * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_1_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_1_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t3   = xo__h * xo__retrograde;
    real_type t4   = sin(xo__L);
    real_type t11  = 1.0 / (t1 * xo__f + t4 * xo__g + 1);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t25  = t15 * t15;
    real_type result__ = 1.0 / t15 * t11 * (-2 * t1 * xo__k + 2 * t4 * t3) - 2 * xo__k / t25 * t11 * (t1 * (t13 - t14 + 1) + 2 * t4 * xo__k * t3);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_1_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_1_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t3   = t1 - t2 + 1;
    real_type t4   = sin(xo__L);
    real_type t6   = xo__h * xo__retrograde;
    real_type t7   = cos(xo__L);
    real_type t14  = t4 * xo__g + t7 * xo__f + 1;
    real_type t18  = 1.0 / (t1 + t2 + 1);
    real_type t25  = t14 * t14;
    real_type result__ = t18 / t14 * (2 * t6 * t7 * xo__k - t3 * t4) - (-t4 * xo__f + t7 * xo__g) * t18 / t25 * (2 * t4 * xo__k * t6 + t3 * t7);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_1_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_1_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type result__ = 2 / (t9 + t10 + 1) / (t2 * xo__g + t3 * xo__f + 1) * t2 * xo__h * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_1_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type t16  = pow(t4 * xo__f + t7 * xo__g + 1, 2);
    real_type result__ = -t4 / (t1 + t2 + 1) / t16 * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_2( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_2_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type t15  = t4 * xo__f + t7 * xo__g + 1;
    real_type t16  = t15 * t15;
    real_type t22  = t4 * t4;
    real_type result__ = 2 * t22 / (t1 + t2 + 1) / t16 / t15 * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_2_2( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_2_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type t15  = t4 * xo__f + t7 * xo__g + 1;
    real_type t16  = t15 * t15;
    real_type result__ = 2 * t7 * t4 / (t1 + t2 + 1) / t16 / t15 * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_2_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_2_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t12  = pow(t1 * xo__f + t4 * xo__g + 1, 2);
    real_type t13  = 1.0 / t12;
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t16  = t14 + t15 + 1;
    real_type t30  = t16 * t16;
    real_type result__ = -t1 / t16 * t13 * (2 * t4 * xo__k * xo__retrograde + 2 * t1 * xo__h) * xo__p + 2 * xo__h * t1 / t30 * t13 * (t1 * (t14 - t15 + 1) + 2 * t4 * xo__k * xo__retrograde * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_2_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_2_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t3   = xo__h * xo__retrograde;
    real_type t4   = sin(xo__L);
    real_type t12  = pow(t1 * xo__f + t4 * xo__g + 1, 2);
    real_type t13  = 1.0 / t12;
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t16  = t14 + t15 + 1;
    real_type t29  = t16 * t16;
    real_type result__ = -t1 / t16 * t13 * (-2 * t1 * xo__k + 2 * t3 * t4) * xo__p + 2 * xo__k * t1 / t29 * t13 * (t1 * (t14 - t15 + 1) + 2 * t4 * xo__k * t3) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_2_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_2_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t3   = t1 - t2 + 1;
    real_type t4   = sin(xo__L);
    real_type t6   = xo__h * xo__retrograde;
    real_type t7   = cos(xo__L);
    real_type t15  = t4 * xo__g + t7 * xo__f + 1;
    real_type t16  = t15 * t15;
    real_type t19  = 1.0 / (t1 + t2 + 1);
    real_type t20  = t19 / t16;
    real_type t28  = (2 * t4 * t6 * xo__k + t3 * t7) * xo__p;
    real_type result__ = -t7 * t20 * (2 * t6 * t7 * xo__k - t3 * t4) * xo__p + 2 * (-t4 * xo__f + t7 * xo__g) * t7 * t19 / t16 / t15 * t28 + t4 * t20 * t28;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_2_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_2_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t3   = sin(xo__L);
    real_type t4   = cos(xo__L);
    real_type t8   = pow(t3 * xo__g + t4 * xo__f + 1, 2);
    real_type t11  = xo__h * xo__h;
    real_type t12  = xo__k * xo__k;
    real_type result__ = -2 * t4 / (t11 + t12 + 1) / t8 * t3 * xo__p * xo__h * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_2_8( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type t16  = pow(t4 * xo__f + t7 * xo__g + 1, 2);
    real_type result__ = -t7 / (t1 + t2 + 1) / t16 * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_3_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type t15  = t4 * xo__f + t7 * xo__g + 1;
    real_type t16  = t15 * t15;
    real_type t22  = t7 * t7;
    real_type result__ = 2 * t22 / (t1 + t2 + 1) / t16 / t15 * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_3_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_3_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t12  = pow(t1 * xo__f + t4 * xo__g + 1, 2);
    real_type t13  = 1.0 / t12;
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t16  = t14 + t15 + 1;
    real_type t30  = t16 * t16;
    real_type result__ = -t4 / t16 * t13 * (2 * t4 * xo__k * xo__retrograde + 2 * t1 * xo__h) * xo__p + 2 * xo__h * t4 / t30 * t13 * (t1 * (t14 - t15 + 1) + 2 * t4 * xo__k * xo__retrograde * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_3_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_3_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t3   = xo__h * xo__retrograde;
    real_type t4   = sin(xo__L);
    real_type t12  = pow(t1 * xo__f + t4 * xo__g + 1, 2);
    real_type t13  = 1.0 / t12;
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t16  = t14 + t15 + 1;
    real_type t29  = t16 * t16;
    real_type result__ = -t4 / t16 * t13 * (-2 * t1 * xo__k + 2 * t3 * t4) * xo__p + 2 * xo__k * t4 / t29 * t13 * (t1 * (t14 - t15 + 1) + 2 * t4 * xo__k * t3) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_3_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_3_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t3   = t1 - t2 + 1;
    real_type t4   = sin(xo__L);
    real_type t6   = xo__h * xo__retrograde;
    real_type t7   = cos(xo__L);
    real_type t15  = t4 * xo__g + t7 * xo__f + 1;
    real_type t16  = t15 * t15;
    real_type t19  = 1.0 / (t1 + t2 + 1);
    real_type t20  = t19 / t16;
    real_type t28  = (2 * t4 * t6 * xo__k + t3 * t7) * xo__p;
    real_type result__ = -t4 * t20 * (2 * t6 * t7 * xo__k - t3 * t4) * xo__p + 2 * (-t4 * xo__f + t7 * xo__g) * t4 * t19 / t16 / t15 * t28 - t7 * t20 * t28;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_3_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_3_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t3   = sin(xo__L);
    real_type t4   = t3 * t3;
    real_type t5   = cos(xo__L);
    real_type t9   = pow(t3 * xo__g + t5 * xo__f + 1, 2);
    real_type t12  = xo__h * xo__h;
    real_type t13  = xo__k * xo__k;
    real_type result__ = -2 / (t12 + t13 + 1) / t9 * t4 * xo__p * xo__h * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_3_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t12  = 1.0 / (t1 * xo__f + t4 * xo__g + 1);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t27  = t15 * t15;
    real_type result__ = 1.0 / t15 * t12 * (2 * t4 * xo__k * xo__retrograde + 2 * t1 * xo__h) * xo__p - 2 * xo__h / t27 * t12 * (t1 * (t13 - t14 + 1) + 2 * t4 * xo__k * xo__retrograde * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_4_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t7   = 1.0 / (t1 * xo__f + t4 * xo__g + 1);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t21  = t10 * t10;
    real_type t23  = 1.0 / t21 * t7;
    real_type t34  = (t1 * (t8 - t9 + 1) + 2 * t4 * xo__k * xo__retrograde * xo__h) * xo__p;
    real_type result__ = 2 / t10 * t7 * t1 * xo__p - 4 * xo__h * t23 * (2 * t4 * xo__k * xo__retrograde + 2 * t1 * xo__h) * xo__p + 8 * t8 / t21 / t10 * t7 * t34 - 2 * t23 * t34;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_4_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_4_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t7   = 1.0 / (t2 * xo__g + t3 * xo__f + 1);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t22  = t11 * t11;
    real_type t24  = 1.0 / t22 * t7;
    real_type t29  = xo__h * xo__retrograde;
    real_type result__ = 2 / t11 * t7 * t2 * xo__p * xo__retrograde - 2 * xo__k * t24 * (2 * t2 * xo__k * xo__retrograde + 2 * t3 * xo__h) * xo__p - 2 * xo__h * t24 * (2 * t2 * t29 - 2 * t3 * xo__k) * xo__p + 8 * xo__k * xo__h / t22 / t11 * t7 * (t3 * (t9 - t10 + 1) + 2 * t2 * xo__k * t29) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_4_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_4_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = xo__k * xo__retrograde;
    real_type t4   = cos(xo__L);
    real_type t11  = t1 * xo__g + t4 * xo__f + 1;
    real_type t12  = 1.0 / t11;
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t16  = 1.0 / t15;
    real_type t24  = t11 * t11;
    real_type t25  = 1.0 / t24;
    real_type t29  = -t1 * xo__f + t4 * xo__g;
    real_type t32  = t13 - t14 + 1;
    real_type t34  = xo__h * xo__retrograde;
    real_type t40  = t15 * t15;
    real_type t41  = 1.0 / t40;
    real_type result__ = t16 * t12 * (-2 * t1 * xo__h + 2 * t4 * t3) * xo__p - t29 * t16 * t25 * (2 * t1 * t3 + 2 * t4 * xo__h) * xo__p - 2 * xo__h * t41 * t12 * (2 * t4 * xo__k * t34 - t1 * t32) * xo__p + 2 * t29 * xo__h * t41 * t25 * (2 * t1 * xo__k * t34 + t4 * t32) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_4_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_4_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t8   = 1.0 / (t2 * xo__g + t3 * xo__f + 1) * t2;
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t18  = t11 * t11;
    real_type result__ = 2 / t11 * t8 * xo__p * xo__k - 4 / t18 * t8 * xo__k * t9 * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_4_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t3   = xo__h * xo__retrograde;
    real_type t4   = sin(xo__L);
    real_type t12  = 1.0 / (t1 * xo__f + t4 * xo__g + 1);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t26  = t15 * t15;
    real_type result__ = 1.0 / t15 * t12 * (-2 * t1 * xo__k + 2 * t4 * t3) * xo__p - 2 * xo__k / t26 * t12 * (t1 * (t13 - t14 + 1) + 2 * t4 * xo__k * t3) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_5_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t7   = 1.0 / (t1 * xo__f + t4 * xo__g + 1);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t16  = xo__h * xo__retrograde;
    real_type t21  = t10 * t10;
    real_type t23  = 1.0 / t21 * t7;
    real_type t33  = (t1 * (t8 - t9 + 1) + 2 * t4 * xo__k * t16) * xo__p;
    real_type result__ = -2 / t10 * t7 * t1 * xo__p - 4 * xo__k * t23 * (-2 * t1 * xo__k + 2 * t4 * t16) * xo__p + 8 * t9 / t21 / t10 * t7 * t33 - 2 * t23 * t33;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_5_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_5_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t2   = t1 * xo__k;
    real_type t3   = xo__h * xo__retrograde;
    real_type t4   = cos(xo__L);
    real_type t11  = t1 * xo__g + t4 * xo__f + 1;
    real_type t12  = 1.0 / t11;
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t16  = 1.0 / t15;
    real_type t19  = t4 * xo__k;
    real_type t24  = t11 * t11;
    real_type t25  = 1.0 / t24;
    real_type t29  = -t1 * xo__f + t4 * xo__g;
    real_type t32  = t13 - t14 + 1;
    real_type t38  = t15 * t15;
    real_type t39  = 1.0 / t38;
    real_type result__ = t16 * t12 * (2 * t4 * t3 + 2 * t2) * xo__p - t29 * t16 * t25 * (2 * t1 * t3 - 2 * t19) * xo__p - 2 * xo__k * t39 * t12 * (-t1 * t32 + 2 * t19 * t3) * xo__p + 2 * t29 * xo__k * t39 * t25 * (2 * t2 * t3 + t4 * t32) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_5_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_5_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__p * xo__h;
    real_type t2   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t8   = 1.0 / (t2 * xo__g + t3 * xo__f + 1) * t2;
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t17  = t11 * t11;
    real_type result__ = 2 / t11 * t8 * t1 - 4 / t17 * t8 * t10 * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_5_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t3   = t1 - t2 + 1;
    real_type t4   = sin(xo__L);
    real_type t6   = xo__h * xo__retrograde;
    real_type t7   = cos(xo__L);
    real_type t15  = t4 * xo__g + t7 * xo__f + 1;
    real_type t18  = 1.0 / (t1 + t2 + 1);
    real_type t27  = t15 * t15;
    real_type result__ = t18 / t15 * (2 * t7 * xo__k * t6 - t4 * t3) * xo__p - (-t4 * xo__f + t7 * xo__g) * t18 / t27 * (2 * t4 * xo__k * t6 + t7 * t3) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_6_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t3   = t1 - t2 + 1;
    real_type t4   = cos(xo__L);
    real_type t6   = xo__h * xo__retrograde;
    real_type t7   = sin(xo__L);
    real_type t11  = -2 * t7 * xo__k * t6 - t4 * t3;
    real_type t13  = t4 * xo__f;
    real_type t14  = t7 * xo__g;
    real_type t15  = t13 + t14 + 1;
    real_type t18  = 1.0 / (t1 + t2 + 1);
    real_type t27  = t15 * t15;
    real_type t29  = t18 / t27;
    real_type t32  = t4 * xo__g - t7 * xo__f;
    real_type t37  = -t11 * xo__p;
    real_type t41  = t32 * t32;
    real_type result__ = t18 / t15 * t11 * xo__p - 2 * t32 * t29 * (2 * t4 * xo__k * t6 - t7 * t3) * xo__p + 2 * t41 * t18 / t27 / t15 * t37 - (-t13 - t14) * t29 * t37;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_6_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_6_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = xo__p * xo__h * xo__k;
    real_type t3   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t7   = t3 * xo__f + t5 * xo__g + 1;
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type t13  = 1.0 / (t10 + t11 + 1);
    real_type t16  = t7 * t7;
    real_type result__ = 2 * t13 / t7 * t3 * t2 - 2 * (t3 * xo__g - t5 * xo__f) * t13 / t16 * t5 * t2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_6_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t3   = sin(xo__L);
    real_type t4   = cos(xo__L);
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type result__ = 2 / (t10 + t11 + 1) / (t3 * xo__g + t4 * xo__f + 1) * t3 * xo__p * xo__h * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_position__xo_D_7_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_position__xo_D_7_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type result__ = -1.0 / (t1 + t2 + 1) / (t5 * xo__g + t8 * xo__f + 1) * (t5 * (t1 - t2 - 1) * xo__retrograde - 2 * t8 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type result__ = -1.0 / (t1 + t2 + 1) / (t5 * xo__g + t8 * xo__f + 1) * (t5 * (t1 - t2 - 1) * xo__retrograde - 2 * t8 * xo__k * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_1( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_1_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_1_1( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_1_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type t15  = pow(t5 * xo__g + t8 * xo__f + 1, 2);
    real_type result__ = t8 / (t1 + t2 + 1) / t15 * (t5 * (t1 - t2 - 1) * xo__retrograde - 2 * t8 * xo__k * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_1_2( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_1_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type t15  = pow(t5 * xo__g + t8 * xo__f + 1, 2);
    real_type result__ = t5 / (t1 + t2 + 1) / t15 * (t5 * (t1 - t2 - 1) * xo__retrograde - 2 * t8 * xo__k * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_1_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_1_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t11  = 1.0 / (t1 * xo__f + t4 * xo__g + 1);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t26  = t15 * t15;
    real_type result__ = -1.0 / t15 * t11 * (2 * t4 * xo__h * xo__retrograde - 2 * t1 * xo__k) + 2 * xo__h / t26 * t11 * (t4 * (t13 - t14 - 1) * xo__retrograde - 2 * t1 * xo__k * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_1_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_1_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t4   = cos(xo__L);
    real_type t11  = 1.0 / (t2 * xo__g + t4 * xo__f + 1);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t26  = t15 * t15;
    real_type result__ = -1.0 / t15 * t11 * (-2 * t2 * xo__k * xo__retrograde - 2 * t4 * xo__h) + 2 * xo__k / t26 * t11 * (t2 * (t13 - t14 - 1) * xo__retrograde - 2 * t4 * xo__k * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_1_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_1_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = (t1 - t2 - 1) * xo__retrograde;
    real_type t5   = cos(xo__L);
    real_type t7   = xo__k * xo__h;
    real_type t8   = sin(xo__L);
    real_type t14  = t5 * xo__f + t8 * xo__g + 1;
    real_type t18  = 1.0 / (t1 + t2 + 1);
    real_type t24  = t14 * t14;
    real_type result__ = -t18 / t14 * (t5 * t4 + 2 * t8 * t7) + (t5 * xo__g - t8 * xo__f) * t18 / t24 * (t8 * t4 - 2 * t5 * t7);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_1_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_1_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = sin(xo__L);
    real_type t6   = cos(xo__L);
    real_type result__ = -1.0 / (t1 + t2 + 1) / (t4 * xo__g + t6 * xo__f + 1) * t4 * (t1 - t2 - 1);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_1_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type t16  = pow(t5 * xo__g + t8 * xo__f + 1, 2);
    real_type result__ = t8 / (t1 + t2 + 1) / t16 * (t5 * (t1 - t2 - 1) * xo__retrograde - 2 * t8 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_2( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_2_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type t15  = t5 * xo__g + t8 * xo__f + 1;
    real_type t16  = t15 * t15;
    real_type t22  = t8 * t8;
    real_type result__ = -2 * t22 / (t1 + t2 + 1) / t16 / t15 * (t5 * (t1 - t2 - 1) * xo__retrograde - 2 * t8 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_2_2( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_2_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type t15  = t5 * xo__g + t8 * xo__f + 1;
    real_type t16  = t15 * t15;
    real_type result__ = -2 * t5 * t8 / (t1 + t2 + 1) / t16 / t15 * (t5 * (t1 - t2 - 1) * xo__retrograde - 2 * t8 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_2_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_2_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t12  = pow(t1 * xo__f + t4 * xo__g + 1, 2);
    real_type t13  = 1.0 / t12;
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t16  = t14 + t15 + 1;
    real_type t30  = t16 * t16;
    real_type result__ = t1 / t16 * t13 * (2 * t4 * xo__h * xo__retrograde - 2 * t1 * xo__k) * xo__p - 2 * xo__h * t1 / t30 * t13 * (t4 * (t14 - t15 - 1) * xo__retrograde - 2 * t1 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_2_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_2_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t4   = cos(xo__L);
    real_type t12  = pow(t2 * xo__g + t4 * xo__f + 1, 2);
    real_type t13  = 1.0 / t12;
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t16  = t14 + t15 + 1;
    real_type t30  = t16 * t16;
    real_type result__ = t4 / t16 * t13 * (-2 * t2 * xo__k * xo__retrograde - 2 * t4 * xo__h) * xo__p - 2 * xo__k * t4 / t30 * t13 * (t2 * (t14 - t15 - 1) * xo__retrograde - 2 * t4 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_2_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_2_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = (t1 - t2 - 1) * xo__retrograde;
    real_type t5   = cos(xo__L);
    real_type t7   = xo__k * xo__h;
    real_type t8   = sin(xo__L);
    real_type t15  = t5 * xo__f + t8 * xo__g + 1;
    real_type t16  = t15 * t15;
    real_type t19  = 1.0 / (t1 + t2 + 1);
    real_type t20  = t19 / t16;
    real_type t27  = (t4 * t8 - 2 * t5 * t7) * xo__p;
    real_type result__ = t5 * t20 * (t4 * t5 + 2 * t7 * t8) * xo__p - 2 * (t5 * xo__g - t8 * xo__f) * t5 * t19 / t16 / t15 * t27 - t8 * t20 * t27;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_2_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_2_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t7   = cos(xo__L);
    real_type t11  = pow(t5 * xo__g + t7 * xo__f + 1, 2);
    real_type result__ = t7 / (t1 + t2 + 1) / t11 * t5 * (t1 - t2 - 1) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_2_8( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type t16  = pow(t5 * xo__g + t8 * xo__f + 1, 2);
    real_type result__ = t5 / (t1 + t2 + 1) / t16 * (t5 * (t1 - t2 - 1) * xo__retrograde - 2 * t8 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_3_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type t15  = t5 * xo__g + t8 * xo__f + 1;
    real_type t16  = t15 * t15;
    real_type t22  = t5 * t5;
    real_type result__ = -2 * t22 / (t1 + t2 + 1) / t16 / t15 * (t5 * (t1 - t2 - 1) * xo__retrograde - 2 * t8 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_3_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_3_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t12  = pow(t1 * xo__f + t4 * xo__g + 1, 2);
    real_type t13  = 1.0 / t12;
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t16  = t14 + t15 + 1;
    real_type t30  = t16 * t16;
    real_type result__ = t4 / t16 * t13 * (2 * t4 * xo__h * xo__retrograde - 2 * t1 * xo__k) * xo__p - 2 * xo__h * t4 / t30 * t13 * (t4 * (t14 - t15 - 1) * xo__retrograde - 2 * t1 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_3_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_3_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t4   = cos(xo__L);
    real_type t12  = pow(t2 * xo__g + t4 * xo__f + 1, 2);
    real_type t13  = 1.0 / t12;
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t16  = t14 + t15 + 1;
    real_type t30  = t16 * t16;
    real_type result__ = t2 / t16 * t13 * (-2 * t2 * xo__k * xo__retrograde - 2 * t4 * xo__h) * xo__p - 2 * xo__k * t2 / t30 * t13 * (t2 * (t14 - t15 - 1) * xo__retrograde - 2 * t4 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_3_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_3_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = (t1 - t2 - 1) * xo__retrograde;
    real_type t5   = cos(xo__L);
    real_type t7   = xo__k * xo__h;
    real_type t8   = sin(xo__L);
    real_type t15  = t5 * xo__f + t8 * xo__g + 1;
    real_type t16  = t15 * t15;
    real_type t19  = 1.0 / (t1 + t2 + 1);
    real_type t20  = t19 / t16;
    real_type t27  = (t4 * t8 - 2 * t5 * t7) * xo__p;
    real_type result__ = t8 * t20 * (t4 * t5 + 2 * t7 * t8) * xo__p - 2 * (t5 * xo__g - t8 * xo__f) * t8 * t19 / t16 / t15 * t27 + t5 * t20 * t27;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_3_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_3_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t6   = t5 * t5;
    real_type t7   = cos(xo__L);
    real_type t11  = pow(t5 * xo__g + t7 * xo__f + 1, 2);
    real_type result__ = 1.0 / (t1 + t2 + 1) / t11 * t6 * (t1 - t2 - 1) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_3_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t12  = 1.0 / (t1 * xo__f + t4 * xo__g + 1);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t27  = t15 * t15;
    real_type result__ = -1.0 / t15 * t12 * (2 * t4 * xo__h * xo__retrograde - 2 * t1 * xo__k) * xo__p + 2 * xo__h / t27 * t12 * (t4 * (t13 - t14 - 1) * xo__retrograde - 2 * t1 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_4_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t7   = 1.0 / (t2 * xo__g + t3 * xo__f + 1);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t22  = t11 * t11;
    real_type t24  = 1.0 / t22 * t7;
    real_type t35  = (t2 * (t9 - t10 - 1) * xo__retrograde - 2 * t3 * xo__k * xo__h) * xo__p;
    real_type result__ = -2 / t11 * t7 * t2 * xo__p * xo__retrograde + 4 * xo__h * t24 * (2 * t2 * xo__h * xo__retrograde - 2 * t3 * xo__k) * xo__p - 8 * t9 / t22 / t11 * t7 * t35 + 2 * t24 * t35;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_4_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_4_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t7   = 1.0 / (t1 * xo__f + t4 * xo__g + 1);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t21  = t10 * t10;
    real_type t23  = 1.0 / t21 * t7;
    real_type result__ = 2 / t10 * t7 * t1 * xo__p + 2 * xo__k * t23 * (2 * t4 * xo__h * xo__retrograde - 2 * t1 * xo__k) * xo__p + 2 * xo__h * t23 * (-2 * t4 * xo__k * xo__retrograde - 2 * t1 * xo__h) * xo__p - 8 * xo__k * xo__h / t21 / t10 * t7 * (t4 * (t8 - t9 - 1) * xo__retrograde - 2 * t1 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_4_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_4_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = xo__h * xo__retrograde;
    real_type t4   = cos(xo__L);
    real_type t11  = t1 * xo__g + t4 * xo__f + 1;
    real_type t12  = 1.0 / t11;
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t16  = 1.0 / t15;
    real_type t24  = t11 * t11;
    real_type t25  = 1.0 / t24;
    real_type t29  = -t1 * xo__f + t4 * xo__g;
    real_type t33  = (t13 - t14 - 1) * xo__retrograde;
    real_type t35  = xo__k * xo__h;
    real_type t40  = t15 * t15;
    real_type t41  = 1.0 / t40;
    real_type result__ = -t16 * t12 * (2 * t1 * xo__k + 2 * t4 * t3) * xo__p + t29 * t16 * t25 * (2 * t1 * t3 - 2 * t4 * xo__k) * xo__p + 2 * xo__h * t41 * t12 * (2 * t1 * t35 + t4 * t33) * xo__p - 2 * t29 * xo__h * t41 * t25 * (t1 * t33 - 2 * t4 * t35) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_4_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_4_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t7   = 1.0 / (t2 * xo__g + t3 * xo__f + 1);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t18  = t11 * t11;
    real_type result__ = -2 / t11 * t7 * t2 * xo__p * xo__h + 2 * xo__h / t18 * t7 * t2 * (t9 - t10 - 1) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_4_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t4   = cos(xo__L);
    real_type t12  = 1.0 / (t2 * xo__g + t4 * xo__f + 1);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t27  = t15 * t15;
    real_type result__ = -1.0 / t15 * t12 * (-2 * t2 * xo__k * xo__retrograde - 2 * t4 * xo__h) * xo__p + 2 * xo__k / t27 * t12 * (t2 * (t13 - t14 - 1) * xo__retrograde - 2 * t4 * xo__k * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_5_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t7   = 1.0 / (t2 * xo__g + t3 * xo__f + 1);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t22  = t11 * t11;
    real_type t24  = 1.0 / t22 * t7;
    real_type t35  = (t2 * (t9 - t10 - 1) * xo__retrograde - 2 * t3 * xo__k * xo__h) * xo__p;
    real_type result__ = 2 / t11 * t7 * t2 * xo__p * xo__retrograde + 4 * xo__k * t24 * (-2 * t2 * xo__k * xo__retrograde - 2 * t3 * xo__h) * xo__p - 8 * t10 / t22 / t11 * t7 * t35 + 2 * t24 * t35;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_5_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_5_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__k * xo__retrograde;
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t11  = t2 * xo__f + t4 * xo__g + 1;
    real_type t12  = 1.0 / t11;
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t15  = t13 + t14 + 1;
    real_type t16  = 1.0 / t15;
    real_type t24  = t11 * t11;
    real_type t25  = 1.0 / t24;
    real_type t29  = t2 * xo__g - t4 * xo__f;
    real_type t33  = (t13 - t14 - 1) * xo__retrograde;
    real_type t35  = xo__k * xo__h;
    real_type t40  = t15 * t15;
    real_type t41  = 1.0 / t40;
    real_type result__ = -t16 * t12 * (-2 * t2 * t1 + 2 * t4 * xo__h) * xo__p + t29 * t16 * t25 * (-2 * t4 * t1 - 2 * t2 * xo__h) * xo__p + 2 * xo__k * t41 * t12 * (t2 * t33 + 2 * t4 * t35) * xo__p - 2 * t29 * xo__k * t41 * t25 * (-2 * t2 * t35 + t4 * t33) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_5_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_5_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t7   = 1.0 / (t2 * xo__g + t3 * xo__f + 1);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t18  = t11 * t11;
    real_type result__ = 2 / t11 * t7 * t2 * xo__p * xo__k + 2 * xo__k / t18 * t7 * t2 * (t9 - t10 - 1) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_5_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = (t1 - t2 - 1) * xo__retrograde;
    real_type t5   = cos(xo__L);
    real_type t7   = xo__k * xo__h;
    real_type t8   = sin(xo__L);
    real_type t15  = t5 * xo__f + t8 * xo__g + 1;
    real_type t18  = 1.0 / (t1 + t2 + 1);
    real_type t26  = t15 * t15;
    real_type result__ = -t18 / t15 * (t4 * t5 + 2 * t7 * t8) * xo__p + (t5 * xo__g - t8 * xo__f) * t18 / t26 * (t4 * t8 - 2 * t5 * t7) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_6_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = (t1 - t2 - 1) * xo__retrograde;
    real_type t5   = sin(xo__L);
    real_type t7   = xo__k * xo__h;
    real_type t8   = cos(xo__L);
    real_type t11  = -t4 * t5 + 2 * t7 * t8;
    real_type t13  = t8 * xo__f;
    real_type t14  = t5 * xo__g;
    real_type t15  = 1 + t13 + t14;
    real_type t18  = 1.0 / (t1 + t2 + 1);
    real_type t26  = t15 * t15;
    real_type t28  = t18 / t26;
    real_type t31  = -t5 * xo__f + t8 * xo__g;
    real_type t36  = -t11 * xo__p;
    real_type t40  = t31 * t31;
    real_type result__ = -t18 / t15 * t11 * xo__p + 2 * t31 * t28 * (t4 * t8 + 2 * t5 * t7) * xo__p - 2 * t40 * t18 / t26 / t15 * t36 + (-t13 - t14) * t28 * t36;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_6_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_6_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = (t1 - t2 - 1) * xo__p;
    real_type t5   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type t9   = t5 * xo__f + t7 * xo__g + 1;
    real_type t13  = 1.0 / (t1 + t2 + 1);
    real_type t17  = t9 * t9;
    real_type result__ = -t13 / t9 * t5 * t4 + (t5 * xo__g - t7 * xo__f) * t13 / t17 * t7 * t4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_6_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t5   = sin(xo__L);
    real_type t6   = cos(xo__L);
    real_type result__ = -1.0 / (t1 + t2 + 1) / (t5 * xo__g + t6 * xo__f + 1) * t5 * (t1 - t2 - 1) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_position__xo_D_7_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_position__xo_D_7_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t12  = xo__h * xo__h;
    real_type t13  = xo__k * xo__k;
    real_type result__ = 2 / (t12 + t13 + 1) / (t2 * xo__f + t4 * xo__g + 1) * (-t2 * xo__k * xo__retrograde + t4 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t12  = xo__h * xo__h;
    real_type t13  = xo__k * xo__k;
    real_type result__ = 2 / (t12 + t13 + 1) / (t2 * xo__f + t4 * xo__g + 1) * (-t2 * xo__k * xo__retrograde + t4 * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_1( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_1_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_1_1( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_1_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t10  = pow(t2 * xo__f + t4 * xo__g + 1, 2);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type result__ = -2 * t2 / (t13 + t14 + 1) / t10 * (-t2 * xo__k * xo__retrograde + t4 * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_1_2( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_1_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t10  = pow(t2 * xo__f + t4 * xo__g + 1, 2);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type result__ = -2 * t4 / (t13 + t14 + 1) / t10 * (-t2 * xo__k * xo__retrograde + t4 * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_1_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_1_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t2   = cos(xo__L);
    real_type t6   = 1.0 / (t1 * xo__g + t2 * xo__f + 1);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t19  = t10 * t10;
    real_type result__ = 2 / t10 * t6 * t1 - 4 * xo__h / t19 * t6 * (-t2 * xo__k * xo__retrograde + t1 * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_1_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_1_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t7   = 1.0 / (t1 * xo__f + t4 * xo__g + 1);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t20  = t10 * t10;
    real_type result__ = -2 / t10 * t7 * t1 * xo__retrograde - 4 * xo__k / t20 * t7 * (-t1 * xo__k * xo__retrograde + t4 * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_1_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_1_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__k * xo__retrograde;
    real_type t2   = sin(xo__L);
    real_type t4   = cos(xo__L);
    real_type t9   = t2 * xo__g + t4 * xo__f + 1;
    real_type t12  = xo__h * xo__h;
    real_type t13  = xo__k * xo__k;
    real_type t15  = 1.0 / (t12 + t13 + 1);
    real_type t20  = t9 * t9;
    real_type result__ = 2 * t15 / t9 * (t2 * t1 + t4 * xo__h) - 2 * (-t2 * xo__f + t4 * xo__g) * t15 / t20 * (-t4 * t1 + t2 * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_1_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_1_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type result__ = -2 / (t8 + t9 + 1) / (t1 * xo__f + t4 * xo__g + 1) * t1 * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_1_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t11  = pow(t2 * xo__f + t4 * xo__g + 1, 2);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type result__ = -2 * t2 / (t13 + t14 + 1) / t11 * (-t2 * xo__k * xo__retrograde + t4 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_2( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_2_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t10  = t2 * xo__f + t4 * xo__g + 1;
    real_type t11  = t10 * t10;
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t19  = t2 * t2;
    real_type result__ = 4 * t19 / (t14 + t15 + 1) / t11 / t10 * (-t2 * xo__k * xo__retrograde + t4 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_2_2( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_2_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t10  = t2 * xo__f + t4 * xo__g + 1;
    real_type t11  = t10 * t10;
    real_type t15  = xo__h * xo__h;
    real_type t16  = xo__k * xo__k;
    real_type result__ = 4 * t4 * t2 / (t15 + t16 + 1) / t11 / t10 * (-t2 * xo__k * xo__retrograde + t4 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_2_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_2_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t7   = pow(t1 * xo__g + t3 * xo__f + 1, 2);
    real_type t8   = 1.0 / t7;
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t23  = t11 * t11;
    real_type result__ = -2 * t3 / t11 * t8 * t1 * xo__p + 4 * xo__h * t3 / t23 * t8 * (-t3 * xo__k * xo__retrograde + t1 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_2_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_2_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t3   = t2 * t2;
    real_type t5   = sin(xo__L);
    real_type t8   = pow(t2 * xo__f + t5 * xo__g + 1, 2);
    real_type t9   = 1.0 / t8;
    real_type t11  = xo__h * xo__h;
    real_type t12  = xo__k * xo__k;
    real_type t13  = t11 + t12 + 1;
    real_type t24  = t13 * t13;
    real_type result__ = 2 / t13 * t9 * t3 * xo__p * xo__retrograde + 4 * xo__k * t2 / t24 * t9 * (-t2 * xo__k * xo__retrograde + t5 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_2_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_2_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__k * xo__retrograde;
    real_type t2   = sin(xo__L);
    real_type t4   = cos(xo__L);
    real_type t10  = t2 * xo__g + t4 * xo__f + 1;
    real_type t11  = t10 * t10;
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t16  = 1.0 / (t13 + t14 + 1);
    real_type t17  = t16 / t11;
    real_type t24  = (-t4 * t1 + t2 * xo__h) * xo__p;
    real_type result__ = -2 * t4 * t17 * (t2 * t1 + t4 * xo__h) * xo__p + 4 * (-t2 * xo__f + t4 * xo__g) * t4 * t16 / t11 / t10 * t24 + 2 * t2 * t17 * t24;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_2_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_2_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t3   = t2 * t2;
    real_type t5   = sin(xo__L);
    real_type t8   = pow(t2 * xo__f + t5 * xo__g + 1, 2);
    real_type t11  = xo__h * xo__h;
    real_type t12  = xo__k * xo__k;
    real_type result__ = 2 / (t11 + t12 + 1) / t8 * t3 * xo__p * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_2_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t11  = pow(t2 * xo__f + t4 * xo__g + 1, 2);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type result__ = -2 * t4 / (t13 + t14 + 1) / t11 * (-t2 * xo__k * xo__retrograde + t4 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_3_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t10  = t2 * xo__f + t4 * xo__g + 1;
    real_type t11  = t10 * t10;
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t19  = t4 * t4;
    real_type result__ = 4 * t19 / (t14 + t15 + 1) / t11 / t10 * (-t2 * xo__k * xo__retrograde + t4 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_3_3( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_3_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t2   = t1 * t1;
    real_type t4   = cos(xo__L);
    real_type t8   = pow(t1 * xo__g + t4 * xo__f + 1, 2);
    real_type t9   = 1.0 / t8;
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type t12  = t10 + t11 + 1;
    real_type t23  = t12 * t12;
    real_type result__ = -2 / t12 * t9 * t2 * xo__p + 4 * xo__h * t1 / t23 * t9 * (-t4 * xo__k * xo__retrograde + t1 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_3_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_3_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t8   = pow(t2 * xo__f + t5 * xo__g + 1, 2);
    real_type t9   = 1.0 / t8;
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type t12  = t10 + t11 + 1;
    real_type t24  = t12 * t12;
    real_type result__ = 2 * t5 / t12 * t9 * t2 * xo__retrograde * xo__p + 4 * xo__k * t5 / t24 * t9 * (-t2 * xo__k * xo__retrograde + t5 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_3_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_3_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__k * xo__retrograde;
    real_type t2   = sin(xo__L);
    real_type t4   = cos(xo__L);
    real_type t10  = t2 * xo__g + t4 * xo__f + 1;
    real_type t11  = t10 * t10;
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type t16  = 1.0 / (t13 + t14 + 1);
    real_type t17  = t16 / t11;
    real_type t24  = (-t4 * t1 + t2 * xo__h) * xo__p;
    real_type result__ = -2 * t2 * t17 * (t2 * t1 + t4 * xo__h) * xo__p + 4 * (-t2 * xo__f + t4 * xo__g) * t2 * t16 / t11 / t10 * t24 - 2 * t4 * t17 * t24;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_3_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_3_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t8   = pow(t2 * xo__f + t5 * xo__g + 1, 2);
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type result__ = 2 * t5 / (t10 + t11 + 1) / t8 * t2 * xo__k * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_3_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t7   = 1.0 / (t1 * xo__g + t3 * xo__f + 1);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t20  = t10 * t10;
    real_type result__ = 2 / t10 * t7 * t1 * xo__p - 4 * xo__h / t20 * t7 * (-t3 * xo__k * xo__retrograde + t1 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_4_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t7   = 1.0 / (t1 * xo__g + t3 * xo__f + 1);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t11  = t10 * t10;
    real_type t13  = 1.0 / t11 * t7;
    real_type t21  = (-t3 * xo__k * xo__retrograde + t1 * xo__h) * xo__p;
    real_type result__ = -8 * xo__h * t13 * t1 * xo__p + 16 * t8 / t11 / t10 * t7 * t21 - 4 * t13 * t21;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_4_4( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_4_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t7   = 1.0 / (t1 * xo__g + t3 * xo__f + 1);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t11  = t10 * t10;
    real_type t13  = 1.0 / t11 * t7;
    real_type result__ = -4 * xo__k * t13 * t1 * xo__p + 4 * xo__h * t13 * t3 * xo__retrograde * xo__p + 16 * xo__k * xo__h / t11 / t10 * t7 * (-t3 * xo__k * xo__retrograde + t1 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_4_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_4_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t6   = t1 * xo__f + t4 * xo__g + 1;
    real_type t7   = 1.0 / t6;
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t11  = 1.0 / t10;
    real_type t16  = t6 * t6;
    real_type t17  = 1.0 / t16;
    real_type t21  = t1 * xo__g - t4 * xo__f;
    real_type t25  = xo__k * xo__retrograde;
    real_type t30  = t10 * t10;
    real_type t31  = 1.0 / t30;
    real_type result__ = 2 * t11 * t7 * t1 * xo__p - 2 * t21 * t11 * t17 * t4 * xo__p - 4 * xo__h * t31 * t7 * (t1 * xo__h + t4 * t25) * xo__p + 4 * t21 * xo__h * t31 * t17 * (-t1 * t25 + t4 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_4_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_4_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t12  = pow(t9 + t10 + 1, 2);
    real_type result__ = 4 * xo__h / t12 / (t2 * xo__f + t5 * xo__g + 1) * t2 * xo__k * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_4_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t7   = 1.0 / (t2 * xo__f + t4 * xo__g + 1);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t21  = t11 * t11;
    real_type result__ = -2 / t11 * t7 * t2 * xo__p * xo__retrograde - 4 * xo__k / t21 * t7 * (-t2 * xo__k * xo__retrograde + t4 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_5_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t8   = 1.0 / (t2 * xo__f + t5 * xo__g + 1);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t12  = t11 * t11;
    real_type t14  = 1.0 / t12 * t8;
    real_type t22  = (-t2 * xo__k * xo__retrograde + t5 * xo__h) * xo__p;
    real_type result__ = 8 * xo__k * t14 * t2 * xo__retrograde * xo__p + 16 * t10 / t12 / t11 * t8 * t22 - 4 * t14 * t22;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_5_5( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_5_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__p * xo__retrograde;
    real_type t2   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t6   = t2 * xo__g + t3 * xo__f + 1;
    real_type t7   = 1.0 / t6;
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t12  = 1.0 / t11;
    real_type t17  = t6 * t6;
    real_type t18  = 1.0 / t17;
    real_type t22  = -t2 * xo__f + t3 * xo__g;
    real_type t26  = xo__k * xo__retrograde;
    real_type t31  = t11 * t11;
    real_type t32  = 1.0 / t31;
    real_type result__ = 2 * t12 * t7 * t2 * t1 + 2 * t22 * t12 * t18 * t3 * t1 - 4 * xo__k * t32 * t7 * (t2 * t26 + t3 * xo__h) * xo__p + 4 * t22 * xo__k * t32 * t18 * (t2 * xo__h - t26 * t3) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_5_6( p={}, f={}, g={}, h={}, k={}, L={},  retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_5_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t7   = 1.0 / (t1 * xo__f + t4 * xo__g + 1);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t17  = t10 * t10;
    real_type result__ = -2 / t10 * t7 * t1 * xo__p + 4 / t17 * t7 * t1 * t9 * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_5_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__k * xo__retrograde;
    real_type t2   = sin(xo__L);
    real_type t4   = cos(xo__L);
    real_type t10  = t2 * xo__g + t4 * xo__f + 1;
    real_type t12  = xo__h * xo__h;
    real_type t13  = xo__k * xo__k;
    real_type t15  = 1.0 / (t12 + t13 + 1);
    real_type t22  = t10 * t10;
    real_type result__ = 2 * t15 / t10 * (t2 * t1 + t4 * xo__h) * xo__p - 2 * (-t2 * xo__f + t4 * xo__g) * t15 / t22 * (-t4 * t1 + t2 * xo__h) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_6_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__k * xo__retrograde;
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t6   = t2 * t1 - t4 * xo__h;
    real_type t8   = t2 * xo__f;
    real_type t9   = t4 * xo__g;
    real_type t10  = t8 + t9 + 1;
    real_type t12  = xo__h * xo__h;
    real_type t13  = xo__k * xo__k;
    real_type t15  = 1.0 / (t12 + t13 + 1);
    real_type t23  = t10 * t10;
    real_type t25  = t15 / t23;
    real_type t28  = t2 * xo__g - t4 * xo__f;
    real_type t33  = -t6 * xo__p;
    real_type t37  = t28 * t28;
    real_type result__ = 2 * t15 / t10 * t6 * xo__p - 4 * t28 * t25 * (t4 * t1 + t2 * xo__h) * xo__p + 4 * t37 * t15 / t23 / t10 * t33 - 2 * (-t8 - t9) * t25 * t33;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_6_6( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_6_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t1   = xo__p * xo__k;
    real_type t2   = sin(xo__L);
    real_type t3   = cos(xo__L);
    real_type t6   = t2 * xo__g + t3 * xo__f + 1;
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t12  = 1.0 / (t9 + t10 + 1);
    real_type t16  = t6 * t6;
    real_type result__ = 2 * t12 / t6 * t2 * t1 + 2 * (-t2 * xo__f + t3 * xo__g) * t12 / t16 * t3 * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_6_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type result__ = -2 / (t9 + t10 + 1) / (t2 * xo__f + t4 * xo__g + 1) * t2 * xo__p * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_position__xo_D_7_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_position__xo_D_7_7( p={}, f={}, g={}, h={}, k={}, L={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = sin(xo__L);
    real_type t7   = cos(xo__L);
    real_type t18  = sqrt(xo__p);
    real_type t21  = sqrt(xo__muS);
    real_type result__ = -1.0 / (t1 + t2 + 1) * t21 / t18 * (t4 * (t1 - t2 + 1) - 2 * t7 * xo__k * xo__retrograde * xo__h - 2 * xo__f * xo__h * xo__k * xo__retrograde + t1 * xo__g - t2 * xo__g + xo__g);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = sin(xo__L);
    real_type t7   = cos(xo__L);
    real_type t18  = sqrt(xo__p);
    real_type t22  = sqrt(xo__muS);
    real_type result__ = 1.0 / (t1 + t2 + 1) * t22 / t18 / xo__p * (t4 * (t1 - t2 + 1) - 2 * t7 * xo__k * xo__retrograde * xo__h - 2 * xo__f * xo__h * xo__k * xo__retrograde + t1 * xo__g - t2 * xo__g + xo__g) / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_1( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_1_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = sin(xo__L);
    real_type t7   = cos(xo__L);
    real_type t18  = xo__p * xo__p;
    real_type t19  = sqrt(xo__p);
    real_type t23  = sqrt(xo__muS);
    real_type result__ = -3.0 / 4.0 / (t1 + t2 + 1) * t23 / t19 / t18 * (t4 * (t1 - t2 + 1) - 2 * t7 * xo__k * xo__retrograde * xo__h - 2 * xo__f * xo__h * xo__k * xo__retrograde + t1 * xo__g - t2 * xo__g + xo__g);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_1_1( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_1_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t3   = sqrt(xo__p);
    real_type t6   = sqrt(xo__muS);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type result__ = -1.0 / (t8 + t9 + 1) * t6 / t3 / xo__p * xo__h * xo__k * xo__retrograde;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_1_2( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_1_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = sqrt(xo__p);
    real_type t8   = sqrt(xo__muS);
    real_type result__ = 1.0 / (t1 + t2 + 1) * t8 / t4 / xo__p * (t1 - t2 + 1) / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_1_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_1_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = xo__k * xo__retrograde;
    real_type t4   = cos(xo__L);
    real_type t11  = sqrt(xo__p);
    real_type t13  = 1.0 / t11 / xo__p;
    real_type t15  = sqrt(xo__muS);
    real_type t16  = xo__h * xo__h;
    real_type t17  = xo__k * xo__k;
    real_type t18  = t16 + t17 + 1;
    real_type t36  = t18 * t18;
    real_type result__ = 1.0 / t18 * t15 * t13 * (-2 * xo__f * xo__k * xo__retrograde + 2 * t1 * xo__h - 2 * t4 * t3 + 2 * xo__g * xo__h) / 2 - xo__h / t36 * t15 * t13 * (t1 * (t16 - t17 + 1) - 2 * t4 * xo__k * xo__retrograde * xo__h - 2 * t3 * xo__f * xo__h + t16 * xo__g - t17 * xo__g + xo__g);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_1_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_1_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = xo__h * xo__retrograde;
    real_type t4   = cos(xo__L);
    real_type t6   = xo__f * xo__h;
    real_type t11  = sqrt(xo__p);
    real_type t13  = 1.0 / t11 / xo__p;
    real_type t15  = sqrt(xo__muS);
    real_type t16  = xo__h * xo__h;
    real_type t17  = xo__k * xo__k;
    real_type t18  = t16 + t17 + 1;
    real_type t35  = t18 * t18;
    real_type result__ = 1.0 / t18 * t15 * t13 * (-2 * t1 * xo__k - 2 * t4 * t3 - 2 * xo__retrograde * t6 - 2 * xo__g * xo__k) / 2 - xo__k / t35 * t15 * t13 * (t1 * (t16 - t17 + 1) - 2 * t4 * xo__k * t3 - 2 * xo__k * xo__retrograde * t6 + t16 * xo__g - t17 * xo__g + xo__g);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_1_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_1_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type t12  = sqrt(xo__p);
    real_type t16  = sqrt(xo__muS);
    real_type result__ = 1.0 / (t1 + t2 + 1) * t16 / t12 / xo__p * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h) / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_1_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_1_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = sin(xo__L);
    real_type t7   = cos(xo__L);
    real_type t18  = sqrt(xo__p);
    real_type t22  = sqrt(xo__muS);
    real_type result__ = 1.0 / (t1 + t2 + 1) / t22 / t18 / xo__p * (t4 * (t1 - t2 + 1) - 2 * t7 * xo__k * xo__retrograde * xo__h - 2 * xo__f * xo__h * xo__k * xo__retrograde + t1 * xo__g - t2 * xo__g + xo__g) / 4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_1_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_1_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t8   = sqrt(xo__p);
    real_type t12  = sqrt(xo__muS);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type result__ = 1.0 / (t13 + t14 + 1) * t12 / t8 / xo__p * (-2 * t2 * xo__k * xo__h - 2 * xo__f * xo__h * xo__k) / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_1_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t3   = sqrt(xo__p);
    real_type t5   = sqrt(xo__muS);
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type result__ = 2 / (t7 + t8 + 1) * t5 / t3 * xo__h * xo__k * xo__retrograde;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_2( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_2_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_2_2( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_2_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_2_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_2_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t5   = t4 / t2;
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t8   = t6 + t7 + 1;
    real_type t15  = t8 * t8;
    real_type result__ = 2 / t8 * t5 * xo__k * xo__retrograde - 4 / t15 * t5 * xo__retrograde * xo__k * t6;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_2_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_2_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t5   = t4 / t2;
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t8   = t6 + t7 + 1;
    real_type t15  = t8 * t8;
    real_type result__ = 2 / t8 * t5 * xo__h * xo__retrograde - 4 / t15 * t5 * xo__retrograde * t7 * xo__h;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_2_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_2_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_2_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_2_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t3   = sqrt(xo__p);
    real_type t5   = sqrt(xo__muS);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type result__ = 1.0 / (t8 + t9 + 1) / t5 / t3 * xo__h * xo__k * xo__retrograde;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_2_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_2_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type result__ = 2 / (t6 + t7 + 1) * t4 / t2 * xo__h * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_2_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = sqrt(xo__p);
    real_type t7   = sqrt(xo__muS);
    real_type result__ = -1.0 / (t1 + t2 + 1) * t7 / t4 * (t1 - t2 + 1);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_3_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_3_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_3_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t2   = 1.0 / t1;
    real_type t4   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t7   = t5 + t6 + 1;
    real_type t13  = t7 * t7;
    real_type result__ = -2 / t7 * t4 * t2 * xo__h + 2 * xo__h / t13 * t4 * t2 * (t5 - t6 + 1);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_3_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_3_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t2   = 1.0 / t1;
    real_type t4   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t7   = t5 + t6 + 1;
    real_type t13  = t7 * t7;
    real_type result__ = 2 / t7 * t4 * t2 * xo__k + 2 * xo__k / t13 * t4 * t2 * (t5 - t6 + 1);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_3_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_3_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_3_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_3_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = sqrt(xo__p);
    real_type t7   = sqrt(xo__muS);
    real_type result__ = -1.0 / (t1 + t2 + 1) / t7 / t4 * (t1 - t2 + 1) / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_3_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_3_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_3_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = xo__k * xo__retrograde;
    real_type t4   = cos(xo__L);
    real_type t11  = sqrt(xo__p);
    real_type t12  = 1.0 / t11;
    real_type t14  = sqrt(xo__muS);
    real_type t15  = xo__h * xo__h;
    real_type t16  = xo__k * xo__k;
    real_type t17  = t15 + t16 + 1;
    real_type t34  = t17 * t17;
    real_type result__ = -1.0 / t17 * t14 * t12 * (-2 * xo__f * xo__k * xo__retrograde + 2 * t1 * xo__h - 2 * t4 * t3 + 2 * xo__g * xo__h) + 2 * xo__h / t34 * t14 * t12 * (t1 * (t15 - t16 + 1) - 2 * t4 * xo__k * xo__retrograde * xo__h - 2 * t3 * xo__f * xo__h + t15 * xo__g - t16 * xo__g + xo__g);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_4_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t4   = sqrt(xo__p);
    real_type t5   = 1.0 / t4;
    real_type t7   = sqrt(xo__muS);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t15  = xo__k * xo__retrograde;
    real_type t16  = cos(xo__L);
    real_type t24  = t10 * t10;
    real_type t26  = 1.0 / t24 * t7;
    real_type t42  = t5 * (t1 * (t8 - t9 + 1) - 2 * t16 * xo__k * xo__retrograde * xo__h - 2 * t15 * xo__f * xo__h + t8 * xo__g - t9 * xo__g + xo__g);
    real_type result__ = -1.0 / t10 * t7 * t5 * (2 * t1 + 2 * xo__g) + 4 * xo__h * t26 * t5 * (-2 * xo__f * xo__k * xo__retrograde + 2 * t1 * xo__h - 2 * t16 * t15 + 2 * xo__g * xo__h) - 8 * t8 / t24 / t10 * t7 * t42 + 2 * t26 * t42;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_4_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_4_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t6   = sqrt(xo__p);
    real_type t7   = 1.0 / t6;
    real_type t9   = sqrt(xo__muS);
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type t12  = t10 + t11 + 1;
    real_type t16  = sin(xo__L);
    real_type t18  = xo__k * xo__retrograde;
    real_type t26  = t12 * t12;
    real_type t28  = 1.0 / t26 * t9;
    real_type t33  = xo__h * xo__retrograde;
    real_type t35  = xo__f * xo__h;
    real_type result__ = -1.0 / t12 * t9 * t7 * (-2 * t1 * xo__retrograde - 2 * xo__f * xo__retrograde) + 2 * xo__k * t28 * t7 * (-2 * xo__f * xo__k * xo__retrograde - 2 * t1 * t18 + 2 * t16 * xo__h + 2 * xo__g * xo__h) + 2 * xo__h * t28 * t7 * (-2 * t1 * t33 - 2 * t16 * xo__k - 2 * xo__retrograde * t35 - 2 * xo__g * xo__k) - 8 * xo__k * xo__h / t26 / t12 * t9 * t7 * (t16 * (t10 - t11 + 1) - 2 * t1 * xo__k * t33 - 2 * t18 * t35 + t10 * xo__g - t11 * xo__g + xo__g);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_4_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_4_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t4   = sin(xo__L);
    real_type t8   = sqrt(xo__p);
    real_type t9   = 1.0 / t8;
    real_type t11  = sqrt(xo__muS);
    real_type t12  = xo__h * xo__h;
    real_type t13  = xo__k * xo__k;
    real_type t14  = t12 + t13 + 1;
    real_type t26  = t14 * t14;
    real_type result__ = -1.0 / t14 * t11 * t9 * (2 * t4 * xo__k * xo__retrograde + 2 * t1 * xo__h) + 2 * xo__h / t26 * t11 * t9 * (t1 * (t12 - t13 + 1) + 2 * t4 * xo__k * xo__retrograde * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_4_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_4_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = xo__k * xo__retrograde;
    real_type t4   = cos(xo__L);
    real_type t11  = sqrt(xo__p);
    real_type t12  = 1.0 / t11;
    real_type t14  = sqrt(xo__muS);
    real_type t15  = 1.0 / t14;
    real_type t16  = xo__h * xo__h;
    real_type t17  = xo__k * xo__k;
    real_type t18  = t16 + t17 + 1;
    real_type t36  = t18 * t18;
    real_type result__ = -1.0 / t18 * t15 * t12 * (-2 * xo__f * xo__k * xo__retrograde + 2 * t1 * xo__h - 2 * t4 * t3 + 2 * xo__g * xo__h) / 2 + xo__h / t36 * t15 * t12 * (t1 * (t16 - t17 + 1) - 2 * t4 * xo__k * xo__retrograde * xo__h - 2 * t3 * xo__f * xo__h + t16 * xo__g - t17 * xo__g + xo__g);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_4_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_4_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t6   = sqrt(xo__p);
    real_type t7   = 1.0 / t6;
    real_type t9   = sqrt(xo__muS);
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type t12  = t10 + t11 + 1;
    real_type t23  = t12 * t12;
    real_type result__ = -1.0 / t12 * t9 * t7 * (-2 * t1 * xo__k - 2 * xo__f * xo__k) + 2 * xo__h / t23 * t9 * t7 * (-2 * t1 * xo__k * xo__h - 2 * xo__f * xo__h * xo__k);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_4_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = xo__h * xo__retrograde;
    real_type t4   = cos(xo__L);
    real_type t6   = xo__f * xo__h;
    real_type t11  = sqrt(xo__p);
    real_type t12  = 1.0 / t11;
    real_type t14  = sqrt(xo__muS);
    real_type t15  = xo__h * xo__h;
    real_type t16  = xo__k * xo__k;
    real_type t17  = t15 + t16 + 1;
    real_type t33  = t17 * t17;
    real_type result__ = -1.0 / t17 * t14 * t12 * (-2 * t1 * xo__k - 2 * t4 * t3 - 2 * xo__retrograde * t6 - 2 * xo__g * xo__k) + 2 * xo__k / t33 * t14 * t12 * (t1 * (t15 - t16 + 1) - 2 * t4 * xo__k * t3 - 2 * xo__k * xo__retrograde * t6 + t15 * xo__g - t16 * xo__g + xo__g);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_5_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t4   = sqrt(xo__p);
    real_type t5   = 1.0 / t4;
    real_type t7   = sqrt(xo__muS);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t15  = xo__h * xo__retrograde;
    real_type t16  = cos(xo__L);
    real_type t18  = xo__f * xo__h;
    real_type t24  = t10 * t10;
    real_type t26  = 1.0 / t24 * t7;
    real_type t41  = t5 * (t1 * (t8 - t9 + 1) - 2 * t16 * xo__k * t15 - 2 * xo__k * xo__retrograde * t18 + t8 * xo__g - t9 * xo__g + xo__g);
    real_type result__ = -1.0 / t10 * t7 * t5 * (-2 * t1 - 2 * xo__g) + 4 * xo__k * t26 * t5 * (-2 * t1 * xo__k - 2 * t16 * t15 - 2 * xo__retrograde * t18 - 2 * xo__g * xo__k) - 8 * t9 / t24 / t10 * t7 * t41 + 2 * t26 * t41;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_5_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_5_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t3   = xo__h * xo__retrograde;
    real_type t4   = sin(xo__L);
    real_type t8   = sqrt(xo__p);
    real_type t9   = 1.0 / t8;
    real_type t11  = sqrt(xo__muS);
    real_type t12  = xo__h * xo__h;
    real_type t13  = xo__k * xo__k;
    real_type t14  = t12 + t13 + 1;
    real_type t25  = t14 * t14;
    real_type result__ = -1.0 / t14 * t11 * t9 * (-2 * t1 * xo__k + 2 * t4 * t3) + 2 * xo__k / t25 * t11 * t9 * (t1 * (t12 - t13 + 1) + 2 * t4 * xo__k * t3);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_5_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_5_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sin(xo__L);
    real_type t3   = xo__h * xo__retrograde;
    real_type t4   = cos(xo__L);
    real_type t6   = xo__f * xo__h;
    real_type t11  = sqrt(xo__p);
    real_type t12  = 1.0 / t11;
    real_type t14  = sqrt(xo__muS);
    real_type t15  = 1.0 / t14;
    real_type t16  = xo__h * xo__h;
    real_type t17  = xo__k * xo__k;
    real_type t18  = t16 + t17 + 1;
    real_type t35  = t18 * t18;
    real_type result__ = -1.0 / t18 * t15 * t12 * (-2 * t1 * xo__k - 2 * t4 * t3 - 2 * xo__retrograde * t6 - 2 * xo__g * xo__k) / 2 + xo__k / t35 * t15 * t12 * (t1 * (t16 - t17 + 1) - 2 * t4 * xo__k * t3 - 2 * xo__k * xo__retrograde * t6 + t16 * xo__g - t17 * xo__g + xo__g);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_5_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_5_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = cos(xo__L);
    real_type t3   = xo__f * xo__h;
    real_type t6   = sqrt(xo__p);
    real_type t7   = 1.0 / t6;
    real_type t9   = sqrt(xo__muS);
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type t12  = t10 + t11 + 1;
    real_type t22  = t12 * t12;
    real_type result__ = -1.0 / t12 * t9 * t7 * (-2 * t1 * xo__h - 2 * t3) + 2 * xo__k / t22 * t9 * t7 * (-2 * t1 * xo__k * xo__h - 2 * xo__k * t3);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_5_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type t12  = sqrt(xo__p);
    real_type t15  = sqrt(xo__muS);
    real_type result__ = -1.0 / (t1 + t2 + 1) * t15 / t12 * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_6_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = sin(xo__L);
    real_type t7   = cos(xo__L);
    real_type t12  = sqrt(xo__p);
    real_type t15  = sqrt(xo__muS);
    real_type result__ = -1.0 / (t1 + t2 + 1) * t15 / t12 * (-t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_6_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_6_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = cos(xo__L);
    real_type t7   = sin(xo__L);
    real_type t12  = sqrt(xo__p);
    real_type t15  = sqrt(xo__muS);
    real_type result__ = -1.0 / (t1 + t2 + 1) / t15 / t12 * (t4 * (t1 - t2 + 1) + 2 * t7 * xo__k * xo__retrograde * xo__h) / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_6_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_6_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = sin(xo__L);
    real_type t4   = sqrt(xo__p);
    real_type t6   = sqrt(xo__muS);
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type result__ = -2 / (t8 + t9 + 1) * t6 / t4 * t2 * xo__k * xo__h;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_6_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = sin(xo__L);
    real_type t7   = cos(xo__L);
    real_type t18  = sqrt(xo__p);
    real_type t21  = sqrt(xo__muS);
    real_type result__ = -1.0 / (t1 + t2 + 1) / t21 / t18 * (t4 * (t1 - t2 + 1) - 2 * t7 * xo__k * xo__retrograde * xo__h - 2 * xo__f * xo__h * xo__k * xo__retrograde + t1 * xo__g - t2 * xo__g + xo__g) / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_7_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__h * xo__h;
    real_type t2   = xo__k * xo__k;
    real_type t4   = sin(xo__L);
    real_type t7   = cos(xo__L);
    real_type t18  = sqrt(xo__p);
    real_type t21  = sqrt(xo__muS);
    real_type result__ = 1.0 / (t1 + t2 + 1) / t21 / xo__muS / t18 * (t4 * (t1 - t2 + 1) - 2 * t7 * xo__k * xo__retrograde * xo__h - 2 * xo__f * xo__h * xo__k * xo__retrograde + t1 * xo__g - t2 * xo__g + xo__g) / 4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_7_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_7_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t8   = sqrt(xo__p);
    real_type t11  = sqrt(xo__muS);
    real_type t13  = xo__h * xo__h;
    real_type t14  = xo__k * xo__k;
    real_type result__ = -1.0 / (t13 + t14 + 1) / t11 / t8 * (-2 * t2 * xo__k * xo__h - 2 * xo__f * xo__h * xo__k) / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_7_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = cos(xo__L);
    real_type t8   = sqrt(xo__p);
    real_type t11  = sqrt(xo__muS);
    real_type t12  = xo__h * xo__h;
    real_type t13  = xo__k * xo__k;
    real_type result__ = -1.0 / (t12 + t13 + 1) * t11 / t8 * (-2 * t2 * xo__k * xo__h - 2 * xo__f * xo__h * xo__k);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_x_velocity__xo_D_8_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_x_velocity__xo_D_8_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t7   = t5 - t6 - 1;
    real_type t9   = cos(xo__L);
    real_type t13  = sin(xo__L);
    real_type result__ = -2 / (t5 + t6 + 1) * (t9 * t7 * xo__retrograde / 2 + t13 * xo__k * xo__h + xo__retrograde * t7 * xo__f / 2 + xo__g * xo__h * xo__k) * t3 / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t8   = t6 - t7 - 1;
    real_type t10  = cos(xo__L);
    real_type t14  = sin(xo__L);
    real_type result__ = 1.0 / (t6 + t7 + 1) * (t10 * t8 * xo__retrograde / 2 + t14 * xo__k * xo__h + xo__retrograde * t8 * xo__f / 2 + xo__g * xo__h * xo__k) * t4 / t1 / xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_1( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_1_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = xo__p * xo__p;
    real_type t2   = sqrt(xo__p);
    real_type t5   = sqrt(xo__muS);
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type t9   = t7 - t8 - 1;
    real_type t11  = cos(xo__L);
    real_type t15  = sin(xo__L);
    real_type result__ = -3.0 / 2.0 / (t7 + t8 + 1) * (t11 * t9 * xo__retrograde / 2 + t15 * xo__k * xo__h + xo__retrograde * t9 * xo__f / 2 + xo__g * xo__h * xo__k) * t5 / t2 / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_1_1( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_1_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type result__ = 1.0 / (t6 + t7 + 1) * xo__retrograde * (t6 - t7 - 1) * t4 / t1 / xo__p / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_1_2( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_1_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type result__ = 1.0 / (t7 + t8 + 1) * xo__k * xo__h * t4 / t1 / xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_1_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_1_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t5   = t4 / t1 / xo__p;
    real_type t7   = cos(xo__L);
    real_type t9   = sin(xo__L);
    real_type t15  = xo__h * xo__h;
    real_type t16  = xo__k * xo__k;
    real_type t17  = t15 + t16 + 1;
    real_type t21  = t15 - t16 - 1;
    real_type t33  = t17 * t17;
    real_type result__ = 1.0 / t17 * (t7 * xo__h * xo__retrograde + xo__f * xo__h * xo__retrograde + t9 * xo__k + xo__g * xo__k) * t5 - 2 * xo__h / t33 * (t7 * t21 * xo__retrograde / 2 + t9 * xo__k * xo__h + xo__retrograde * t21 * xo__f / 2 + xo__g * xo__h * xo__k) * t5;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_1_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_1_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t5   = t4 / t1 / xo__p;
    real_type t7   = cos(xo__L);
    real_type t9   = sin(xo__L);
    real_type t13  = xo__g * xo__h;
    real_type t15  = xo__h * xo__h;
    real_type t16  = xo__k * xo__k;
    real_type t17  = t15 + t16 + 1;
    real_type t21  = t15 - t16 - 1;
    real_type t32  = t17 * t17;
    real_type result__ = 1.0 / t17 * (-t7 * xo__k * xo__retrograde - xo__f * xo__k * xo__retrograde + t9 * xo__h + t13) * t5 - 2 * xo__k / t32 * (t7 * t21 * xo__retrograde / 2 + t9 * xo__k * xo__h + xo__retrograde * t21 * xo__f / 2 + xo__k * t13) * t5;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_1_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_1_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t10  = sin(xo__L);
    real_type t14  = cos(xo__L);
    real_type result__ = 1.0 / (t6 + t7 + 1) * (-t10 * (t6 - t7 - 1) * xo__retrograde / 2 + t14 * xo__k * xo__h) * t4 / t1 / xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_1_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_1_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type t9   = t7 - t8 - 1;
    real_type t11  = cos(xo__L);
    real_type t15  = sin(xo__L);
    real_type result__ = 1.0 / (t7 + t8 + 1) * (t11 * t9 * xo__retrograde / 2 + t15 * xo__k * xo__h + xo__retrograde * t9 * xo__f / 2 + xo__g * xo__h * xo__k) / t4 / t1 / xo__p / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_1_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_1_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t8   = t6 - t7 - 1;
    real_type t9   = cos(xo__L);
    real_type result__ = 1.0 / (t6 + t7 + 1) * (t9 * t8 / 2 + t8 * xo__f / 2) * t4 / t1 / xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_1_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type result__ = -1.0 / (t5 + t6 + 1) * xo__retrograde * (t5 - t6 - 1) * t3 / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_2( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_2_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_2_2( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_2_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_2_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_2_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t5   = t4 / t2;
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t8   = t6 + t7 + 1;
    real_type t14  = t8 * t8;
    real_type result__ = -2 / t8 * t5 * xo__h * xo__retrograde + 2 * xo__h / t14 * xo__retrograde * (t6 - t7 - 1) * t5;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_2_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_2_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t5   = t4 / t2;
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t8   = t6 + t7 + 1;
    real_type t14  = t8 * t8;
    real_type result__ = 2 / t8 * t5 * xo__k * xo__retrograde + 2 * xo__k / t14 * xo__retrograde * (t6 - t7 - 1) * t5;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_2_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_2_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_2_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_2_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type result__ = -1.0 / (t6 + t7 + 1) * xo__retrograde * (t6 - t7 - 1) / t3 / t1 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_2_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_2_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type result__ = -1.0 / (t5 + t6 + 1) * (t5 - t6 - 1) * t3 / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_2_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type result__ = -2 / (t6 + t7 + 1) * t4 / t2 * xo__h * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_3_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_3_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_3_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t2   = 1.0 / t1;
    real_type t4   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t7   = t5 + t6 + 1;
    real_type t14  = t7 * t7;
    real_type result__ = -2 / t7 * t4 * t2 * xo__k + 4 / t14 * t4 * t2 * xo__k * t5;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_3_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_3_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t2   = 1.0 / t1;
    real_type t4   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t7   = t5 + t6 + 1;
    real_type t14  = t7 * t7;
    real_type result__ = -2 / t7 * t4 * t2 * xo__h + 4 / t14 * t4 * t2 * t6 * xo__h;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_3_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_3_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_3_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_3_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type result__ = -1.0 / (t7 + t8 + 1) / t4 / t2 * xo__h * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_3_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_3_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_3_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t4   = t3 / t1;
    real_type t6   = cos(xo__L);
    real_type t8   = sin(xo__L);
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t16  = t14 + t15 + 1;
    real_type t21  = t14 - t15 - 1;
    real_type t33  = t16 * t16;
    real_type result__ = -2 / t16 * (t6 * xo__h * xo__retrograde + xo__f * xo__h * xo__retrograde + t8 * xo__k + xo__g * xo__k) * t4 + 4 * xo__h / t33 * (t6 * t21 * xo__retrograde / 2 + t8 * xo__k * xo__h + xo__retrograde * t21 * xo__f / 2 + xo__g * xo__h * xo__k) * t4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_4_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t4   = t3 / t1;
    real_type t5   = cos(xo__L);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t18  = sin(xo__L);
    real_type t24  = t11 * t11;
    real_type t25  = 1.0 / t24;
    real_type t30  = t9 - t10 - 1;
    real_type t41  = t5 * t30 * xo__retrograde / 2 + t18 * xo__k * xo__h + xo__retrograde * t30 * xo__f / 2 + xo__g * xo__h * xo__k;
    real_type result__ = -2 / t11 * (t5 * xo__retrograde + xo__f * xo__retrograde) * t4 + 8 * xo__h * t25 * (t5 * xo__h * xo__retrograde + xo__f * xo__h * xo__retrograde + t18 * xo__k + xo__g * xo__k) * t4 - 16 * t9 / t24 / t11 * t41 * t4 + 4 * t25 * t41 * t4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_4_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_4_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t4   = t3 / t1;
    real_type t5   = sin(xo__L);
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type t9   = t7 + t8 + 1;
    real_type t15  = cos(xo__L);
    real_type t22  = t9 * t9;
    real_type t23  = 1.0 / t22;
    real_type t33  = xo__g * xo__h;
    real_type t39  = t7 - t8 - 1;
    real_type result__ = -2 / t9 * (t5 + xo__g) * t4 + 4 * xo__k * t23 * (t15 * xo__h * xo__retrograde + xo__f * xo__h * xo__retrograde + t5 * xo__k + xo__g * xo__k) * t4 + 4 * xo__h * t23 * (-t15 * xo__k * xo__retrograde - xo__f * xo__k * xo__retrograde + t5 * xo__h + t33) * t4 - 16 * xo__k * xo__h / t22 / t9 * (t15 * t39 * xo__retrograde / 2 + t5 * xo__k * xo__h + xo__retrograde * t39 * xo__f / 2 + xo__k * t33) * t4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_4_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_4_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t4   = t3 / t1;
    real_type t6   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type t11  = xo__h * xo__h;
    real_type t12  = xo__k * xo__k;
    real_type t13  = t11 + t12 + 1;
    real_type t25  = t13 * t13;
    real_type result__ = -2 / t13 * (-t6 * xo__h * xo__retrograde + t8 * xo__k) * t4 + 4 * xo__h / t25 * (-t6 * (t11 - t12 - 1) * xo__retrograde / 2 + t8 * xo__k * xo__h) * t4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_4_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_4_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t5   = 1.0 / t3 / t1;
    real_type t7   = cos(xo__L);
    real_type t9   = sin(xo__L);
    real_type t15  = xo__h * xo__h;
    real_type t16  = xo__k * xo__k;
    real_type t17  = t15 + t16 + 1;
    real_type t21  = t15 - t16 - 1;
    real_type t33  = t17 * t17;
    real_type result__ = -1.0 / t17 * (t7 * xo__h * xo__retrograde + xo__f * xo__h * xo__retrograde + t9 * xo__k + xo__g * xo__k) * t5 + 2 * xo__h / t33 * (t7 * t21 * xo__retrograde / 2 + t9 * xo__k * xo__h + xo__retrograde * t21 * xo__f / 2 + xo__g * xo__h * xo__k) * t5;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_4_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_4_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t4   = t3 / t1;
    real_type t5   = cos(xo__L);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t16  = t9 - t10 - 1;
    real_type t21  = t11 * t11;
    real_type result__ = -2 / t11 * (t5 * xo__h + xo__f * xo__h) * t4 + 4 * xo__h / t21 * (t5 * t16 / 2 + t16 * xo__f / 2) * t4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_4_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t4   = t3 / t1;
    real_type t6   = cos(xo__L);
    real_type t8   = sin(xo__L);
    real_type t12  = xo__g * xo__h;
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type t16  = t14 + t15 + 1;
    real_type t21  = t14 - t15 - 1;
    real_type t32  = t16 * t16;
    real_type result__ = -2 / t16 * (-t6 * xo__k * xo__retrograde - xo__f * xo__k * xo__retrograde + t8 * xo__h + t12) * t4 + 4 * xo__k / t32 * (t6 * t21 * xo__retrograde / 2 + t8 * xo__k * xo__h + xo__retrograde * t21 * xo__f / 2 + xo__k * t12) * t4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_5_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t4   = t3 / t1;
    real_type t5   = cos(xo__L);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t18  = sin(xo__L);
    real_type t22  = xo__g * xo__h;
    real_type t24  = t11 * t11;
    real_type t25  = 1.0 / t24;
    real_type t30  = t9 - t10 - 1;
    real_type t40  = t5 * t30 * xo__retrograde / 2 + t18 * xo__k * xo__h + xo__retrograde * t30 * xo__f / 2 + xo__k * t22;
    real_type result__ = -2 / t11 * (-t5 * xo__retrograde - xo__f * xo__retrograde) * t4 + 8 * xo__k * t25 * (-t5 * xo__k * xo__retrograde - xo__f * xo__k * xo__retrograde + t18 * xo__h + t22) * t4 - 16 * t10 / t24 / t11 * t40 * t4 + 4 * t25 * t40 * t4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_5_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_5_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t4   = t3 / t1;
    real_type t6   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type t11  = xo__h * xo__h;
    real_type t12  = xo__k * xo__k;
    real_type t13  = t11 + t12 + 1;
    real_type t25  = t13 * t13;
    real_type result__ = -2 / t13 * (t6 * xo__k * xo__retrograde + t8 * xo__h) * t4 + 4 * xo__k / t25 * (-t6 * (t11 - t12 - 1) * xo__retrograde / 2 + t8 * xo__k * xo__h) * t4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_5_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_5_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t5   = 1.0 / t3 / t1;
    real_type t7   = cos(xo__L);
    real_type t9   = sin(xo__L);
    real_type t13  = xo__g * xo__h;
    real_type t15  = xo__h * xo__h;
    real_type t16  = xo__k * xo__k;
    real_type t17  = t15 + t16 + 1;
    real_type t21  = t15 - t16 - 1;
    real_type t32  = t17 * t17;
    real_type result__ = -1.0 / t17 * (-t7 * xo__k * xo__retrograde - xo__f * xo__k * xo__retrograde + t9 * xo__h + t13) * t5 + 2 * xo__k / t32 * (t7 * t21 * xo__retrograde / 2 + t9 * xo__k * xo__h + xo__retrograde * t21 * xo__f / 2 + xo__k * t13) * t5;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_5_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_5_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t4   = t3 / t1;
    real_type t5   = cos(xo__L);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t16  = t9 - t10 - 1;
    real_type t21  = t11 * t11;
    real_type result__ = -2 / t11 * (-t5 * xo__k - xo__f * xo__k) * t4 + 4 * xo__k / t21 * (t5 * t16 / 2 + t16 * xo__f / 2) * t4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_5_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t9   = sin(xo__L);
    real_type t13  = cos(xo__L);
    real_type result__ = -2 / (t5 + t6 + 1) * (-t9 * (t5 - t6 - 1) * xo__retrograde / 2 + t13 * xo__k * xo__h) * t3 / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_6_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t9   = cos(xo__L);
    real_type t13  = sin(xo__L);
    real_type result__ = -2 / (t5 + t6 + 1) * (-t9 * (t5 - t6 - 1) * xo__retrograde / 2 - t13 * xo__k * xo__h) * t3 / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_6_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_6_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t10  = sin(xo__L);
    real_type t14  = cos(xo__L);
    real_type result__ = -1.0 / (t6 + t7 + 1) * (-t10 * (t6 - t7 - 1) * xo__retrograde / 2 + t14 * xo__k * xo__h) / t3 / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_6_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_6_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t8   = sin(xo__L);
    real_type result__ = 1.0 / (t5 + t6 + 1) * t8 * (t5 - t6 - 1) * t3 / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_6_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t8   = t6 - t7 - 1;
    real_type t10  = cos(xo__L);
    real_type t14  = sin(xo__L);
    real_type result__ = -1.0 / (t6 + t7 + 1) * (t10 * t8 * xo__retrograde / 2 + t14 * xo__k * xo__h + xo__retrograde * t8 * xo__f / 2 + xo__g * xo__h * xo__k) / t3 / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_7_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type t9   = t7 - t8 - 1;
    real_type t11  = cos(xo__L);
    real_type t15  = sin(xo__L);
    real_type result__ = 1.0 / (t7 + t8 + 1) * (t11 * t9 * xo__retrograde / 2 + t15 * xo__k * xo__h + xo__retrograde * t9 * xo__f / 2 + xo__g * xo__h * xo__k) / t3 / xo__muS / t1 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_7_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_7_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t8   = t6 - t7 - 1;
    real_type t9   = cos(xo__L);
    real_type result__ = -1.0 / (t6 + t7 + 1) * (t9 * t8 / 2 + t8 * xo__f / 2) / t3 / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_7_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t3   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t7   = t5 - t6 - 1;
    real_type t8   = cos(xo__L);
    real_type result__ = -2 / (t5 + t6 + 1) * (t8 * t7 / 2 + t7 * xo__f / 2) * t3 / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_y_velocity__xo_D_8_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_y_velocity__xo_D_8_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t3   = sin(xo__L);
    real_type t7   = cos(xo__L);
    real_type t12  = sqrt(xo__p);
    real_type t14  = xo__h * xo__h;
    real_type t15  = xo__k * xo__k;
    real_type result__ = 2 / (t14 + t15 + 1) / t12 * (t3 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t7 * xo__h + xo__f * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t3   = sin(xo__L);
    real_type t7   = cos(xo__L);
    real_type t12  = sqrt(xo__p);
    real_type t15  = xo__h * xo__h;
    real_type t16  = xo__k * xo__k;
    real_type result__ = -1.0 / (t15 + t16 + 1) / t12 / xo__p * (t3 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t7 * xo__h + xo__f * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_1( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_1_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t3   = sin(xo__L);
    real_type t7   = cos(xo__L);
    real_type t12  = xo__p * xo__p;
    real_type t13  = sqrt(xo__p);
    real_type t16  = xo__h * xo__h;
    real_type t17  = xo__k * xo__k;
    real_type result__ = 3.0 / 2.0 / (t16 + t17 + 1) / t13 / t12 * (t3 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t7 * xo__h + xo__f * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_1_1( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_1_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t3   = sqrt(xo__p);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type result__ = -1.0 / (t6 + t7 + 1) / t3 / xo__p * xo__h * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_1_2( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_1_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t3   = sqrt(xo__p);
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type result__ = -1.0 / (t7 + t8 + 1) / t3 / xo__p * xo__retrograde * xo__k * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_1_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_1_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = cos(xo__L);
    real_type t5   = sqrt(xo__p);
    real_type t7   = 1.0 / t5 / xo__p;
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t15  = sin(xo__L);
    real_type t23  = t10 * t10;
    real_type result__ = -1.0 / t10 * t7 * (t2 + xo__f) * t1 + 2 * xo__h / t23 * t7 * (t15 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t2 * xo__h + xo__f * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_1_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_1_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = sin(xo__L);
    real_type t7   = sqrt(xo__p);
    real_type t9   = 1.0 / t7 / xo__p;
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type t12  = t10 + t11 + 1;
    real_type t20  = cos(xo__L);
    real_type t25  = t12 * t12;
    real_type result__ = -1.0 / t12 * t9 * (t2 * xo__retrograde + xo__g * xo__retrograde) * t1 + 2 * xo__k / t25 * t9 * (t2 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t20 * xo__h + xo__f * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_1_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_1_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t3   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = sqrt(xo__p);
    real_type t12  = xo__h * xo__h;
    real_type t13  = xo__k * xo__k;
    real_type result__ = -1.0 / (t12 + t13 + 1) / t9 / xo__p * (t3 * xo__k * xo__retrograde - t5 * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_1_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_1_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t4   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type t13  = sqrt(xo__p);
    real_type t16  = xo__h * xo__h;
    real_type t17  = xo__k * xo__k;
    real_type result__ = -1.0 / (t16 + t17 + 1) / t13 / xo__p * (t4 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t8 * xo__h + xo__f * xo__h) / t1 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_1_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_1_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = sin(xo__L);
    real_type t7   = sqrt(xo__p);
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type result__ = -1.0 / (t10 + t11 + 1) / t7 / xo__p * (t2 * xo__k + xo__g * xo__k) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_1_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type result__ = 2 / (t5 + t6 + 1) * t4 / t1 * xo__h;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_2( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_2_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_2_2( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_2_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_2_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_2_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t2   = 1.0 / t1;
    real_type t3   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t7   = t5 + t6 + 1;
    real_type t12  = t7 * t7;
    real_type result__ = 2 / t7 * t3 * t2 - 4 / t12 * t3 * t2 * t5;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_2_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_2_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t8   = pow(t5 + t6 + 1, 2);
    real_type result__ = -4 * xo__k / t8 * t4 / t1 * xo__h;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_2_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_2_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_2_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_2_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type result__ = 1.0 / (t6 + t7 + 1) / t4 / t1 * xo__h;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_2_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_2_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_2_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type result__ = 2 / (t6 + t7 + 1) * t4 / t2 * xo__k * xo__retrograde;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_3_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_3_3( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_3_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = sqrt(xo__p);
    real_type t5   = sqrt(xo__muS);
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t9   = pow(t6 + t7 + 1, 2);
    real_type result__ = -4 * xo__h / t9 * t5 / t2 * xo__retrograde * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_3_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_3_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t2   = 1.0 / t1;
    real_type t4   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type t7   = t5 + t6 + 1;
    real_type t14  = t7 * t7;
    real_type result__ = 2 / t7 * t4 * t2 * xo__retrograde - 4 / t14 * t4 * t2 * xo__retrograde * t6;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_3_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_3_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_3_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_3_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t2   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type result__ = 1.0 / (t7 + t8 + 1) / t4 / t2 * xo__k * xo__retrograde;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_3_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_3_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t4   = sqrt(xo__muS);
    real_type t5   = xo__h * xo__h;
    real_type t6   = xo__k * xo__k;
    real_type result__ = 2 / (t5 + t6 + 1) * t4 / t1 * xo__k;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_3_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = cos(xo__L);
    real_type t5   = sqrt(xo__p);
    real_type t6   = 1.0 / t5;
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type t9   = t7 + t8 + 1;
    real_type t15  = sin(xo__L);
    real_type t23  = t9 * t9;
    real_type result__ = 2 / t9 * t6 * (t2 + xo__f) * t1 - 4 * xo__h / t23 * t6 * (t15 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t2 * xo__h + xo__f * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_4_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = cos(xo__L);
    real_type t5   = sqrt(xo__p);
    real_type t6   = 1.0 / t5;
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type t9   = t7 + t8 + 1;
    real_type t10  = t9 * t9;
    real_type t12  = 1.0 / t10 * t6;
    real_type t17  = sin(xo__L);
    real_type t24  = (t17 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t2 * xo__h + xo__f * xo__h) * t1;
    real_type result__ = -8 * xo__h * t12 * (t2 + xo__f) * t1 + 16 * t7 / t10 / t9 * t6 * t24 - 4 * t12 * t24;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_4_4( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_4_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = cos(xo__L);
    real_type t5   = sqrt(xo__p);
    real_type t6   = 1.0 / t5;
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type t9   = t7 + t8 + 1;
    real_type t10  = t9 * t9;
    real_type t12  = 1.0 / t10 * t6;
    real_type t16  = sin(xo__L);
    real_type result__ = -4 * xo__k * t12 * (t2 + xo__f) * t1 - 4 * xo__h * t12 * (t16 * xo__retrograde + xo__g * xo__retrograde) * t1 + 16 * xo__k * xo__h / t10 / t9 * t6 * (t16 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t2 * xo__h + xo__f * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_4_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_4_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = sin(xo__L);
    real_type t4   = sqrt(xo__p);
    real_type t5   = 1.0 / t4;
    real_type t6   = xo__h * xo__h;
    real_type t7   = xo__k * xo__k;
    real_type t8   = t6 + t7 + 1;
    real_type t14  = cos(xo__L);
    real_type t19  = t8 * t8;
    real_type result__ = -2 / t8 * t5 * t2 * t1 - 4 * xo__h / t19 * t5 * (t14 * xo__k * xo__retrograde - t2 * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_4_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_4_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = 1.0 / t1;
    real_type t3   = cos(xo__L);
    real_type t6   = sqrt(xo__p);
    real_type t7   = 1.0 / t6;
    real_type t8   = xo__h * xo__h;
    real_type t9   = xo__k * xo__k;
    real_type t10  = t8 + t9 + 1;
    real_type t15  = sin(xo__L);
    real_type t23  = t10 * t10;
    real_type result__ = 1.0 / t10 * t7 * (t3 + xo__f) * t2 - 2 * xo__h / t23 * t7 * (t15 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t3 * xo__h + xo__f * xo__h) * t2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_4_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_4_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = sin(xo__L);
    real_type t7   = sqrt(xo__p);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t12  = pow(t9 + t10 + 1, 2);
    real_type result__ = -4 * xo__h / t12 / t7 * (t2 * xo__k + xo__g * xo__k) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_4_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = sin(xo__L);
    real_type t7   = sqrt(xo__p);
    real_type t8   = 1.0 / t7;
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t20  = cos(xo__L);
    real_type t25  = t11 * t11;
    real_type result__ = 2 / t11 * t8 * (t2 * xo__retrograde + xo__g * xo__retrograde) * t1 - 4 * xo__k / t25 * t8 * (t2 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t20 * xo__h + xo__f * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_5_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = sin(xo__L);
    real_type t7   = sqrt(xo__p);
    real_type t8   = 1.0 / t7;
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type t11  = t9 + t10 + 1;
    real_type t12  = t11 * t11;
    real_type t14  = 1.0 / t12 * t8;
    real_type t22  = cos(xo__L);
    real_type t26  = (t2 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t22 * xo__h + xo__f * xo__h) * t1;
    real_type result__ = -8 * xo__k * t14 * (t2 * xo__retrograde + xo__g * xo__retrograde) * t1 + 16 * t10 / t12 / t11 * t8 * t26 - 4 * t14 * t26;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_5_5( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_5_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t3   = cos(xo__L);
    real_type t4   = sqrt(xo__p);
    real_type t5   = 1.0 / t4;
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type t9   = t7 + t8 + 1;
    real_type t16  = sin(xo__L);
    real_type t20  = t9 * t9;
    real_type result__ = 2 / t9 * t5 * t3 * xo__retrograde * t1 - 4 * xo__k / t20 * t5 * (t3 * xo__k * xo__retrograde - t16 * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_5_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_5_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = 1.0 / t1;
    real_type t3   = sin(xo__L);
    real_type t8   = sqrt(xo__p);
    real_type t9   = 1.0 / t8;
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type t12  = t10 + t11 + 1;
    real_type t20  = cos(xo__L);
    real_type t25  = t12 * t12;
    real_type result__ = 1.0 / t12 * t9 * (t3 * xo__retrograde + xo__g * xo__retrograde) * t2 - 2 * xo__k / t25 * t9 * (t3 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t20 * xo__h + xo__f * xo__h) * t2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_5_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_5_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__p);
    real_type t2   = 1.0 / t1;
    real_type t3   = sqrt(xo__muS);
    real_type t5   = sin(xo__L);
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type t9   = t7 + t8 + 1;
    real_type t18  = t9 * t9;
    real_type result__ = 2 / t9 * (t5 + xo__g) * t3 * t2 - 4 * xo__k / t18 * t2 * (t5 * xo__k + xo__g * xo__k) * t3;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_5_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t3   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = sqrt(xo__p);
    real_type t11  = xo__h * xo__h;
    real_type t12  = xo__k * xo__k;
    real_type result__ = 2 / (t11 + t12 + 1) / t9 * (t3 * xo__k * xo__retrograde - t5 * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_6_6( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t3   = sin(xo__L);
    real_type t5   = cos(xo__L);
    real_type t9   = sqrt(xo__p);
    real_type t11  = xo__h * xo__h;
    real_type t12  = xo__k * xo__k;
    real_type result__ = 2 / (t11 + t12 + 1) / t9 * (-t3 * xo__k * xo__retrograde - t5 * xo__h) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_6_6( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_6_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t4   = cos(xo__L);
    real_type t6   = sin(xo__L);
    real_type t10  = sqrt(xo__p);
    real_type t12  = xo__h * xo__h;
    real_type t13  = xo__k * xo__k;
    real_type result__ = 1.0 / (t12 + t13 + 1) / t10 * (t4 * xo__k * xo__retrograde - t6 * xo__h) / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_6_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_6_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t3   = cos(xo__L);
    real_type t4   = sqrt(xo__p);
    real_type t7   = xo__h * xo__h;
    real_type t8   = xo__k * xo__k;
    real_type result__ = 2 / (t7 + t8 + 1) / t4 * t3 * xo__k * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_6_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t4   = sin(xo__L);
    real_type t8   = cos(xo__L);
    real_type t13  = sqrt(xo__p);
    real_type t15  = xo__h * xo__h;
    real_type t16  = xo__k * xo__k;
    real_type result__ = 1.0 / (t15 + t16 + 1) / t13 * (t4 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t8 * xo__h + xo__f * xo__h) / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_7_7( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t5   = sin(xo__L);
    real_type t9   = cos(xo__L);
    real_type t14  = sqrt(xo__p);
    real_type t16  = xo__h * xo__h;
    real_type t17  = xo__k * xo__k;
    real_type result__ = -1.0 / (t16 + t17 + 1) / t14 * (t5 * xo__k * xo__retrograde + xo__g * xo__k * xo__retrograde + t9 * xo__h + xo__f * xo__h) / t1 / xo__muS / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_7_7( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_7_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t3   = sin(xo__L);
    real_type t8   = sqrt(xo__p);
    real_type t10  = xo__h * xo__h;
    real_type t11  = xo__k * xo__k;
    real_type result__ = 1.0 / (t10 + t11 + 1) / t8 * (t3 * xo__k + xo__g * xo__k) / t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_7_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type t1   = sqrt(xo__muS);
    real_type t2   = sin(xo__L);
    real_type t7   = sqrt(xo__p);
    real_type t9   = xo__h * xo__h;
    real_type t10  = xo__k * xo__k;
    real_type result__ = 2 / (t9 + t10 + 1) / t7 * (t2 * xo__k + xo__g * xo__k) * t1;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_z_velocity__xo_D_8_8( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__h, real_type xo__k, real_type xo__L, real_type xo__muS, real_type xo__retrograde ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_z_velocity__xo_D_8_8( p={}, f={}, g={}, h={}, k={}, L={}, muS={}, retrograde={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__h, xo__k, xo__L, xo__muS, xo__retrograde, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type result__ = 1.0 / (t1 * xo__f + t3 * xo__g + 1) * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "ray( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type result__ = 1.0 / (t1 * xo__f + t3 * xo__g + 1);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_1( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_1_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type result__ = 0;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_1_1( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_1_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type t6   = pow(t1 * xo__f + t3 * xo__g + 1, 2);
    real_type result__ = -t1 / t6;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_1_2( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_1_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type t6   = pow(t1 * xo__f + t3 * xo__g + 1, 2);
    real_type result__ = -t3 / t6;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_1_3( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_1_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type t6   = pow(t1 * xo__f + t3 * xo__g + 1, 2);
    real_type result__ = -(t1 * xo__g - t3 * xo__f) / t6;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_1_4( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type t6   = pow(t1 * xo__f + t3 * xo__g + 1, 2);
    real_type result__ = -t1 / t6 * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_2( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_2_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type t5   = t1 * xo__f + t3 * xo__g + 1;
    real_type t6   = t5 * t5;
    real_type t10  = t1 * t1;
    real_type result__ = 2 * t10 / t6 / t5 * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_2_2( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_2_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type t5   = t1 * xo__f + t3 * xo__g + 1;
    real_type t6   = t5 * t5;
    real_type result__ = 2 * t3 * t1 / t6 / t5 * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_2_3( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_2_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type t5   = t1 * xo__f + t3 * xo__g + 1;
    real_type t6   = t5 * t5;
    real_type result__ = 2 * (t1 * xo__g - t3 * xo__f) * t1 / t6 / t5 * xo__p + t3 / t6 * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_2_4( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type t6   = pow(t1 * xo__f + t3 * xo__g + 1, 2);
    real_type result__ = -t3 / t6 * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_3( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_3_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type t5   = t1 * xo__f + t3 * xo__g + 1;
    real_type t6   = t5 * t5;
    real_type t10  = t3 * t3;
    real_type result__ = 2 * t10 / t6 / t5 * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_3_3( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_3_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type t5   = t1 * xo__f + t3 * xo__g + 1;
    real_type t6   = t5 * t5;
    real_type result__ = 2 * (t1 * xo__g - t3 * xo__f) * t3 / t6 / t5 * xo__p - t1 / t6 * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_3_4( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t3   = sin(xo__L);
    real_type t6   = pow(t1 * xo__f + t3 * xo__g + 1, 2);
    real_type result__ = -(t1 * xo__g - t3 * xo__f) / t6 * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_4( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_ray__xo_D_4_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L ) {
    real_type t1   = cos(xo__L);
    real_type t2   = t1 * xo__f;
    real_type t3   = sin(xo__L);
    real_type t4   = t3 * xo__g;
    real_type t5   = 1 + t2 + t4;
    real_type t6   = t5 * t5;
    real_type t13  = pow(t1 * xo__g - t3 * xo__f, 2);
    real_type result__ = 2 * t13 / t6 / t5 * xo__p - (-t2 - t4) / t6 * xo__p;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_ray__xo_D_4_4( p={}, f={}, g={}, L={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type result__ = sqrt(1.0 / xo__p * ((2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1) * xo__muS);
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "vel( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = (2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1;
    real_type t13  = sqrt(1.0 / xo__p * t9 * xo__muS);
    real_type t16  = xo__p * xo__p;
    real_type result__ = -1.0 / t16 * t9 * xo__muS / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_1( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_1_1( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = (2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1;
    real_type t12  = 1.0 / xo__p * t9 * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t16  = xo__muS * xo__muS;
    real_type t18  = t9 * t9;
    real_type t19  = xo__p * xo__p;
    real_type t20  = t19 * t19;
    real_type result__ = -1.0 / t20 * t18 * t16 / t13 / t12 / 4 + 1.0 / t19 / xo__p * t9 * xo__muS / t13;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_1_1( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_1_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = (2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1;
    real_type t12  = 1.0 / xo__p * t9 * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t16  = xo__muS * xo__muS;
    real_type t18  = xo__p * xo__p;
    real_type t23  = 2 * t1 + 2 * xo__f;
    real_type result__ = t23 / t18 / xo__p * t9 * t16 / t13 / t12 / 4 - 1.0 / t18 * t23 * xo__muS / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_1_2( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_1_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = (2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1;
    real_type t12  = 1.0 / xo__p * t9 * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t16  = xo__muS * xo__muS;
    real_type t18  = xo__p * xo__p;
    real_type t23  = 2 * t5 + 2 * xo__g;
    real_type result__ = t23 / t18 / xo__p * t9 * t16 / t13 / t12 / 4 - 1.0 / t18 * t23 * xo__muS / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_1_3( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_1_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = (2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1;
    real_type t12  = 1.0 / xo__p * t9 * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t16  = xo__muS * xo__muS;
    real_type t18  = xo__p * xo__p;
    real_type t25  = 2 * t1 * xo__g - 2 * t5 * xo__f;
    real_type result__ = t25 / t18 / xo__p * t9 * t16 / t13 / t12 / 4 - 1.0 / t18 * t25 * xo__muS / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_1_4( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_1_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = (2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1;
    real_type t12  = 1.0 / xo__p * t9 * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t17  = t9 * t9;
    real_type t18  = xo__p * xo__p;
    real_type result__ = 1.0 / t18 / xo__p * t17 * xo__muS / t13 / t12 / 4 - 1.0 / t18 * t9 / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_1_5( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t11  = 1.0 / xo__p;
    real_type t13  = sqrt(t11 * ((2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1) * xo__muS);
    real_type result__ = t11 * (2 * t1 + 2 * xo__f) * xo__muS / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_2( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_2_2( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t11  = 1.0 / xo__p;
    real_type t12  = t11 * ((2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1) * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t16  = xo__muS * xo__muS;
    real_type t20  = pow(2 * t1 + 2 * xo__f, 2);
    real_type t21  = xo__p * xo__p;
    real_type result__ = -1.0 / t21 * t20 * t16 / t13 / t12 / 4 + t11 * xo__muS / t13;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_2_2( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_2_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t12  = 1.0 / xo__p * ((2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1) * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t16  = xo__muS * xo__muS;
    real_type t20  = xo__p * xo__p;
    real_type result__ = -(2 * t5 + 2 * xo__g) / t20 * (2 * t1 + 2 * xo__f) * t16 / t13 / t12 / 4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_2_3( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_2_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t11  = 1.0 / xo__p;
    real_type t12  = t11 * ((2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1) * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t16  = xo__muS * xo__muS;
    real_type t20  = xo__p * xo__p;
    real_type result__ = -(2 * t1 * xo__g - 2 * t5 * xo__f) / t20 * (2 * t1 + 2 * xo__f) * t16 / t13 / t12 / 4 - t11 * t5 * xo__muS / t13;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_2_4( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_2_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = (2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1;
    real_type t11  = 1.0 / xo__p;
    real_type t12  = t11 * t9 * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t18  = 2 * t1 + 2 * xo__f;
    real_type t19  = xo__p * xo__p;
    real_type result__ = -t9 / t19 * t18 * xo__muS / t13 / t12 / 4 + t11 * t18 / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_2_5( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t11  = 1.0 / xo__p;
    real_type t13  = sqrt(t11 * ((2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1) * xo__muS);
    real_type result__ = t11 * (2 * t5 + 2 * xo__g) * xo__muS / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_3( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_3_3( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t11  = 1.0 / xo__p;
    real_type t12  = t11 * ((2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1) * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t16  = xo__muS * xo__muS;
    real_type t20  = pow(2 * t5 + 2 * xo__g, 2);
    real_type t21  = xo__p * xo__p;
    real_type result__ = -1.0 / t21 * t20 * t16 / t13 / t12 / 4 + t11 * xo__muS / t13;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_3_3( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_3_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t11  = 1.0 / xo__p;
    real_type t12  = t11 * ((2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1) * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t16  = xo__muS * xo__muS;
    real_type t20  = xo__p * xo__p;
    real_type result__ = -(2 * t1 * xo__g - 2 * t5 * xo__f) / t20 * (2 * t5 + 2 * xo__g) * t16 / t13 / t12 / 4 + t11 * t1 * xo__muS / t13;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_3_4( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_3_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = (2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1;
    real_type t11  = 1.0 / xo__p;
    real_type t12  = t11 * t9 * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t18  = 2 * t5 + 2 * xo__g;
    real_type t19  = xo__p * xo__p;
    real_type result__ = -t9 / t19 * t18 * xo__muS / t13 / t12 / 4 + t11 * t18 / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_3_5( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t11  = 1.0 / xo__p;
    real_type t13  = sqrt(t11 * ((2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1) * xo__muS);
    real_type result__ = t11 * (2 * t1 * xo__g - 2 * t5 * xo__f) * xo__muS / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_4( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_4_4( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t11  = 1.0 / xo__p;
    real_type t12  = t11 * ((2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1) * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t16  = xo__muS * xo__muS;
    real_type t22  = pow(2 * t1 * xo__g - 2 * t5 * xo__f, 2);
    real_type t23  = xo__p * xo__p;
    real_type result__ = -1.0 / t23 * t22 * t16 / t13 / t12 / 4 + t11 * (-2 * t1 * xo__f - 2 * t5 * xo__g) * xo__muS / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_4_4( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_4_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = (2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1;
    real_type t11  = 1.0 / xo__p;
    real_type t12  = t11 * t9 * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t20  = 2 * t1 * xo__g - 2 * t5 * xo__f;
    real_type t21  = xo__p * xo__p;
    real_type result__ = -t9 / t21 * t20 * xo__muS / t13 / t12 / 4 + t11 * t20 / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_4_5( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = (2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1;
    real_type t11  = 1.0 / xo__p;
    real_type t13  = sqrt(t11 * t9 * xo__muS);
    real_type result__ = t11 * t9 / t13 / 2;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_5( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

  real_type
  astro_vel__xo_D_5_5( real_type xo__p, real_type xo__f, real_type xo__g, real_type xo__L, real_type xo__muS ) {
    real_type t1   = cos(xo__L);
    real_type t5   = sin(xo__L);
    real_type t9   = (2 * t1 + xo__f) * xo__f + (2 * t5 + xo__g) * xo__g + 1;
    real_type t12  = 1.0 / xo__p * t9 * xo__muS;
    real_type t13  = sqrt(t12);
    real_type t16  = t9 * t9;
    real_type t18  = xo__p * xo__p;
    real_type result__ = -1.0 / t18 * t16 / t13 / t12 / 4;
    if ( m_debug ) {
      UTILS_ASSERT(
        Utils::is_finite(result__),
        "astro_vel__xo_D_5_5( p={}, f={}, g={}, L={}, muS={} ) return {}\n",
        xo__p, xo__f, xo__g, xo__L, xo__muS, result__
      );
    }
    return result__;
  }

}

// EOF: gtoc11_2burn_pars_T_fixed_Methods_UserFunctions.cc
