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
  ) {

    real_type Dr[3], Dt[3], Dn[3];
    point_and_velocity_to_Frenet_RTN( Z, Z+3, Dr, Dt, Dn );

    real_type x        = Z[0];
    real_type y        = Z[1];
    real_type z        = Z[2];
    real_type vx       = Z[3];
    real_type vy       = Z[4];
    real_type vz       = Z[5];
    real_type Thrust_x = Thrust_r*Dr[0]+Thrust_t*Dt[0]+Thrust_n*Dn[0];
    real_type Thrust_y = Thrust_r*Dr[1]+Thrust_t*Dt[1]+Thrust_n*Dn[1];
    real_type Thrust_z = Thrust_r*Dr[2]+Thrust_t*Dt[2]+Thrust_n*Dn[2];
    real_type r2       = x*x+y*y+z*z;
    real_type r        = sqrt(r2);
    real_type mur3     = (mu/r2)/r;

    res[0] = vx;
    res[1] = vy;
    res[2] = vz;
    res[3] = Thrust_x/m - mur3*x;
    res[4] = Thrust_y/m - mur3*y;
    res[5] = Thrust_z/m - mur3*z;
  }

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
  ) {

    real_type p = Z[0];
    real_type f = Z[1];
    real_type g = Z[2];
    real_type h = Z[3];
    real_type k = Z[4];
    real_type L = Z[5];

    real_type sinL = sin(L);
    real_type cosL = cos(L);

    real_type w  = 1+(f*cosL+g*sinL);
    real_type s2 = 1+(h*h+k*k);
    real_type bf = sqrt(p/mu);

    real_type C_r = bf*Thrust_r/m;
    real_type C_t = bf*Thrust_t/m/w;
    real_type C_n = bf*Thrust_n/m/w;

    res[0] = 0;
    res[1] = 0;
    res[2] = 0;
    res[3] = 0;
    res[4] = 0;
    res[5] = 0;

    res[1] += C_r*sinL;
    res[2] -= C_r*cosL;

    res[0] += C_t*2*p;
    res[1] += C_t*((w+1)*cosL+f);
    res[2] += C_t*((w+1)*sinL+g);

    res[1] -= C_n*(h*sinL-k*cosL)*g;
    res[2] += C_n*(h*sinL-k*cosL)*f;
    res[3] += C_n*s2*cosL/2;
    res[4] += C_n*s2*sinL/2;
    res[5] += C_n*(h*sinL-k*cosL);

    res[5] += sqrt(p*mu)*power2(w/p);
  }
}
