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
 |      Università degli Studi di Trento                                    |
 |      Via Mesiano 77, I-38050 Trento, Italy                               |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Astro.hh"

namespace AstroLib {

  using autodiff::detail::val;

  template <size_t N>
  DualN<N>
  Astro::x_position( DualN<N> const & t ) const {
    real_type L[4];
    eval_L( val(t), L, 0 );
    real_type p    = m_EQ.p;
    real_type f    = m_EQ.f;
    real_type g    = m_EQ.g;
    real_type h    = m_EQ.h;
    real_type k    = m_EQ.k;
    real_type cosL = cos(L[0]);
    real_type sinL = sin(L[0]);
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type hk   = h*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;
    real_type I    = m_EQ.retrograde ? -1 : 1;

    real_type x_pos = (1+h2-k2)*X+2*I*hk*Y;

    DualN<N> res(x_pos);

    if constexpr ( N > 0 ) {
      real_type bf1   = sqrt(m_muS/p)/(1+h2+k2);
      real_type cosLf = bf1 * (cosL+f);
      real_type sinLg = bf1 * (sinL+g);
      real_type x_vel = I*2*hk*cosLf - (1+h2-k2)*sinLg;

      if constexpr ( N > 1 ) {
        real_type cosL_D = -bf1*sinL*L[1];
        real_type sinL_D =  bf1*cosL*L[1];
        real_type x_acc  = I*2*hk*cosL_D - (1+h2-k2)*sinL_D;
        real_type t_grad = val(t.grad);

        res.grad      = x_vel * t_grad;
        res.grad.grad = x_acc*(t_grad*t_grad)+x_vel*t.grad.grad;
      } else {
        res.grad = x_vel * t.grad;
      }
    }
    return res;
  }

  template <size_t N>
  DualN<N>
  Astro::y_position( DualN<N> const & t ) const {
    real_type L[4];
    eval_L( val(t), L, 0 );
    real_type cosL = cos(L[0]);
    real_type sinL = sin(L[0]);
    real_type p    = m_EQ.p;
    real_type f    = m_EQ.f;
    real_type g    = m_EQ.g;
    real_type h    = m_EQ.h;
    real_type k    = m_EQ.k;
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type hk   = h*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;
    real_type I    = m_EQ.retrograde ? -1 : 1;

    real_type y_pos = I*(1-h2+k2)*Y+2*hk*X;

    DualN<N> res(y_pos);

    if constexpr ( N > 0 ) {
      real_type bf1   = sqrt(m_muS/p)/(1+h2+k2);
      real_type cosLf = bf1 * (cosL+f);
      real_type sinLg = bf1 * (sinL+g);
      real_type y_vel = I*(1-h2+k2)*cosLf - 2*hk*sinLg;

      if constexpr ( N > 1 ) {
        real_type cosL_D = -bf1*sinL*L[1];
        real_type sinL_D =  bf1*cosL*L[1];
        real_type y_acc  = I*(1-h2+k2)*cosL_D - 2*hk*sinL_D;
        real_type t_grad = val(t.grad);

        res.grad      = y_vel * t_grad;
        res.grad.grad = y_acc*(t_grad*t_grad)+y_vel*t.grad.grad;

      } else {
        res.grad = y_vel * t.grad;
      }
    }

    return res;

  }
  
  template <size_t N>
  DualN<N> 
  Astro::z_position( DualN<N> const & t ) const {
    real_type L[4];
    eval_L( val(t), L, 0 );
    real_type cosL = cos(L[0]);
    real_type sinL = sin(L[0]);
    real_type p    = m_EQ.p;
    real_type f    = m_EQ.f;
    real_type g    = m_EQ.g;
    real_type h    = m_EQ.h;
    real_type k    = m_EQ.k;
    real_type h2   = h*h;
    real_type k2   = k*k;
    real_type bf   = p/((1+f*cosL+g*sinL)*(1+h2+k2));
    real_type X    = bf*cosL;
    real_type Y    = bf*sinL;
    real_type I    = m_EQ.retrograde ? -1 : 1;

    real_type z_pos = 2*(h*Y-I*k*X);

    DualN<N> res(z_pos);

    if constexpr ( N > 0 ) {

      real_type bf1   = sqrt(m_muS/p)/(1+h2+k2);
      real_type cosLf = bf1 * (cosL+f);
      real_type sinLg = bf1 * (sinL+g);
      real_type z_vel = 2 * ( h*cosLf + I*k*sinLg );

      if constexpr ( N > 1 ) {

        real_type cosL_D = -bf1*sinL*L[1];
        real_type sinL_D =  bf1*cosL*L[1];
        real_type z_acc  = 2 * ( h*cosL_D + I*k*sinL_D );
        real_type t_grad = val(t.grad);

        res.grad      = z_vel * t_grad;
        res.grad.grad = z_acc*(t_grad*t_grad)+z_vel*t.grad.grad;

      } else {
        res.grad = z_vel * t.grad;
      }
    }

    return res;

  }

  template <size_t N>
  DualN<N>
  Astro::x_velocity( DualN<N> const & t ) const {
    real_type L[4];
    eval_L( val(t), L, 0 );
    real_type cosL  = cos(L[0]);
    real_type sinL  = sin(L[0]);
    real_type p     = m_EQ.p;
    real_type f     = m_EQ.f;
    real_type g     = m_EQ.g;
    real_type h     = m_EQ.h;
    real_type k     = m_EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type hk    = h*k;
    real_type bf    = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cosL+f);
    real_type sinLg = bf * (sinL+g);
    real_type I     = m_EQ.retrograde ? -1 : 1;

    real_type x_vel = I*2*hk*cosLf - (1+h2-k2)*sinLg;

    DualN<N> res(x_vel);

    if constexpr ( N > 0 ) {

      real_type cosL_D = -bf*sinL*L[1];
      real_type sinL_D =  bf*cosL*L[1];
      real_type x_acc  = I*2*hk*cosL_D - (1+h2-k2)*sinL_D;

      if constexpr ( N > 1 ) {

        real_type cosL_2 = bf * ( -sinL*L[2] - cosL * power2(L[1]) );
        real_type sinL_2 = bf * (  cosL*L[2] - sinL * power2(L[1]) );
        real_type x_jerk = I*2*hk*cosL_2 - (1+h2-k2)*sinL_2;
        real_type t_grad = val(t.grad);

        res.grad      = x_acc * t_grad;
        res.grad.grad = x_jerk*(t_grad*t_grad)+x_acc*t.grad.grad;

      } else {
        res.grad = x_acc * t.grad;
      }
    }

    return res;
  }
  
  template <size_t N>
  DualN<N>
  Astro::y_velocity( DualN<N> const & t ) const {
    real_type L[4];
    eval_L( val(t), L, 0 );
    real_type cosL  = cos(L[0]);
    real_type sinL  = sin(L[0]);
    real_type p     = m_EQ.p;
    real_type f     = m_EQ.f;
    real_type g     = m_EQ.g;
    real_type h     = m_EQ.h;
    real_type k     = m_EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type hk    = h*k;
    real_type bf    = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cosL+f);
    real_type sinLg = bf * (sinL+g);
    real_type I     = m_EQ.retrograde ? -1 : 1;

    real_type y_vel = I*(1-h2+k2)*cosLf - 2*hk*sinLg;

    DualN<N> res(y_vel);

    if constexpr ( N > 0 ) {

      real_type bf     = sqrt(m_muS/p)/(1+h2+k2);
      real_type cosL_D = -bf*sinL*L[1];
      real_type sinL_D =  bf*cosL*L[1];
      real_type y_acc  = I*(1-h2+k2)*cosL_D - 2*hk*sinL_D;

      if constexpr ( N > 1 ) {

        real_type cosL_2 = bf * ( -sinL*L[2] - cosL * power2(L[1]) );
        real_type sinL_2 = bf * (  cosL*L[2] - sinL * power2(L[1]) );
        real_type y_jerk = I*(1-h2+k2)*cosL_2 - 2*hk*sinL_2;
        real_type t_grad = val(t.grad);

        res.grad      = y_acc * t_grad;
        res.grad.grad = y_jerk*(t_grad*t_grad)+y_acc*t.grad.grad;
      } else {
        res.grad = y_acc * t.grad;
      }
    }

    return res;
  }
  
  template <size_t N>
  DualN<N>
  Astro::z_velocity( DualN<N> const & t ) const {
    real_type L[4];
    eval_L( val(t), L, 0 );
    real_type cosL  = cos(L[0]);
    real_type sinL  = sin(L[0]);
    real_type p     = m_EQ.p;
    real_type f     = m_EQ.f;
    real_type g     = m_EQ.g;
    real_type h     = m_EQ.h;
    real_type k     = m_EQ.k;
    real_type h2    = h*h;
    real_type k2    = k*k;
    real_type bf    = sqrt(m_muS/p)/(1+h2+k2);
    real_type cosLf = bf * (cosL+f);
    real_type sinLg = bf * (sinL+g);
    real_type I     = m_EQ.retrograde ? -1 : 1;

    real_type z_vel = 2 * ( h*cosLf + I*k*sinLg );

    DualN<N> res(z_vel);

    if constexpr ( N > 0 ) {

      real_type cosL_D = -bf*sinL*L[1];
      real_type sinL_D =  bf*cosL*L[1];
      real_type z_acc  = 2 * ( h*cosL_D + I*k*sinL_D );

      if constexpr ( N > 1 ) {

        real_type cosL_2 = bf * ( -sinL*L[2] - cosL * power2(L[1]) );
        real_type sinL_2 = bf * (  cosL*L[2] - sinL * power2(L[1]) );
        real_type z_jerk = 2 * ( h*cosL_2 + I*k*sinL_2 );
        real_type t_grad = val(t.grad);

        res.grad      = z_acc * t_grad;
        res.grad.grad = z_jerk*(t_grad*t_grad)+z_acc*t.grad.grad;

      } else {
        res.grad = z_acc * t.grad;
      }
    }

    return res;
  }

  template <size_t N>
  DualN<N>
  Astro::radius_by_L( DualN<N> const & L ) const {
    DualN<N> res{ this->radius_by_L( val(L) ) };
    if constexpr ( N > 0 ) {
      real_type r_D{ radius_by_L_D( val(L) ) };
      if constexpr ( N > 1 ) {
        real_type r_DD{ radius_by_L_DD( val(L) ) };
        real_type L_grad = val(L.grad);

        res.grad      = r_D * L_grad;
        res.grad.grad = r_DD*(L_grad*L_grad)+r_D*L.grad.grad;
      } else {
        res.grad = r_D * L.grad;
      }
    }
    return res;
  }

  template DualN<0> Astro::x_position<0>( DualN<0> const & t ) const;
  template DualN<1> Astro::x_position<1>( DualN<1> const & t ) const;
  template DualN<2> Astro::x_position<2>( DualN<2> const & t ) const;

  template DualN<0> Astro::y_position<0>( DualN<0> const & t ) const;
  template DualN<1> Astro::y_position<1>( DualN<1> const & t ) const;
  template DualN<2> Astro::y_position<2>( DualN<2> const & t ) const;

  template DualN<0> Astro::z_position<0>( DualN<0> const & t ) const;
  template DualN<1> Astro::z_position<1>( DualN<1> const & t ) const;
  template DualN<2> Astro::z_position<2>( DualN<2> const & t ) const;

  template DualN<0> Astro::x_velocity<0>( DualN<0> const & t ) const;
  template DualN<1> Astro::x_velocity<1>( DualN<1> const & t ) const;
  template DualN<2> Astro::x_velocity<2>( DualN<2> const & t ) const;

  template DualN<0> Astro::y_velocity<0>( DualN<0> const & t ) const;
  template DualN<1> Astro::y_velocity<1>( DualN<1> const & t ) const;
  template DualN<2> Astro::y_velocity<2>( DualN<2> const & t ) const;

  template DualN<0> Astro::z_velocity<0>( DualN<0> const & t ) const;
  template DualN<1> Astro::z_velocity<1>( DualN<1> const & t ) const;
  template DualN<2> Astro::z_velocity<2>( DualN<2> const & t ) const;

  template DualN<0> Astro::radius_by_L<0>( DualN<0> const & L ) const;
  template DualN<1> Astro::radius_by_L<1>( DualN<1> const & L ) const;
  template DualN<2> Astro::radius_by_L<2>( DualN<2> const & L ) const;

}
