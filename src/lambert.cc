#include "astro.hh"

#include <cmath>
#include <limits>
#include <algorithm> // for min, max

namespace astro {

  using namespace std;

  static real_type eps = std::numeric_limits<real_type>::epsilon() ;

  template <typename T>
  int
  sign(T val) {
    if      ( val < 0 ) return -1 ;
    else if ( val > 0 ) return  1 ;
    return 0 ;
  }

  // series approximation to T(x) and its derivatives (used for near-parabolic cases)
  static
  void
  sigmax( real_type const y,
          real_type     & sig,
          real_type     & dsigdx,
          real_type     & d2sigdx2,
          real_type     & d3sigdx3 ) {

    // preload the factors [an]
    // (25 factors is more than enough for 16-digit accuracy)
    static real_type an[] = {
      4./3.,
      4.000000000000000e-001, 2.142857142857143e-001, 4.629629629629630e-002,
      6.628787878787879e-003, 7.211538461538461e-004, 6.365740740740740e-005,
      4.741479925303455e-006, 3.059406328320802e-007, 1.742836409255060e-008,
      8.892477331109578e-010, 4.110111531986532e-011, 1.736709384841458e-012,
      6.759767240041426e-014, 2.439123386614026e-015, 8.203411614538007e-017,
      2.583771576869575e-018, 7.652331327976716e-020, 2.138860629743989e-021,
      5.659959451165552e-023, 1.422104833817366e-024, 3.401398483272306e-026,
      7.762544304774155e-028, 1.693916882090479e-029, 3.541295006766860e-031,
      7.105336187804402e-033
    };
      
    // powers of y
    // sigma
    sig = an[25] ;
    for ( int i = 24 ; i >= 0 ; --i ) sig = sig*y+an[i] ;

    // dsigma / dx
    dsigdx = 25*an[25] ;
    for ( int i = 24 ; i > 0 ; --i ) dsigdx = dsigdx*y+i*an[i] ;

    // d2sigma / dx2
    d2sigdx2 = (25*24)*an[25] ;
    for ( int i = 24 ; i > 1 ; --i ) d2sigdx2 = d2sigdx2*y+(i*(i-1))*an[i] ;

    // d3sigma / dx3
    d3sigdx3 = (25*24*23)*an[25] ;
    for ( int i = 24 ; i > 2 ; --i ) d3sigdx3 = d3sigdx3*y+(i*(i-1)*(i-2))*an[i] ;

  }

  /*\
   | Lancaster & Blanchard's function, and three derivatives thereof
  \*/
  static
  void
  LancasterBlanchard( real_type const x_in,
                      real_type const q,
                      int       const m,
                      real_type       T[4],
                      int             nd ) {
    
    // protection against idiotic input
    real_type x = x_in ;
    if ( x < -1 ) { //  impossible; negative eccentricity
      x = abs(x) - 2 ;
    } else if ( x == -1 ) { // impossible; offset x slightly
      x += eps;
    }

    // compute parameter E
    real_type E  = x*x - 1 ;
    real_type q2 = q*q ;
    real_type q3 = q2*q ;
    real_type x2 = x*x ;

    // T(x), T'(x), T''(x)
    if ( x == 1 ) { // exactly parabolic; solutions known exactly
      // T(x)
      real_type q3 = q2*q ;
      T[0] = (4./3.)*(1-q3) ;
      // T'(x)
      if ( nd < 1 ) return ;
      real_type q5 = q2*q3 ;
      T[1] = (4./5.)*(q5 - 1) ;
      // T''(x)
      if ( nd < 2 ) return ;
      real_type q7 = q2*q5 ;
      T[2] = T[1] + (120./70.)*(1 - q7) ;
      // T'''(x)
      if ( nd < 3 ) return ;
      real_type q9 = q2*q7 ;
      T[3] = 3*(T[2] - T[1]) + (2400./1080.)*(q9 - 1) ;
    } else if ( abs(x-1) < 1e-2 ) { // near-parabolic; compute with series
      // evaluate sigma
      real_type sig1, dsigdx1, d2sigdx21, d3sigdx31 ;
      real_type sig2, dsigdx2, d2sigdx22, d3sigdx32 ;
      sigmax(-E,sig1, dsigdx1, d2sigdx21, d3sigdx31) ;
      sigmax(-E*q2,sig2, dsigdx2, d2sigdx22, d3sigdx32) ;
      // T(x)
      T[0] = sig1 - q3*sig2;
      // T'(x)
      if ( nd < 1 ) return ;
      real_type q5 = q2*q3 ;
      T[1] = 2*x*(q5*dsigdx2 - dsigdx1);
      // T''(x)
      if ( nd < 2 ) return ;
      real_type q7 = q2*q5 ;
      T[2] = T[1]/x + 4*x2*(d2sigdx21 - q7*d2sigdx22);
      // T'''(x)
      if ( nd < 3 ) return ;
      real_type q9 = q2*q7 ;
      T[3] = 3*(T[2]-T[1]/x)/x + 8*x2*(q9*d3sigdx32 - d3sigdx31);
    } else { // all other cases
      // compute all substitution functions
      real_type y  = sqrt(abs(E));
      real_type z  = sqrt(1 + q2*E);
      real_type f  = y*(z - q*x);
      real_type g  = x*z - q*E;
      real_type d  = 0 ;
      if      ( E < 0 ) d = atan2(f, g) + m_pi*m ;
      else if ( E > 0 ) d = log( pos(f+g) ) ;
      // %%%%%%d  = (E < 0)*(atan2(f, g) + pi*m) + (E > 0)*log( max(0, f + g) );
      // T(x)
      T[0] = 2*(x - q*z - d/y)/E;
      // T'(x)
      if ( nd < 1 ) return ;
      real_type q3 = q2*q ;
      T[1] = (4 - 4*q3*x/z - 3*x*T[0])/E;
      // T''(x)
      if ( nd < 2 ) return ;
      real_type z2 = z*z ;
      T[2] = (-4*q3/z * (1 - q2*x2/z2) - 3*T[0] - 3*x*T[1])/E;
      // T'''(x)
      if ( nd < 3 ) return ;
      T[3] = (4*q3/z2*((1 - q2*x2/z2) + 2*q2*x/z2*(z - x)) - 8*T[1] - 7*x*T[2])/E;
    }
  }

  /*\
   | -----------------------------------------------------------------
   | Lancaster & Blanchard version, with improvements by Gooding
   | Very reliable, moderately fast for both simple and complicated cases
   | -----------------------------------------------------------------
  \*/

  /*\
   | LAMBERT_LANCASTERBLANCHARD       High-Thrust Lambert-targeter
   |
   | lambert_LancasterBlanchard() uses the method developed by
   | Lancaster & Blancard, as described in their 1969 paper. Initial values,
   | and several details of the procedure, are provided by R.H. Gooding,
   | as described in his 1990 paper.
  \*/

  int_type
  lambert_LancasterBlanchard( real_type const r1vec[3],
                              real_type const r2vec[3],
                              real_type const tf_in,
                              int       const m_in,
                              real_type const muC,
                              real_type       V1[3],
                              real_type       V2[3] ) {

    // manipulate input
    real_type tol       = 1e-12 ;               // optimum for numerical noise v.s. actual precision
    real_type r1        = norm3(r1vec) ;        // magnitude of r1vec
    real_type r2        = norm3(r2vec) ;        // magnitude of r2vec
    real_type r1unit[3] = { r1vec[0]/r1, r1vec[1]/r1, r1vec[2]/r1 } ; // unit vector of r1vec
    real_type r2unit[3] = { r2vec[0]/r2, r2vec[1]/r2, r2vec[2]/r2 } ; // unit vector of r2vec
    real_type crsprod[3] ;
    cross( r1vec, r2vec, crsprod ) ;             // cross product of r1vec and r2vec
    real_type mcrsprd = norm3(crsprod) ;      // magnitude of that cross product
    real_type tmp[3] = { crsprod[0]/mcrsprd,
                         crsprod[1]/mcrsprd,
                         crsprod[2]/mcrsprd } ;
    real_type th1unit[3], th2unit[3] ;
    cross(tmp, r1unit, th1unit) ; // unit vectors in the tangential-directions
    cross(tmp, r2unit, th2unit) ;
    // make 100.4% sure it's in (-1 <= x <= +1)
    real_type dth = dot3(r1vec,r2vec)/r1/r2 ; // turn angle
    if      ( dth <= -1 ) dth = m_pi ;
    else if ( dth >=  1 ) dth = 0 ;
    else                  dth = acos( dth ) ;

    // if the long way was selected, the turn-angle must be negative
    // to take care of the direction of final velocity
    int_type  longway = sign(tf_in) ;
    real_type tf      = abs(tf_in) ;
    if (longway < 0) dth = dth-2*m_pi ;

    // left-branch
    int_type leftbranch = sign(m_in) ;
    int_type m          = m_in > 0 ? m_in : -m_in ;
      
    // define constants
    real_type c = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cos(dth));
    real_type s = (r1 + r2 + c) / 2;
    real_type T = sqrt(8*muC/(s*s*s)) * tf ;
    real_type q = sqrt(r1*r2)/s * cos(dth/2);

    // general formulae for the initial values (Gooding)
    // -------------------------------------------------

    // some initial values
    real_type T0 ;
    LancasterBlanchard(0, q, m, &T0, 0 );
    real_type Td  = T0 - T;
    real_type phr = 2*atan2(1 - q*q, 2*q) ;
    phr -= floor(phr/(2*m_pi))*2*m_pi ;

    // initial output is pessimistic
    //V1 = NaN(1,3);
    //V2 = V1;
    //extremal_distances = [NaN, NaN];
      
    // single-revolution case
    real_type x0 ;
    if (m == 0) {
      real_type x01, x02, x03, W, lambda ;
      x01 = T0*Td/4/T;
      if ( Td > 0 ) x0 = x01;
      else {
        x01 = Td/(4 - Td);
        x02 = -sqrt( -Td/(T+T0/2) );
        W   = x01 + 1.7*sqrt(2 - phr/m_pi);
        if (W >= 0) x03 = x01;
        else        x03 = x01 + pow(-W,1./16.)*(x02 - x01);
        lambda = 1 + x03*(1 + x01)/2 - 0.03*power2(x03)*sqrt(1 + x01);
        x0     = lambda*x03;
      }
      // this estimate might not give a solution
      if (x0 < -1) return -1 ;

    } else {
      // multi-revolution case
      // determine minimum Tp(x)
      real_type xMpi = 4/(3*m_pi*(2*m + 1));
      real_type xM0  = 0 ; // EMLMEX requires this one
      if      ( phr < m_pi ) xM0 = xMpi*pow(phr/m_pi,1./8.);
      else if ( phr > m_pi ) xM0 = xMpi*(2 - pow(2 - phr/m_pi,1./8.));

      // use Halley's method
      real_type   xM   = xM0 ;
      real_type   TT[] = { 1e100, 1e100, 1e100, 1e100 } ;
      real_type & Tp   = TT[1] ;
      real_type & Tpp  = TT[2] ;
      real_type & Tppp = TT[3] ;
      for ( int_type iterations = 0 ; abs(TT[1]) > tol ; ) {
        // iterations
        ++iterations ;
        // compute first three derivatives
        LancasterBlanchard(xM, q, m, TT, 3 );
        // new value of xM
        real_type xMp = xM;
        xM -= 2*Tp*Tpp/ (2*power2(Tpp) - Tp*Tppp);
        // escape clause
        if ( (iterations%7) == 0 ) xM = (xMp+xM)/2 ;
        // the method might fail. Exit in that case
        if ( iterations > 25 ) return -2 ;
      }
    
      // xM should be elliptic (-1 < x < 1)
      // (this should be impossible to go wrong)
      if ( (xM < -1) || (xM > 1) ) return -1 ;

      // corresponding time
      real_type TM ;
      LancasterBlanchard( xM, q, m, &TM, 0 );
          
      // T should lie above the minimum T
      if ( TM > T ) return -1 ;
          
      // find two initial values for second solution (again with lambda-type patch)
      // --------------------------------------------------------------------------
          
      // some initial values
      real_type TmTM  = T - TM ;
      real_type T0mTM = T0 - TM;
      LancasterBlanchard( xM, q, m, TT, 2) ;

      // first estimate (only if m > 0)
      if ( leftbranch > 0 ) {
        real_type x  = sqrt( TmTM / (Tpp/2 + TmTM/power2(1-xM)) );
        real_type W  = xM + x;
        W  = 4*W/(4 + TmTM) + power2(1 - W);
        x0 = x*(1 - (1 + m + (dth - 1/2)) / (1 + 0.15*m)*x*(W/2 + 0.03*x*sqrt(W))) + xM;
              
        // first estimate might not be able to yield possible solution
        if ( x0 > 1 ) return -1 ;
      } else {
        // second estimate (only if m > 0)
        if ( Td > 0 )
          x0 = xM - sqrt(TM/(Tpp/2 - TmTM*(Tpp/2/T0mTM - 1/power2(xM))));
        else {
          real_type x00 = Td / (4 - Td);
          real_type W   = x00 + 1.7*sqrt(2 - phr/m_pi);
          real_type x03 = x00 ;
          if ( W < 0 ) x03 -= sqrt(pow(-W,1./8.))*(x00 + sqrt(-Td/(1.5*T0 - Td))) ;
          W  = 4/(4 - Td);
          x0 = x03*(1 + (1 + m + 0.24*(dth - 1/2)) / (1 + 0.15*m)*x03*(W/2 - 0.03*x03*sqrt(W)));
        }
        // estimate might not give solutions
        if ( x0 < -1 ) return -1 ;
      }
    }

    // find root of Lancaster & Blancard's function
    // --------------------------------------------
    // (Halley's method)
    real_type x = x0 ;
    {
      real_type   TT[] = { 1./0., 1e100, 1e100 } ;
      real_type & Tx   = TT[0] ;
      real_type & Tp   = TT[1] ;
      real_type & Tpp  = TT[2] ;
      for ( int iterations = 0 ; abs(Tx) > tol ; ) {
        // iterations
        ++iterations ;
        // compute function value, and first two derivatives
        LancasterBlanchard(x, q, m, TT, 2) ;
        // find the root of the *difference* between the function value [T_x] and the required time [T]
        Tx -= T ;
        // new value of x
        real_type xp = x;
        x -= 2*Tx*Tp / (2*power2(Tp) - Tx*Tpp);
        // escape clause
        if ( (iterations%7) == 0 ) x = (xp+x)/2 ;
        // Halley's method might fail
        if ( iterations > 25 || isnan(Tx) ) return -2 ;
      }
    }
      
    // calculate terminal velocities
    // -----------------------------
      
    // constants required for this calculation
    real_type sigma, rho, z, gamma = sqrt(muC*s/2);
    if ( c == 0 ) {
      sigma = 1 ;
      rho   = 0 ;
      z     = abs(x);
    } else {
      sigma = 2*sqrt(r1*r2/(c*c)) * sin(dth/2);
      rho   = (r1 - r2)/c;
      z     = sqrt(1 + q*q*(x*x - 1));
    }
      
    // radial component
    real_type Vr1       = +gamma*((q*z - x) - rho*(q*z + x)) / r1;
    real_type Vr1vec[3] = { Vr1*r1unit[0], Vr1*r1unit[1], Vr1*r1unit[2] } ;
    real_type Vr2       = -gamma*((q*z - x) + rho*(q*z + x)) / r2;
    real_type Vr2vec[3] = { Vr2*r2unit[0], Vr2*r2unit[1], Vr2*r2unit[2] } ;
      
    // tangential component
    real_type Vtan1       = sigma * gamma * (z + q*x) / r1;
    real_type Vtan2       = sigma * gamma * (z + q*x) / r2;
    real_type Vtan1vec[3] = { Vtan1 * th1unit[0], Vtan1 * th1unit[1], Vtan1 * th1unit[2] } ;
    real_type Vtan2vec[3] = { Vtan2 * th2unit[0], Vtan2 * th2unit[1], Vtan2 * th2unit[2] } ;

    // Cartesian velocity
    V1[0] = Vtan1vec[0] + Vr1vec[0];
    V1[1] = Vtan1vec[1] + Vr1vec[1];
    V1[2] = Vtan1vec[2] + Vr1vec[2];
    V2[0] = Vtan2vec[0] + Vr2vec[0];
    V2[1] = Vtan2vec[1] + Vr2vec[1];
    V2[2] = Vtan2vec[2] + Vr2vec[2];

    return 1 ; // (success)
  }

  /*\
   % -----------------------------------------------------------------
   % Izzo's version:
   % Very fast, but not very robust for more complicated cases
   % -----------------------------------------------------------------
  \*/
  
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  int
  lambert( real_type const R1[3],
           real_type const R2[3],
           real_type const tf_in,
           int       const m_in,
           real_type const muC,
           real_type       V1[3],
           real_type       V2[3] ) {

  /* original documentation:
  %{
   This routine implements a new algorithm that solves Lambert's problem. The
   algorithm has two major characteristics that makes it favorable to other
   existing ones.

   1) It describes the generic orbit solution of the boundary condition
   problem through the variable X=log(1+cos(alpha/2)). By doing so the
   graph of the time of flight become defined in the entire real axis and
   resembles a straight line. Convergence is granted within few iterations
   for all the possible geometries (except, of course, when the transfer
   angle is zero). When multiple revolutions are considered the variable is
   X=tan(cos(alpha/2)*pi/2).

   2) Once the orbit has been determined in the plane, this routine
   evaluates the velocity vectors at the two points in a way that is not
   singular for the transfer angle approaching to pi (Lagrange coefficient
   based methods are numerically not well suited for this purpose).

   As a result Lambert's problem is solved (with multiple revolutions
   being accounted for) with the same computational effort for all
   possible geometries. The case of near 180 transfers is also solved
   efficiently.

   We note here that even when the transfer angle is exactly equal to pi
   the algorithm does solve the problem in the plane (it finds X), but it
   is not able to evaluate the plane in which the orbit lies. A solution
   to this would be to provide the direction of the plane containing the
   transfer orbit from outside. This has not been implemented in this
   routine since such a direction would depend on which application the
   transfer is going to be used in.

   please report bugs to dario.izzo@esa.int    
  %}

  % adjusted documentation:
  %{
   By default, the short-way solution is computed. The long way solution
   may be requested by giving a negative value to the corresponding
   time-of-flight [tf].

   For problems with |m| > 0, there are generally two solutions. By
   default, the right branch solution will be returned. The left branch
   may be requested by giving a negative value to the corresponding
   number of complete revolutions [m].
  %}

  % Authors
  % .-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.
  % Name       : Dr. Dario Izzo
  % E-mail     : dario.izzo@esa.int
  % Affiliation: ESA / Advanced Concepts Team (ACT)

  % Made more readible and optimized for speed by Rody P.S. Oldenhuis
  % Code available in MGA.M on   http://www.esa.int/gsp/ACT/inf/op/globopt.htm

  % last edited 12/Dec/2009

  % ADJUSTED FOR EML-COMPILATION 24/Dec/2009
  */
    // initial values
    real_type tol  = 1e-14 ;
    bool      bad  = false ;

    // work with non-dimensional units
    real_type r1       = norm3(R1) ;
    real_type V        = sqrt(muC/r1) ;
    real_type r1vec[3] = { R1[0]/r1, R1[1]/r1, R1[2]/r1 } ;
    real_type r2vec[3] = { R2[0]/r1, R2[1]/r1, R2[2]/r1 } ;
    real_type T        = r1/V ;
    real_type tf       = tf_in/T ; // also transform to seconds

    // relevant geometry parameters (non dimensional)
    real_type mr2vec = norm3(r2vec) ;
    // make 100% sure it's in (-1 <= dth <= +1)
    real_type dth = dot3(r1vec,r2vec)/mr2vec ;
    if      ( dth <= -1 ) dth = m_pi ;
    else if ( dth >=  1 ) dth = 0 ;
    else                  dth = acos( dth ) ;

    // decide whether to use the left or right branch
    // (for multi-revolution problems), and the long- or short way
    int_type m          = m_in ;
    int_type leftbranch = sign(m) ;
    int_type longway    = sign(tf) ;
    m  = m > 0 ? m : -m ;
    tf = abs(tf) ;
    if ( longway < 0 ) dth = 2*m_pi - dth ;

    // derived quantities
    //real_type c      = sqrt(1 + mr2vec^2 - 2*mr2vec*cos(dth)) ; % non-dimensional chord
    real_type tmp[3]   = { r2vec[0]-2*r1vec[0], r2vec[1]-2*r1vec[1], r2vec[2]-2*r1vec[2] } ;
    real_type c        = sqrt(1 + dot3(tmp,r2vec)) ;          // non-dimensional chord
    real_type s        = (1 + mr2vec + c)/2 ;                 // non-dimensional semi-perimeter
    real_type a_min    = s/2 ;                                // minimum energy ellipse semi major axis
    real_type Lambda   = sqrt(mr2vec)*cos(dth/2)/s ;          // lambda parameter (from BATTIN's book)
    real_type crossprd[3], nrmunit[3] ;
    cross(r1vec,r2vec,crossprd) ;                             // non-dimensional normal vectors
    real_type mcr      = norm3(crossprd) ;                    // magnitues thereof
    nrmunit[0] = crossprd[0]/mcr ;                            // unit vector thereof
    nrmunit[1] = crossprd[1]/mcr ;
    nrmunit[2] = crossprd[2]/mcr ;

    // Initial values
    // ---------------------------------------------------------

    // ELMEX requires this variable to be declared OUTSIDE the IF-statement
    real_type logt = log(tf) ; // avoid re-computing the same value
      

    real_type inn1 ; // first initial guess
    real_type inn2 ; // second initial guess
    real_type x1   ; // transformed first initial guess
    real_type x2   ; // transformed first second guess

    // single revolution (1 solution)
    if ( m == 0 ) {
      // initial values
      inn1 = -0.5233 ;       // first initial guess
      inn2 = +0.5233 ;       // second initial guess
      x1   = log(1 + inn1) ; // transformed first initial guess
      x2   = log(1 + inn2) ; // transformed first second guess
      // multiple revolutions (0, 1 or 2 solutions)
      // the returned soltuion depends on the sign of [m]
    } else {
      // select initial values
      if ( leftbranch < 0 ) {
        inn1 = -0.5234 ; // first initial guess, left branch
        inn2 = -0.2234 ; // second initial guess, left branch
      } else {
        inn1 = +0.7234 ; // first initial guess, right branch
        inn2 = +0.5234 ; // second initial guess, right branch
      }
      x1 = tan(inn1*m_pi_2) ; // transformed first initial guess
      x2 = tan(inn2*m_pi_2) ; // transformed first second guess
    }

    // since (inn1, inn2) < 0, initial estimate is always ellipse
    real_type xx[2], aa[2], bbeta[2], aalfa[2], y12[2] ;
    xx[0] = inn1 ;
    xx[1] = inn2 ;
    aa[0] = a_min/(1-inn1*inn1) ;
    aa[1] = a_min/(1-inn2*inn2) ;
    bbeta[0] = 2*longway*asin(sqrt((s-c)/(2*aa[0]))) ;
    bbeta[1] = 2*longway*asin(sqrt((s-c)/(2*aa[1]))) ;

    // make 100.4% sure it's in (-1 <= xx <= +1)
    if      ( xx[0] <= -1 ) aalfa[0] = 2*m_pi ;
    else if ( xx[0] >=  1 ) aalfa[0] = 0 ;
    else                    aalfa[0] = 2*acos( xx[0] ) ;
    if      ( xx[1] <= -1 ) aalfa[1] = 2*m_pi ;
    else if ( xx[1] >=  1 ) aalfa[1] = 0 ;
    else                    aalfa[1] = 2*acos( xx[1] ) ;

    // evaluate the time of flight via Lagrange expression
    //y12 = aa*sqrt(aa)*((aalfa - sin(aalfa)) - (bbeta-sin(bbeta)) + 2*m_pi*m) ;
    y12[0] = aa[0]*sqrt(aa[0])*((aalfa[0] - sin(aalfa[0])) - (bbeta[0]-sin(bbeta[0])) + 2*m_pi*m) ;
    y12[1] = aa[1]*sqrt(aa[1])*((aalfa[1] - sin(aalfa[1])) - (bbeta[1]-sin(bbeta[1])) + 2*m_pi*m) ;
    
    // initial estimates for y
    real_type y1, y2 ;
    if ( m == 0 ) {
      y1 = log(y12[0]) - logt;
      y2 = log(y12[1]) - logt;
    } else {
      y1 = y12[0] - tf;
      y2 = y12[1] - tf;
    }

    // Solve for x
    // ---------------------------------------------------------

    // Newton-Raphson iterations
    // NOTE - the number of iterations will go to infinity in case
    // m > 0  and there is no solution. Start the other routine in
    // that case
    real_type err        = 1E100 ;
    int_type  iterations = 0 ;
    real_type xnew       = 0 ;
    while ( err > tol ) {
      // increment number of iterations
      ++iterations ;
      // new x
      xnew = (x1*y2 - y1*x2) / (y2-y1);
      // copy-pasted code (for performance)
      real_type x = m == 0 ? exp(xnew) - 1 : atan(xnew)/m_pi_2 ;
      real_type a = a_min/(1 - x*x) ;
      real_type alfa, beta ;
      if (x < 1) { // ellipse
        beta = longway * 2*asin(sqrt((s-c)/2/a));
        // make 100.4% sure it's in (-1 <= xx <= +1)
        if      ( x <= -1 ) alfa = 2*m_pi ;
        else if ( x >=  1 ) alfa = 0 ;
        else                alfa = 2*acos( x ) ;
      } else { // hyperbola
        alfa = 2*acosh(x);
        beta = longway * 2*asinh(sqrt((s-c)/(-2*a)));
      }
      // evaluate the time of flight via Lagrange expression
      real_type tof ;
      if (a > 0) tof = a*sqrt(a)*((alfa - sin(alfa)) - (beta-sin(beta)) + 2*m_pi*m);
      else       tof = -a*sqrt(-a)*((sinh(alfa) - alfa) - (sinh(beta) - beta));

      // new value of y
      real_type ynew ;
      if ( m ==0 ) ynew = log(tof) - logt ;
      else         ynew = tof - tf ;

      // save previous and current values for the next iterarion
      // (prevents getting stuck between two values)
      x1 = x2 ; x2 = xnew ;
      y1 = y2 ; y2 = ynew ;
      // update error
      err = abs(x1 - xnew) ;

      // escape clause
      if ( iterations > 15 ) {
        bad = true;
        break;
      }
    }
    
    // If the Newton-Raphson scheme failed, try to solve the problem
    // with the other Lambert targeter.
    if ( bad || isnan(err) || isinf(err) ) {
      // NOTE: use the original, UN-normalized quantities
      return lambert_LancasterBlanchard( R1, R2, tf_in, m_in, muC, V1, V2 ) ;
    }

    // convert converged value of x
    real_type x = m==0 ? exp(xnew) - 1 : atan(xnew)/m_pi_2 ;

    /*
    // The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
    // now need the conic. As for transfer angles near to pi the Lagrange-
    // coefficients technique goes singular (dg approaches a zero/zero that is
    // numerically bad) we here use a different technique for those cases. When
    // the transfer angle is exactly equal to pi, then the ih unit vector is not
    // determined. The remaining equations, though, are still valid.
    */

    // Solution for the semi-major axis
    real_type a = a_min/(1-x*x);

    // Calculate psi
    real_type beta, alfa, psi, eta, eta2 ;
    if (x < 1) { // ellipse
      beta = longway * 2*asin(sqrt((s-c)/2/a));
      // make 100.4% sure it's in (-1 <= xx <= +1)
      if      ( x <= -1 ) alfa = 2*m_pi ;
      else if ( x >=  1 ) alfa = 0 ;
      else                alfa = 2*acos( x ) ;
      psi  = (alfa-beta)/2;
      eta2 = 2*a*power2(sin(psi))/s;
      eta  = sqrt(eta2);
    } else { // hyperbola
      beta = longway * 2*asinh(sqrt((c-s)/2/a));
      alfa = 2*acosh(x);
      psi  = (alfa-beta)/2;
      eta2 = -2*a*power2(sinh(psi))/s;
      eta  = sqrt(eta2);
    }

    // unit of the normalized normal vector
    real_type ih[3] = { longway * nrmunit[0],
                        longway * nrmunit[1],
                        longway * nrmunit[2] } ;

    // unit vector for normalized [r2vec]
    real_type r2n[3] = { r2vec[0]/mr2vec,
                         r2vec[1]/mr2vec,
                         r2vec[2]/mr2vec } ;

    // cross-products
    // don't use cross() (emlmex() would try to compile it, and this way it
    // also does not create any additional overhead)
    real_type crsprd1[3], crsprd2[3] ;
    cross(ih,r1vec,crsprd1) ;
    cross(ih,r2n,crsprd2) ;

    // radial and tangential directions for departure velocity
    real_type Vr1 = 1/eta/sqrt(a_min) * (2*Lambda*a_min - Lambda - x*eta);
    real_type Vt1 = sqrt(mr2vec/a_min/eta2 * power2(sin(dth/2)));

    // radial and tangential directions for arrival velocity
    real_type Vt2 = Vt1/mr2vec;
    real_type Vr2 = (Vt1 - Vt2)/tan(dth/2) - Vr1;
      
    // terminal velocities
    V1[0] = (Vr1*r1vec[0] + Vt1*crsprd1[0])*V;
    V1[1] = (Vr1*r1vec[1] + Vt1*crsprd1[1])*V;
    V1[2] = (Vr1*r1vec[2] + Vt1*crsprd1[2])*V;
    V2[0] = (Vr2*r2n[0] + Vt2*crsprd2[0])*V;
    V2[1] = (Vr2*r2n[1] + Vt2*crsprd2[1])*V;
    V2[2] = (Vr2*r2n[2] + Vt2*crsprd2[2])*V;

    return 1 ; // (success)
  }

}
