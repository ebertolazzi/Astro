/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  file: Spacecraft_Main_noRoad.cc                                         |
 |                                                                          |
 |  version: 1.0   date 5/10/2010                                           |
 |                                                                          |
 |  Copyright (C) 2010                                                      |
 |                                                                          |
 |      Enrico Bertolazzi and Francesco Biral                               |
 |      Dipartimento di Ingegneria Meccanica e Strutturale                  |
 |      Universita` degli Studi di Trento                                   |
 |      Via Mesiano 77, I-38050 Trento, Italy                               |
 |      email: enrico.bertolazzi@ing.unitn.it                               |
 |             francesco.biral@ing.unitn.it                                 |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#define BUILDING_DLL

#include <MechatronixInterfaceMatlab/MechatronixInterfaceMatlab.hh>
#include "MechatronixAstro.hh"
#include <string.h>

#include "SpacecraftCompetition.hh"

using namespace CommonLoad ;
using namespace MechatronixAstro ;

/*\
 *                      _____                 _   _             
 *  _ __ ___   _____  _|  ___|   _ _ __   ___| |_(_) ___  _ __  
 * | '_ ` _ \ / _ \ \/ / |_ | | | | '_ \ / __| __| |/ _ \| '_ \ 
 * | | | | | |  __/>  <|  _|| |_| | | | | (__| |_| | (_) | | | |
 * |_| |_| |_|\___/_/\_\_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|
 *                                                              
\*/

static vector<Orbits::Asteroid> asteroids ;

#define GET_ID_AST \
 indexType id = Matlab::getIndex( prhs[1] ) ; \
 ASSERT( id < asteroids.size() && id >= 0, "bad asteroid number:" << id << " must be less than " << asteroids.size() ) ; \
 Orbits::Asteroid & A = asteroids[id]

static
void
getPositions( Orbits::Asteroid const & A,
              VectorOfValues   const & t,
              VectorOfValues         & x,
              VectorOfValues         & y,
              VectorOfValues         & z ) {
  indexType sz = t . size() ;
  x . resize(sz) ;
  y . resize(sz) ;
  z . resize(sz) ;
  for ( indexType i = 0 ; i < sz ; ++i ) 
    A . position( t[i], x[i], y[i], z[i] ) ;
}

static
void
getVelocities( Orbits::Asteroid const & A,
               VectorOfValues   const & t,
               VectorOfValues         & vx,
               VectorOfValues         & vy,
               VectorOfValues         & vz ) {
  indexType sz = t . size() ;
  vx . resize(sz) ;
  vy . resize(sz) ;
  vz . resize(sz) ;
  for ( indexType i = 0 ; i < sz ; ++i ) 
    A . velocity( t[i], vx[i], vy[i], vz[i] ) ;
}

static
void
getAccelerations( Orbits::Asteroid const & A,
                  VectorOfValues   const & t,
                  VectorOfValues         & ax,
                  VectorOfValues         & ay,
                  VectorOfValues         & az ) {
  indexType sz = t . size() ;
  ax . resize(sz) ;
  ay . resize(sz) ;
  az . resize(sz) ;
  for ( indexType i = 0 ; i < sz ; ++i ) 
    A . acceleration( t[i], ax[i], ay[i], az[i] ) ;
}

static
void
getJerks( Orbits::Asteroid const & A,
          VectorOfValues   const & t,
          VectorOfValues         & jx,
          VectorOfValues         & jy,
          VectorOfValues         & jz ) {
  indexType sz = t . size() ;
  jx . resize(sz) ;
  jy . resize(sz) ;
  jz . resize(sz) ;
  for ( indexType i = 0 ; i < sz ; ++i ) 
    A . jerk( t[i], jx[i], jy[i], jz[i] ) ;
}

void
mexFunction( int             nargout,
             mxArray       * plhs[], 
             int             nargin, 
             mxArray const * prhs[] ) {

  static Matlab::mexStreambuf<char> mexBuffer ;
  cout.rdbuf(&mexBuffer);

  // Read first argument: must be a string
  if ( nargin == 0 || !mxIsChar(prhs[0]) )
    mexErrMsgTxt( "First argument should be a string." ) ;

  try {

    string what = mxArrayToString(prhs[0]);

    for ( string::iterator i = what.begin(); i != what.end(); ++i)
      *i = tolower(*i);

    if ( what == "load" ) {

      Matlab::expected( what, 1, nargin, 0, nargout ) ;
      if ( !mxIsChar(prhs[1]) ) mexErrMsgTxt( "Second argument should be a string." ) ;
      string fname = mxArrayToString(prhs[1]);
      Orbits::loadAsteroids( fname.c_str(), asteroids ) ;

    } else if ( what == "numberofasteroids" ) {

      Matlab::expected( what, 0, nargin, 1, nargout ) ;
      Matlab::toMatlab(asteroids.size(),plhs[0]) ;

    } else if ( what == "getname" ) {

      Matlab::expected( what, 1, nargin, 1, nargout ) ;
      GET_ID_AST ;

      /* allocate memory for output string */
      char * output_buf = (char*)mxCalloc(A.theName().size()+1, sizeof(char)) ;
      strcpy( output_buf, A.theName().c_str() ) ;
      /* set C-style string output_buf to MATLAB mexFunction output*/
      plhs[0] = mxCreateString(output_buf);

    } else if ( what == "period" ) {

      Matlab::expected( what, 1, nargin, 1, nargout ) ;
      GET_ID_AST ;
      Matlab::toMatlab(A.period(),plhs[0]) ;

    } else if ( what == "pos+vel" ) {

      Matlab::expected( what, 2, nargin, 6, nargout ) ;
      GET_ID_AST ;
      VectorOfValues t, x, y, z, vx, vy, vz ;
      Matlab::fromMatlab(prhs[2],t) ;
      getPositions ( A, t, x, y, z ) ;
      getVelocities( A, t, vx, vy, vz ) ;
      Matlab::toMatlab(x,plhs[0]) ;
      Matlab::toMatlab(y,plhs[1]) ;
      Matlab::toMatlab(z,plhs[2]) ;
      Matlab::toMatlab(vx,plhs[3]) ;
      Matlab::toMatlab(vy,plhs[4]) ;
      Matlab::toMatlab(vz,plhs[5]) ;

    } else if ( what == "position" ) {

      Matlab::expected( what, 2, nargin, 3, nargout ) ;
      GET_ID_AST ;
      VectorOfValues t, x, y, z ;
      Matlab::fromMatlab(prhs[2],t) ;
      getPositions( A, t, x, y, z ) ;
      Matlab::toMatlab(x,plhs[0]) ;
      Matlab::toMatlab(y,plhs[1]) ;
      Matlab::toMatlab(z,plhs[2]) ;

    } else if ( what == "velocity" ) {

      Matlab::expected( what, 2, nargin, 3, nargout ) ;
      GET_ID_AST ;
      VectorOfValues t, vx, vy, vz ;
      Matlab::fromMatlab(prhs[2],t) ;
      getVelocities( A, t, vx, vy, vz ) ; 
      Matlab::toMatlab(vx,plhs[0]) ;
      Matlab::toMatlab(vy,plhs[1]) ;
      Matlab::toMatlab(vz,plhs[2]) ;

    } else if ( what == "acceleration" ) {

      Matlab::expected( what, 2, nargin, 3, nargout ) ;
      GET_ID_AST ;
      VectorOfValues t, ax, ay, az ;
      Matlab::fromMatlab(prhs[2],t) ;
      getAccelerations( A, t, ax, ay, az ) ; 
      Matlab::toMatlab(ax,plhs[0]) ;
      Matlab::toMatlab(ay,plhs[1]) ;
      Matlab::toMatlab(az,plhs[2]) ;

    } else if ( what == "jerk" ) {

      Matlab::expected( what, 2, nargin, 3, nargout ) ;
      GET_ID_AST ;
      VectorOfValues t, jx, jy, jz ;
      Matlab::fromMatlab(prhs[2],t) ;
      getJerks( A, t, jx, jy, jz ) ; 
      Matlab::toMatlab(jx,plhs[0]) ;
      Matlab::toMatlab(jy,plhs[1]) ;
      Matlab::toMatlab(jz,plhs[2]) ;

    } else if ( what == "orbitalclassic" ) {

      Matlab::expected( what, 1, nargin, 5, nargout ) ;
      GET_ID_AST ;
      valueType e, a, i, Omega, omega ;
      A . getClassicalOrbital( e, a, i, Omega, omega ) ;
      Matlab::toMatlab(e,plhs[0]) ;
      Matlab::toMatlab(a,plhs[1]) ;
      Matlab::toMatlab(i,plhs[2]) ;
      Matlab::toMatlab(Omega,plhs[3]) ;
      Matlab::toMatlab(omega,plhs[4]) ;

    } else if ( what == "invariantl" ) {

      Matlab::expected( what, 5, nargin, 1, nargout ) ;
      valueType const p = Matlab::getValue( prhs[1] ) ;
      valueType const f = Matlab::getValue( prhs[2] ) ;
      valueType const g = Matlab::getValue( prhs[3] ) ;
      valueType const h = Matlab::getValue( prhs[4] ) ;
      valueType const k = Matlab::getValue( prhs[5] ) ;
      VectorOfValues Linvariant(3) ;
      invariantL( p, f, g, h, k, &Linvariant.front() ) ;
      Matlab::toMatlab(Linvariant,plhs[0]) ;

    } else if ( what == "invarianta" ) {

      Matlab::expected( what, 5, nargin, 1, nargout ) ;
      valueType const p = Matlab::getValue( prhs[1] ) ;
      valueType const f = Matlab::getValue( prhs[2] ) ;
      valueType const g = Matlab::getValue( prhs[3] ) ;
      valueType const h = Matlab::getValue( prhs[4] ) ;
      valueType const k = Matlab::getValue( prhs[5] ) ;
      VectorOfValues Ainvariant(3) ;
      invariantA( p, f, g, h, k, &Ainvariant.front() ) ;
      Matlab::toMatlab(Ainvariant,plhs[0]) ;

    } else if ( what == "orbitenergy" ) {

      Matlab::expected( what, 5, nargin, 1, nargout ) ;
      valueType const p = Matlab::getValue( prhs[1] ) ;
      valueType const f = Matlab::getValue( prhs[2] ) ;
      valueType const g = Matlab::getValue( prhs[3] ) ;
      valueType const h = Matlab::getValue( prhs[4] ) ;
      valueType const k = Matlab::getValue( prhs[5] ) ;
      valueType E = orbitEnergy( p, f, g, h, k ) ;
      Matlab::toMatlab(E,plhs[0]) ;

    } else if ( what == "orbitalequinoctial" ) {

      Matlab::expected( what, 1, nargin, 5, nargout ) ;
      GET_ID_AST ;
      valueType p, f, g, h, k ;
      A . getEquinoctialOrbital( p, f, g, h, k ) ;
      Matlab::toMatlab(p,plhs[0]) ;
      Matlab::toMatlab(f,plhs[1]) ;
      Matlab::toMatlab(g,plhs[2]) ;
      Matlab::toMatlab(h,plhs[3]) ;
      Matlab::toMatlab(k,plhs[4]) ;

    } else if ( what == "langle" ) {

      Matlab::expected( what, 2, nargin, 1, nargout ) ;
      GET_ID_AST ;
      valueType t = Matlab::getValue( prhs[2] ) ;
      valueType L = A . Langle( t ) ;
      Matlab::toMatlab(L,plhs[0]) ;
      
    } else if ( what == "classictoequinoctial" ) {

      Matlab::expected( what, 6, nargin, 6, nargout ) ;
      valueType const e     = Matlab::getValue( prhs[1] ) ;
      valueType const a     = Matlab::getValue( prhs[2] ) ;
      valueType const i     = Matlab::getValue( prhs[3] ) ;
      valueType const Omega = Matlab::getValue( prhs[4] ) ;
      valueType const omega = Matlab::getValue( prhs[5] ) ;
      valueType const nu    = Matlab::getValue( prhs[6] ) ;
      valueType p, f, g, h, k, L ;
      Orbits::fromClassicToEquinoctial( e, a, i, Omega, omega, nu,
                                        p, f, g, h, k, L ) ;
      Matlab::toMatlab(p,plhs[0]) ;
      Matlab::toMatlab(f,plhs[1]) ;
      Matlab::toMatlab(g,plhs[2]) ;
      Matlab::toMatlab(h,plhs[3]) ;
      Matlab::toMatlab(k,plhs[4]) ;
      Matlab::toMatlab(L,plhs[5]) ;

    } else if ( what == "equinoctialtoclassic" ) {

      Matlab::expected( what, 6, nargin, 6, nargout ) ;
      valueType const p = Matlab::getValue( prhs[1] ) ;
      valueType const f = Matlab::getValue( prhs[2] ) ;
      valueType const g = Matlab::getValue( prhs[3] ) ;
      valueType const h = Matlab::getValue( prhs[4] ) ;
      valueType const k = Matlab::getValue( prhs[5] ) ;
      valueType const L = Matlab::getValue( prhs[6] ) ;
      valueType e, a, i, Omega, omega, nu ;
      Orbits::fromEquinoctialToClassic( p, f, g, h, k, L,
                                        e, a, i, Omega, omega, nu ) ;
      Matlab::toMatlab(e,    plhs[0]) ;
      Matlab::toMatlab(a,    plhs[1]) ;
      Matlab::toMatlab(i,    plhs[2]) ;
      Matlab::toMatlab(Omega,plhs[3]) ;
      Matlab::toMatlab(omega,plhs[4]) ;
      Matlab::toMatlab(nu,   plhs[5]) ;

    } else if ( what == "equinoctialtocartesian" ) {

      Matlab::expected( what, 6, nargin, 6, nargout ) ;
      valueType const p  = Matlab::getValue( prhs[1] ) ;
      valueType const f  = Matlab::getValue( prhs[2] ) ;
      valueType const g  = Matlab::getValue( prhs[3] ) ;
      valueType const h  = Matlab::getValue( prhs[4] ) ;
      valueType const k  = Matlab::getValue( prhs[5] ) ;
      valueType const L  = Matlab::getValue( prhs[6] ) ;     
      valueType const x  = Orbits::evalX( p, f, g, h, k, L ) ;
      valueType const y  = Orbits::evalY( p, f, g, h, k, L ) ;
      valueType const z  = Orbits::evalZ( p, f, g, h, k, L ) ;
      valueType const vx = Orbits::evalVX( p, f, g, h, k, L ) ;
      valueType const vy = Orbits::evalVY( p, f, g, h, k, L ) ;
      valueType const vz = Orbits::evalVZ( p, f, g, h, k, L ) ;
      Matlab::toMatlab(x,  plhs[0]) ;
      Matlab::toMatlab(y,  plhs[1]) ;
      Matlab::toMatlab(z,  plhs[2]) ;
      Matlab::toMatlab(vx, plhs[3]) ;
      Matlab::toMatlab(vy, plhs[4]) ;
      Matlab::toMatlab(vz, plhs[5]) ;

    } else if ( what == "cartesiantoequinoctial" ) {

      Matlab::expected( what, 6, nargin, 6, nargout ) ;
      valueType const x  = Matlab::getValue( prhs[1] ) ;
      valueType const y  = Matlab::getValue( prhs[2] ) ;
      valueType const z  = Matlab::getValue( prhs[3] ) ;
      valueType const vx = Matlab::getValue( prhs[4] ) ;
      valueType const vy = Matlab::getValue( prhs[5] ) ;
      valueType const vz = Matlab::getValue( prhs[6] ) ;
      
      valueType p, f, g, h, k, L ;
      valueType P[3] = {x,y,z} ;
      valueType V[3] = {vx,vy,vz} ;
      Orbits::evalOrbitalEquinoctial( P, V, p, f, g, h, k, L ) ;
      Matlab::toMatlab(p, plhs[0]) ;
      Matlab::toMatlab(f, plhs[1]) ;
      Matlab::toMatlab(g, plhs[2]) ;
      Matlab::toMatlab(h, plhs[3]) ;
      Matlab::toMatlab(k, plhs[4]) ;
      Matlab::toMatlab(L, plhs[5]) ;

    } else if ( what == "evalclassic" ) {

      Matlab::expected( what, 6, nargin, 6, nargout ) ;
      valueType P[3], V[3] ;
      P[0] = Matlab::getValue( prhs[1] ) ;
      P[1] = Matlab::getValue( prhs[2] ) ;
      P[2] = Matlab::getValue( prhs[3] ) ;
      V[0] = Matlab::getValue( prhs[4] ) ;
      V[1] = Matlab::getValue( prhs[5] ) ;
      V[2] = Matlab::getValue( prhs[6] ) ;
      valueType e, a, i, Omega, omega, nu ;
      Orbits::evalOrbitalClassic( P, V, e, a, i, Omega, omega, nu ) ;
      Matlab::toMatlab(e,    plhs[0]) ;
      Matlab::toMatlab(a,    plhs[1]) ;
      Matlab::toMatlab(i,    plhs[2]) ;
      Matlab::toMatlab(Omega,plhs[3]) ;
      Matlab::toMatlab(omega,plhs[4]) ;
      Matlab::toMatlab(nu,   plhs[5]) ;

    } else if ( what == "evalequinoctial" ) {

      Matlab::expected( what, 6, nargin, 6, nargout ) ;
      valueType P[3], V[3] ;
      P[0] = Matlab::getValue( prhs[1] ) ;
      P[1] = Matlab::getValue( prhs[2] ) ;
      P[2] = Matlab::getValue( prhs[3] ) ;
      V[0] = Matlab::getValue( prhs[4] ) ;
      V[1] = Matlab::getValue( prhs[5] ) ;
      V[2] = Matlab::getValue( prhs[6] ) ;
      valueType p, f, g, h, k, L ;
      Orbits::evalOrbitalEquinoctial( P, V, p, f, g, h, k, L ) ;
      Matlab::toMatlab(p,plhs[0]) ;
      Matlab::toMatlab(f,plhs[1]) ;
      Matlab::toMatlab(g,plhs[2]) ;
      Matlab::toMatlab(h,plhs[3]) ;
      Matlab::toMatlab(k,plhs[4]) ;
      Matlab::toMatlab(L,plhs[5]) ;

    } else if ( what == "evalfrenet" ) {

      VectorOfValues x, y, z, vx, vy, vz ;
      VectorOfValues Drx, Dry, Drz ;
      VectorOfValues Dtx, Dty, Dtz ;
      VectorOfValues Dnx, Dny, Dnz ;

      Matlab::expected( what, 6, nargin, 9, nargout ) ;
      Matlab::fromMatlab(prhs[1],x) ;
      Matlab::fromMatlab(prhs[2],y) ;
      Matlab::fromMatlab(prhs[3],z) ;
      Matlab::fromMatlab(prhs[4],vx) ;
      Matlab::fromMatlab(prhs[5],vy) ;
      Matlab::fromMatlab(prhs[6],vz) ;
      indexType sz = x . size() ;
      if ( y.size()  != sz ||
           z.size()  != sz ||
           vx.size() != sz ||
           vy.size() != sz ||
           vz.size() != sz ) mexErrMsgTxt( "Incompatible dimensions" ) ;
      Drx . reserve(sz) ;
      Dry . reserve(sz) ;
      Drz . reserve(sz) ;
      Dtx . reserve(sz) ;
      Dty . reserve(sz) ;
      Dtz . reserve(sz) ;
      Dnx . reserve(sz) ;
      Dny . reserve(sz) ;
      Dnz . reserve(sz) ;
      for ( indexType i = 0 ; i < sz ; ++i ) {
        valueType P[3] = { x[i], y[i], z[i] } ;
        valueType V[3] = { vx[i], vy[i], vz[i] } ;
        valueType Dr[3], Dt[3], Dn[3] ;
        Orbits::evalFrenet( P, V, Dr, Dt, Dn ) ;
        Drx.push_back(Dr[0]) ;
        Dry.push_back(Dr[1]) ;
        Drz.push_back(Dr[2]) ;
        Dtx.push_back(Dt[0]) ;
        Dty.push_back(Dt[1]) ;
        Dtz.push_back(Dt[2]) ;
        Dnx.push_back(Dn[0]) ;
        Dny.push_back(Dn[1]) ;
        Dnz.push_back(Dn[2]) ;
      }

      Matlab::toMatlab(Drx,plhs[0]) ;
      Matlab::toMatlab(Dry,plhs[1]) ;
      Matlab::toMatlab(Drz,plhs[2]) ;
      Matlab::toMatlab(Dtx,plhs[3]) ;
      Matlab::toMatlab(Dty,plhs[4]) ;
      Matlab::toMatlab(Dtz,plhs[5]) ;
      Matlab::toMatlab(Dnx,plhs[6]) ;
      Matlab::toMatlab(Dny,plhs[7]) ;
      Matlab::toMatlab(Dnz,plhs[8]) ;

    } else if ( what == "integratecartesian" ) {

      VectorOfValues t, Tr, Tt, Tn ;
      VectorOfValues x, y, z, vx, vy, vz, mass ;

      Matlab::expected( what, 11, nargin, 7, nargout ) ;

      valueType Z0[7] ;
      Z0[0] = Matlab::getValue( prhs[1] ) ;
      Z0[1] = Matlab::getValue( prhs[2] ) ;
      Z0[2] = Matlab::getValue( prhs[3] ) ;
      Z0[3] = Matlab::getValue( prhs[4] ) ;
      Z0[4] = Matlab::getValue( prhs[5] ) ;
      Z0[5] = Matlab::getValue( prhs[6] ) ;
      Z0[6] = Matlab::getValue( prhs[7] ) ;
      Matlab::fromMatlab(prhs[8],t) ;
      Matlab::fromMatlab(prhs[9],Tr) ;
      Matlab::fromMatlab(prhs[10],Tt) ;
      Matlab::fromMatlab(prhs[11],Tn) ;

      Orbits::odeTrajSolveXYZ( Z0, t, Tr, Tt, Tn,
                               x, y, z, vx, vy, vz, mass ) ;

      Matlab::toMatlab(x,plhs[0]) ;
      Matlab::toMatlab(y,plhs[1]) ;
      Matlab::toMatlab(z,plhs[2]) ;
      Matlab::toMatlab(vx,plhs[3]) ;
      Matlab::toMatlab(vy,plhs[4]) ;
      Matlab::toMatlab(vz,plhs[5]) ;
      Matlab::toMatlab(mass,plhs[6]) ;

    } else if ( what == "integrateequinoctial" ) {

      VectorOfValues t, Tr, Tt, Tn, p, f, g, h, k, L, mass ;
      Matlab::expected( what, 11, nargin, 7, nargout ) ;
      
      valueType Z0[7] ;
      Z0[0] = Matlab::getValue( prhs[1] ) ;
      Z0[1] = Matlab::getValue( prhs[2] ) ;
      Z0[2] = Matlab::getValue( prhs[3] ) ;
      Z0[3] = Matlab::getValue( prhs[4] ) ;
      Z0[4] = Matlab::getValue( prhs[5] ) ;
      Z0[5] = Matlab::getValue( prhs[6] ) ;
      Z0[6] = Matlab::getValue( prhs[7] ) ;
      Matlab::fromMatlab(prhs[8],t) ;
      Matlab::fromMatlab(prhs[9],Tr) ;
      Matlab::fromMatlab(prhs[10],Tt) ;
      Matlab::fromMatlab(prhs[11],Tn) ;

      Orbits::odeTrajSolveEQU( Z0, t, Tr, Tt, Tn, p, f, g, h, k, L, mass ) ;

      Matlab::toMatlab(p,plhs[0]) ;
      Matlab::toMatlab(f,plhs[1]) ;
      Matlab::toMatlab(g,plhs[2]) ;
      Matlab::toMatlab(h,plhs[3]) ;
      Matlab::toMatlab(k,plhs[4]) ;
      Matlab::toMatlab(L,plhs[5]) ;
      Matlab::toMatlab(mass,plhs[6]) ;

    } else if ( what == "nstepxyz" ) {

      VectorOfValues Z0, ts, t0, Tr, Tt, Tn, Z(7) ;

      Matlab::expected( what, 6, nargin, 7, nargout ) ;

      Matlab::fromMatlab(prhs[1],Z0) ;
      Matlab::fromMatlab(prhs[2],t0) ;
      Matlab::fromMatlab(prhs[3],ts) ;
      Matlab::fromMatlab(prhs[4],Tr) ;
      Matlab::fromMatlab(prhs[5],Tt) ;
      Matlab::fromMatlab(prhs[6],Tn) ;
      
      if ( t0 . size() == 1 ) {
        Orbits::odeNstepXYZ( ts . size(),
                             & ts . front(),
                             & Tr . front(),
                             & Tt . front(),
                             & Tn . front(),
                             t0.front(),
                             & Z0 . front(),
                             & Z . front() ) ;
        Matlab::toMatlab(Z[0],plhs[0]) ;
        Matlab::toMatlab(Z[1],plhs[1]) ;
        Matlab::toMatlab(Z[2],plhs[2]) ;
        Matlab::toMatlab(Z[3],plhs[3]) ;
        Matlab::toMatlab(Z[4],plhs[4]) ;
        Matlab::toMatlab(Z[5],plhs[5]) ;
        Matlab::toMatlab(Z[6],plhs[6]) ;
      } else {
        VectorOfValues x, y, z, vx, vy, vz, mass ;
        Orbits::odeNstepXYZ( ts . size(),
                             ts,
                             Tr,
                             Tt,
                             Tn,
                             t0,
                             & Z0 . front(),
                             x, y, z, vx, vy, vz, mass ) ;
      
        Matlab::toMatlab(x,plhs[0]) ;
        Matlab::toMatlab(y,plhs[1]) ;
        Matlab::toMatlab(z,plhs[2]) ;
        Matlab::toMatlab(vx,plhs[3]) ;
        Matlab::toMatlab(vy,plhs[4]) ;
        Matlab::toMatlab(vz,plhs[5]) ;
        Matlab::toMatlab(mass,plhs[6]) ;
      }

    } else if ( what == "trajectory1" ) {
      
      VectorOfValues PARS0, PARS1, alpha, lambda ;

      Matlab::expected( what, 6, nargin, 1, nargout ) ;
      indexType N     = Matlab::getIndex( prhs[1] ) ;
      valueType Tsize = Matlab::getValue( prhs[2] ) ;
      Matlab::fromMatlab(prhs[3],PARS0) ;
      Matlab::fromMatlab(prhs[4],PARS1) ;
      Matlab::fromMatlab(prhs[5],alpha) ;
      Matlab::fromMatlab(prhs[6],lambda) ;

      MapVectorOfValues res ;
      trajectory( N,
                  Tsize,
                  &PARS0.front(),
                  &PARS1.front(),
                  &alpha.front(), 
                  &lambda.front(),
                  res("time"),
                  res("p"),
                  res("f"),
                  res("g"),
                  res("h"),
                  res("k"),
                  res("L"),
                  res("x"),
                  res("y"),
                  res("z") ) ;

      Matlab::toMatlab(res, plhs[0]) ;

    } else if ( what == "trajectory" ) {
      
      VectorOfValues PARS0, PARS1, alpha, lambda ;

      Matlab::expected( what, 6, nargin, 1, nargout ) ;
      indexType N     = Matlab::getIndex( prhs[1] ) ;
      valueType Tsize = Matlab::getValue( prhs[2] ) ;
      Matlab::fromMatlab(prhs[3],PARS0) ;
      Matlab::fromMatlab(prhs[4],PARS1) ;
      Matlab::fromMatlab(prhs[5],alpha) ;
      Matlab::fromMatlab(prhs[6],lambda) ;

      MapVectorOfValues res ;
      trajectory( N,
                  Tsize,
                  &PARS0.front(),
                  &PARS1.front(),
                  &alpha.front(), 
                  &lambda.front(),
                  res("time"),
                  res("p"),
                  res("f"),
                  res("g"),
                  res("h"),
                  res("k"),
                  res("L"),
                  res("x"),
                  res("y"),
                  res("z"),
                  res("Tx"),
                  res("Ty"),
                  res("Tz") ) ;

      Matlab::toMatlab(res, plhs[0]) ;

    } else if ( what == "trajectorytime" ) {
      
      VectorOfValues PARS0, PARS1, alpha, lambda ;

      Matlab::expected( what, 5, nargin, 1, nargout ) ;
      indexType N = Matlab::getIndex( prhs[1] ) ;
      Matlab::fromMatlab(prhs[2],PARS0) ;
      Matlab::fromMatlab(prhs[3],PARS1) ;
      Matlab::fromMatlab(prhs[4],alpha) ;
      Matlab::fromMatlab(prhs[5],lambda) ;

      MapVectorOfValues res ;
      valueType tres = trajectoryTime( N, 
                                       &PARS0.front(),
                                       &PARS1.front(),
                                       &alpha.front(), 
                                       &lambda.front() ) ;
      Matlab::toMatlab( tres, plhs[0] ) ;

    } else if ( what == "timefromlangle" ) {

      Matlab::expected( what, 3, nargin, 1, nargout ) ;
      GET_ID_AST ;
      valueType tbase = Matlab::getValue( prhs[2] ) ;
      valueType L     = Matlab::getValue( prhs[3] ) ;
      valueType tres  = A . timeFromLangle( tbase, L ) ;
      Matlab::toMatlab( tres, plhs[0] ) ; 
/*
    } else if ( what == "crosstime" ) {

      VectorOfValues PARS0, PARS1, alpha, lambda ;

      Matlab::expected( what, 5, nargin, 1, nargout ) ;

      indexType id1 = Matlab::getIndex( prhs[5] ) ;
      ASSERT( id1 < asteroids.size() && id1 >= 0, "bad asteroid number:" << id1 << " must be less than " << asteroids.size() ) ;
      Orbits::Asteroid & A1 = asteroids[id1] ;
 
      indexType nWind       = Matlab::getIndex( prhs[1] ) ;
      valueType initialMass = Matlab::getIndex( prhs[2] ) ;
      valueType Tbegin      = Matlab::getValue( prhs[3] ) ;
      Matlab::fromMatlab(prhs[4],PARS0) ;
  
      valueType ct = crossTime( nWind, // numero di giri, forward backward
                                initialMass,
                                Tbegin,
                                PARS0[0],
                                PARS0[1],
                                PARS0[2],
                                PARS0[3],
                                PARS0[4],
                                PARS0[5],
                                A1 ) ;

      Matlab::toMatlab( ct, plhs[0] ) ;      
*/
    } else if ( what == "rendezvous" ) {

      Matlab::expected( what, 7, nargin, 2, nargout ) ;

      indexType id1 = Matlab::getIndex( prhs[1] ) ;
      ASSERT( id1 < asteroids.size() && id1 >= 0, "bad asteroid number:" << id1 << " must be less than " << asteroids.size() ) ;
      Orbits::Asteroid & A1 = asteroids[id1] ;

      VectorOfValues PARS ;
 
      valueType Tbegin          = Matlab::getValue( prhs[2] ) ;
      indexType N               = Matlab::getIndex( prhs[3] ) ;
      indexType verbosity       = Matlab::getIndex( prhs[4] ) ;
      Matlab::fromMatlab(prhs[5],PARS) ;
      valueType initialMass     = Matlab::getValue( prhs[6] ) ;
      valueType massConsumption = Matlab::getValue( prhs[7] ) ;

      MapVectorOfValues solution ;
      GenericContainer  info ;
      info . clear() ;
      bool ok = rendezvous( N,
                            verbosity,
                            Tbegin,
                            PARS[0], PARS[1], PARS[2], PARS[3], PARS[4], PARS[5],
                            initialMass, massConsumption,
                            100,
                            500,
                            A1, solution, info ) ;
      //if ( !ok ) mexErrMsgTxt("ASTRO: bad computation.") ;

      Matlab::toMatlab(solution, plhs[0]) ;
      Matlab::toMatlab(info,     plhs[1]) ; 

    } else if ( what == "flyby" ) {

      Matlab::expected( what, 7, nargin, 2, nargout ) ;

      indexType id1 = Matlab::getIndex( prhs[1] ) ;
      ASSERT( id1 < asteroids.size() && id1 >= 0, "bad asteroid number:" << id1 << " must be less than " << asteroids.size() ) ;
      Orbits::Asteroid & A1 = asteroids[id1] ;
 
      VectorOfValues PARS ;
 
      valueType Tbegin          = Matlab::getValue( prhs[2] ) ;
      indexType N               = Matlab::getIndex( prhs[3] ) ;
      indexType verbosity       = Matlab::getIndex( prhs[4] ) ;
      Matlab::fromMatlab(prhs[5],PARS) ;
      valueType initialMass     = Matlab::getValue( prhs[6] ) ;
      valueType massConsumption = Matlab::getValue( prhs[7] ) ;

      MapVectorOfValues solution ;
      GenericContainer  info ;
      info . clear() ;
      bool ok = flyby( N,
                       verbosity,
                       Tbegin,
                       PARS[0], PARS[1], PARS[2], PARS[3], PARS[4], PARS[5],
                       initialMass, massConsumption,
                       50,
                       500,
                       A1, solution, info ) ;
      //if ( !ok ) mexErrMsgTxt("ASTRO: bad computation.") ;

      Matlab::toMatlab(solution, plhs[0]) ;      
      Matlab::toMatlab(info,     plhs[1]) ; 

    } else if ( what == "fromearth" ) {

      Matlab::expected( what, 6, nargin, 2, nargout ) ;

      indexType id1 = Matlab::getIndex( prhs[1] ) ;
      ASSERT( id1 < asteroids.size() && id1 >= 0, "bad asteroid number:" << id1 << " must be less than " << asteroids.size() ) ;
      Orbits::Asteroid & earth = asteroids[0] ;
      Orbits::Asteroid & A1    = asteroids[id1] ;

      valueType Tbegin          = Matlab::getValue( prhs[2] ) ;
      indexType N               = Matlab::getIndex( prhs[3] ) ;
      indexType verbosity       = Matlab::getIndex( prhs[4] ) ;
      valueType initialMass     = Matlab::getValue( prhs[5] ) ;
      valueType massConsumption = Matlab::getValue( prhs[6] ) ;

      MapVectorOfValues solution ;
      GenericContainer  info ;
      info . clear() ;
      bool ok = fromEarth( N,
                           verbosity,
                           initialMass,
                           massConsumption,
                           Tbegin,
                           earth,
                           A1,
                           solution,
                           info ) ;
      //if ( !ok ) mexErrMsgTxt("ASTRO: bad computation.") ;

      Matlab::toMatlab(solution, plhs[0]) ;      
      Matlab::toMatlab(info,     plhs[1]) ; 

    } else if ( what == "fromearth1" ) {

      Matlab::expected( what, 6, nargin, 2, nargout ) ;

      indexType id1 = Matlab::getIndex( prhs[1] ) ;
      ASSERT( id1 < asteroids.size() && id1 >= 0, "bad asteroid number:" << id1 << " must be less than " << asteroids.size() ) ;
      Orbits::Asteroid & earth = asteroids[0] ;
      Orbits::Asteroid & A1    = asteroids[id1] ;

      valueType Tbegin          = Matlab::getValue( prhs[2] ) ;
      indexType N               = Matlab::getIndex( prhs[3] ) ;
      indexType verbosity       = Matlab::getIndex( prhs[4] ) ;
      valueType initialMass     = Matlab::getValue( prhs[5] ) ;
      valueType massConsumption = Matlab::getValue( prhs[6] ) ;

      MapVectorOfValues solution ;
      GenericContainer  info ;
      info . clear() ;
      bool ok = fromEarth1( N,
                            verbosity,
                            initialMass,
                            massConsumption,
                            Tbegin,
                            earth,
                            A1,
                            solution,
                            info ) ;
      //if ( !ok ) mexErrMsgTxt("ASTRO: bad computation.") ;

      Matlab::toMatlab(solution, plhs[0]) ;      
      Matlab::toMatlab(info,     plhs[1]) ; 

    } else {

      mexErrMsgTxt("ASTRO: bad call.") ;

    }

  } catch ( exception const & exc ) {
    mexErrMsgTxt(exc.what()) ;  
  } catch ( ... ) {
    mexErrMsgTxt("ASTRO: unknown error") ;    
  }

  cout . flush() ; // non dimenticare di fare il flush!
}
