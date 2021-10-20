AU_to_km = 1.49597870691E8;
Y_to_s   = 86400*365.25;
mu_SUN   = 1.32712440018E11; %km^3/s^2
% to UA^3/Y^2
mu_SUN   = mu_SUN/(AU_to_km^3/Y_to_s^2);

%
% Data taken from:
%
% https://ssd.jpl.nasa.gov/txt/p_elem_t1.txt

data.name = 'Earth';
data.t0   = 54000; %MJD
data.muS  = mu_SUN;
data.p    = 0.999723;
data.f    = -0.0037415;
data.g    = 0.016287;
data.h    = 0,3;
data.k    = -0;
data.L0   = 1.1;
data.retrograde = false;

fprintf('\n\n');
A = Astro();
A.setup_Equinoctial( data );
%A.print();

A1 = Astro();
G  = A.M0_EQ_gradient();
M0 = A.M0();

fprintf('G = %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n', ...
        G(1), G(2), G(3), G(4), G(5), G(6) );

for delta=[1e-3,1e-6,1e-9,1e-12]

  data1   = data;
  data1.p = data.p+delta;
  A1.setup_Equinoctial( data1 );
  pp = A1.M0;

  data1   = data;
  data1.p = data.p-delta;
  A1.setup_Equinoctial( data1 );
  pm = A1.M0;

  Dp = (pp-pm)/(2*delta);

  data1   = data;
  data1.f = data.f+delta;
  A1.setup_Equinoctial( data1 );
  fp = A1.M0;

  data1   = data;
  data1.f = data.f-delta;
  A1.setup_Equinoctial( data1 );
  fm = A1.M0;

  Df = (fp-fm)/(2*delta);

  data1   = data;
  data1.g = data.g+delta;
  A1.setup_Equinoctial( data1 );
  gp = A1.M0;

  data1   = data;
  data1.g = data.g-delta;
  A1.setup_Equinoctial( data1 );
  gm = A1.M0;

  Dg = (gp-gm)/(2*delta);

  data1   = data;
  data1.h = data.h+delta;
  A1.setup_Equinoctial( data1 );
  hp = A1.M0;

  data1   = data;
  data1.h = data.h-delta;
  A1.setup_Equinoctial( data1 );
  hm = A1.M0;
 
  Dh = (hp-hm)/(2*delta);

  data1   = data;
  data1.k = data.k+delta;
  A1.setup_Equinoctial( data1 );
  kp = A1.M0;
 
  data1   = data;
  data1.k = data.k-delta;
  A1.setup_Equinoctial( data1 );
  km = A1.M0;

  Dk = (kp-km)/(2*delta);

  data1    = data;
  data1.L0 = data.L0+delta;
  A1.setup_Equinoctial( data1 );
  Lp = A1.M0;

  data1    = data;
  data1.L0 = data.L0-delta;
  A1.setup_Equinoctial( data1 );
  Lm = A1.M0;

  DL = (Lp-Lm)/(2*delta);

  fprintf('G = %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n', Dp, Df, Dg, Dh, Dk, DL );
end
