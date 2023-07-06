addpath('mArrow3')

close all;

MJD_begin = 64328;

AU_to_km        = 1.49597870691E8;
day_to_s        = 86400;
Y_to_day        = 365.25;
mu_SUN_km3_s2   = 1.32712440018E11; %km^3/s^2
mu_SUN_UA3_day2 = mu_SUN_km3_s2/(AU_to_km^3/day_to_s^2); % to UA^3/day^2

fprintf('load asteroids\n');
lst = load_asteroids2( 'asteroids/GTOC12_Asteroids_Data.txt', mu_SUN_UA3_day2 );

fprintf('plot asteroids\n');
hold on;
A1 = lst{5};
A2 = lst{47439}; %{2345}; %{1325}; %{32803}; %lst{80};

fprintf('plot Earth, Venus, Mars\n');

%
% Data taken from:
%
% https://ssd.jpl.nasa.gov/txt/p_elem_t1.txt

data_Venus.name = 'Venus';
data_Earth.name = 'Earth';
data_Mars.name  = 'Mars';

data_Venus.t0 = 54000; %MJD
data_Earth.t0 = 54000; %MJD
data_Mars.t0  = 54000; %MJD

data_Venus.muS = mu_SUN_UA3_day2;
data_Earth.muS = mu_SUN_UA3_day2;
data_Mars.muS  = mu_SUN_UA3_day2;

data_Venus.a = 0.72333566;
data_Earth.a = 1.00000261;
data_Mars.a  = 1.52371034;

data_Venus.e = 0.00677672;
data_Earth.e = 0.01671123;
data_Mars.e  = 0.09339410;

data_Venus.i = 3.39467605*(pi/180);
data_Earth.i = -0.00001531*(pi/180);
data_Mars.i  = 1.84969142*(pi/180);

data_Venus.L = 181.97909950*(pi/180);
data_Earth.L = 100.46457166*(pi/180);
data_Mars.L  = -4.55343205*(pi/180);

% long.peri.
data_Venus.omega_bar = 131.60246718*(pi/180);
data_Earth.omega_bar = 102.93768193*(pi/180);
data_Mars.omega_bar  = -23.94362959*(pi/180);

% long.node.
data_Venus.Omega = 76.67984255*(pi/180);
data_Earth.Omega = 0*(pi/180);
data_Mars.Omega  = 49.55953891*(pi/180);

data_Venus.omega = data_Venus.omega_bar - data_Venus.Omega;
data_Earth.omega = data_Earth.omega_bar - data_Earth.Omega;
data_Mars.omega  = data_Mars.omega_bar  - data_Mars.Omega;

data_Venus.M0 = data_Venus.L - data_Venus.omega_bar;
data_Earth.M0 = data_Earth.L - data_Earth.omega_bar;
data_Mars.M0  = data_Mars.L  - data_Mars.omega_bar;

fprintf('\n\n');
Venus = Astro();
Venus.setup_Keplerian( data_Venus );
Venus.print();

fprintf('\n\n');
Earth = Astro();
Earth.setup_Keplerian( data_Earth );
Earth.print();

fprintf('\n\n');
Mars = Astro();
Mars.setup_Keplerian( data_Mars );
Mars.print();



m   = 0;
muS = A1.muS;
t0  = MJD_begin;
P1  = A1.position(t0); V1 = A1.velocity(t0);
P2  = A2.position(t0); V2 = A2.velocity(t0);
%
%
%
t_begin     = MJD_begin;
t_end       = MJD_begin+15*365;
t_tolerance = 0.01;
maxDV       = 6/(AU_to_km/day_to_s);
[v_t_begin,v_t_end,v_P1,v_V1,v_W1,v_P2,v_V2,v_W2,v_DV1,v_DV2] = globalMinimumDeltaV(muS,t0,P1,V1,P2,V2,t_begin,t_end,t_tolerance,maxDV);

A = Astro();

tiledlayout(1,3, 'Padding', 'none', 'TileSpacing', 'compact');
for kkk=1:length(v_t_begin)
  nexttile;

  hold on
  P={Venus,Earth,Mars};
  hold on
  for k=1:3
    tt      = 0:P{k}.period/400:P{k}.period;
    [x,y,z] = P{k}.position( tt );
    plot3( x, y, z, 'LineWidth', 3 );
  end

  m = 0;
  t1 = v_t_begin(kkk);
  t2 = v_t_end(kkk);
  tt = t1:1:t2;

  [x,y,z] = A1.position( tt );
  plot3( x, y, z, 'LineWidth', 2, 'Color', 'blue' );

  [x,y,z] = A2.position( tt );
  plot3( x, y, z, 'LineWidth', 2, 'Color', 'green' );

  PP1 = v_P1(:,kkk); PP2 = v_P2(:,kkk);
  VV1 = v_V1(:,kkk); VV2 = v_V2(:,kkk);
  WW1 = v_W1(:,kkk); WW2 = v_W2(:,kkk);

  %[V1,V2,ok] = Lambert( PP1, PP2, t2-t1, m, muS );
  %norm(V1-v_W1(:,kkk))
  %norm(V2-v_W2(:,kkk))

  DV1 = v_DV1(kkk);
  DV2 = v_DV2(kkk);

  A.setup_PV( 'prova', PP1, WW1, t1, muS );
  DV1_kms = DV1*(AU_to_km/day_to_s);
  DV2_kms = DV2*(AU_to_km/day_to_s);
  fprintf('[days]   t1 = %g t2 = %g DT = %g\n', t1-MJD_begin, t2-MJD_begin, t2-t1 );
  fprintf('[UA/day] DV1 = %g DV2 = %g\n', DV1, DV2 );
  fprintf('[km/s]   DV1 = %g DV2 = %g\n\n\n', DV1_kms, DV2_kms );

  if true
    mArrow3( PP1, PP1+VV1./norm(VV1),'color','red','stemWidth',0.01,'facealpha',0.5 );
    mArrow3( PP2, PP2+VV2./norm(VV2),'color','red','stemWidth',0.01,'facealpha',0.5 );

    mArrow3( PP1, PP1+WW1./norm(WW1),'color','blue','stemWidth',0.01,'facealpha',0.5 );
    mArrow3( PP2, PP2+WW2./norm(WW2),'color','blue','stemWidth',0.01,'facealpha',0.5 );
  end

  [x,y,z] = A.position( tt );
  plot3( x, y, z, 'LineWidth', 3, 'Color', 'red' );

  axis equal;

end

fprintf('done\n');
