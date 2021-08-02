AU_to_km = 1.49597870691E8;
Y_to_s   = 86400*365.25;
mu_SUN   = 1.32712440018E11; %km^3/s^2
% to UA^3/Y^2
mu_SUN   = mu_SUN/(AU_to_km^3/Y_to_s^2);

%
% Data taken from:
%
% https://ssd.jpl.nasa.gov/txt/p_elem_t1.txt

data_Mercury.name = 'Mercury';
data_Venus.name   = 'Venus';
data_Earth.name   = 'Earth';
data_Mars.name    = 'Mars';
data_Jupiter.name = 'Jupiter';
data_Saturn.name  = 'Saturn';
data_Uranus.name  = 'Uranus';
data_Neptune.name = 'Neptune';
data_Pluto.name   = 'Pluto';

data_Mercury.t0 = 54000; %MJD
data_Venus.t0   = 54000; %MJD
data_Earth.t0   = 54000; %MJD
data_Mars.t0    = 54000; %MJD
data_Jupiter.t0 = 54000; %MJD
data_Saturn.t0  = 54000; %MJD
data_Uranus.t0  = 54000; %MJD
data_Neptune.t0 = 54000; %MJD
data_Pluto.t0   = 54000; %MJD

data_Mercury.muS = mu_SUN;
data_Venus.muS   = mu_SUN;
data_Earth.muS   = mu_SUN;
data_Mars.muS    = mu_SUN;
data_Jupiter.muS = mu_SUN;
data_Saturn.muS  = mu_SUN;
data_Uranus.muS  = mu_SUN;
data_Neptune.muS = mu_SUN;
data_Pluto.muS   = mu_SUN;

data_Mercury.a = 0.38709927;
data_Venus.a   = 0.72333566;
data_Earth.a   = 1.00000261;
data_Mars.a    = 1.52371034;
data_Jupiter.a = 5.20288700;
data_Saturn.a  = 9.53667594;
data_Uranus.a  = 19.18916464;
data_Neptune.a = 30.06992276;
data_Pluto.a   = 39.48211675;

data_Mercury.e = 0.20563593;
data_Venus.e   = 0.00677672;
data_Earth.e   = 0.01671123;
data_Mars.e    = 0.09339410;
data_Jupiter.e = 0.04838624;
data_Saturn.e  = 0.05386179;
data_Uranus.e  = 0.04725744;
data_Neptune.e = 0.00859048;
data_Pluto.e   = 0.24882730;

data_Mercury.i = 7.00497902*(pi/180);
data_Venus.i   = 3.39467605*(pi/180);
data_Earth.i   = -0.00001531*(pi/180);
data_Mars.i    = 1.84969142*(pi/180);
data_Jupiter.i = 1.30439695*(pi/180);
data_Saturn.i  = 2.48599187*(pi/180);
data_Uranus.i  = 0.77263783*(pi/180);
data_Neptune.i = 1.77004347*(pi/180);
data_Pluto.i   = 17.14001206*(pi/180);

data_Mercury.L = 252.25032350*(pi/180);
data_Venus.L   = 181.97909950*(pi/180);
data_Earth.L   = 100.46457166*(pi/180);
data_Mars.L    = -4.55343205*(pi/180);
data_Jupiter.L = 34.39644051*(pi/180);
data_Saturn.L  = 49.95424423*(pi/180);
data_Uranus.L  = 313.23810451*(pi/180);
data_Neptune.L = -55.12002969*(pi/180);
data_Pluto.L   = 238.92903833*(pi/180);

% long.peri.
data_Mercury.omega_bar = 77.45779628*(pi/180);
data_Venus.omega_bar   = 131.60246718*(pi/180);
data_Earth.omega_bar   = 102.93768193*(pi/180);
data_Mars.omega_bar    = -23.94362959*(pi/180);
data_Jupiter.omega_bar = 14.72847983*(pi/180);
data_Saturn.omega_bar  = 92.59887831*(pi/180);
data_Uranus.omega_bar  = 170.95427630*(pi/180);
data_Neptune.omega_bar = 44.96476227*(pi/180);
data_Pluto.omega_bar   = 224.06891629*(pi/180);

% long.node.
data_Mercury.Omega = 48.33076593*(pi/180);
data_Venus.Omega   = 76.67984255*(pi/180);
data_Earth.Omega   = 0*(pi/180);
data_Mars.Omega    = 49.55953891*(pi/180);
data_Jupiter.Omega = 100.47390909*(pi/180);
data_Saturn.Omega  = 113.66242448*(pi/180);
data_Uranus.Omega  = 74.01692503*(pi/180);
data_Neptune.Omega = 131.78422574*(pi/180);
data_Pluto.Omega   = 110.30393684*(pi/180);

data_Mercury.omega = data_Mercury.omega_bar - data_Mercury.Omega;
data_Venus.omega   = data_Venus.omega_bar   - data_Venus.Omega;
data_Earth.omega   = data_Earth.omega_bar   - data_Earth.Omega;
data_Mars.omega    = data_Mars.omega_bar    - data_Mars.Omega;
data_Jupiter.omega = data_Jupiter.omega_bar - data_Jupiter.Omega;
data_Saturn.omega  = data_Saturn.omega_bar  - data_Saturn.Omega;
data_Uranus.omega  = data_Uranus.omega_bar  - data_Uranus.Omega;
data_Neptune.omega = data_Neptune.omega_bar - data_Neptune.Omega;
data_Pluto.omega   = data_Pluto.omega_bar   - data_Pluto.Omega;

data_Mercury.M0 = data_Mercury.L - data_Mercury.omega_bar;
data_Venus.M0   = data_Venus.L   - data_Venus.omega_bar;
data_Earth.M0   = data_Earth.L   - data_Earth.omega_bar;
data_Mars.M0    = data_Mars.L    - data_Mars.omega_bar;
data_Jupiter.M0 = data_Jupiter.L - data_Jupiter.omega_bar;
data_Saturn.M0  = data_Saturn.L  - data_Saturn.omega_bar;
data_Uranus.M0  = data_Uranus.L  - data_Uranus.omega_bar;
data_Neptune.M0 = data_Neptune.L - data_Neptune.omega_bar;
data_Pluto.M0   = data_Pluto.L   - data_Pluto.omega_bar;

fprintf('\n\n');
Mercury = Astro();
Mercury.setup_Keplerian( data_Mercury );
Mercury.print();

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

fprintf('\n\n');
Jupiter = Astro();
Jupiter.setup_Keplerian( data_Jupiter );
Jupiter.print();

fprintf('\n\n');
Saturn = Astro();
Saturn.setup_Keplerian( data_Saturn );
Saturn.print();

fprintf('\n\n');
Uranus = Astro();
Uranus.setup_Keplerian( data_Uranus );
Uranus.print();

fprintf('\n\n');
Neptune = Astro();
Neptune.setup_Keplerian( data_Saturn );
Neptune.print();

fprintf('\n\n');
Pluto = Astro();
Pluto.setup_Keplerian( data_Pluto );
Pluto.print();


hold on
P={Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune,Pluto};
for k=1:9
  Tmax = 100;
  tt   = 0:P{k}.period/1000:P{k}.period;
  [x,y,z] = P{k}.position( tt );
  plot3( x, y, z );
end

Tmax = 0.25;
tt   = 0:Tmax/1000:Tmax;
for k=1:9
  [x,y,z] = P{k}.position( tt );
  plot3( x, y, z, 'LineWidth', 3 );
end

axis equal;