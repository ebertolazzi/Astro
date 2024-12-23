AU_to_km = 1.49597870691E8;
Y_to_s   = 86400*365.25;
mu_SUN   = 1.32712440018E11; %km^3/s^2
% to UA^3/Y^2
mu_SUN   = mu_SUN/(AU_to_km^3/Y_to_s^2);

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

data_Venus.muS = mu_SUN;
data_Earth.muS = mu_SUN;
data_Mars.muS  = mu_SUN;

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

hold on
P={Venus,Earth,Mars};
for k=1:3
  tt      = 0:P{k}.period/400:P{k}.period;
  [x,y,z] = P{k}.position( tt );
  plot3( x, y, z, 'LineWidth', 3 );
end

axis equal;

t0 = 0;
t1 = 300;
dt = 0.5;

[PE,VE] = Earth.get_PV(t0);
[PM,VM] = Mars.get_PV(t0);

[ta,tb,DV] = LambertMinDeltaV( mu_SUN, PE, VE, PM, VM, t0, t1, dt );

fprintf('ta = %g, tb = %g\n', ta, tb);

fprintf('(AU/Y) DV1 = %g DV2 = %g\n',DV(1),DV(2));
DV = DV*AU_to_km/Y_to_s;
fprintf('(km/s) DV1 = %g DV2 = %g\n',DV(1),DV(2));

fprintf('done\n');
