AU_to_km = 1.49597870691E8;
Y_to_s   = 86400*365.25;
mu_SUN   = 1.32712440018E11; %km^3/s^2
% to UA^3/Y^2
mu_SUN   = mu_SUN/(AU_to_km^3/Y_to_s^2);

lst = load_asteroids( 'asteroids/ast_ephem_gtoc5.txt', mu_SUN );

for k=1:length(lst)
  A  = lst{k};
  tt = 0:A.period/1000:A.period;
  [x,y,z] = A.position( tt );
  plot3( x, y, z );
  hold on;
end

axis equal;