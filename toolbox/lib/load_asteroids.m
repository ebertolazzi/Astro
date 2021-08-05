%
%
%
%
function lst = load_asteroids( fname, muS )
  fd = fopen( fname, 'r' );
  % skip first 2 lines
  tline = fgets(fd);
  tline = fgets(fd);
  kkk   = 0;
  while ~feof(fd)
    tline      = fgets(fd);
    ss         = split( tline, '''' );
    pars       = textscan(ss{3},'%f %f %f %f %f %f %f');
    data.name  = ss{2};
    data.t0    = pars{1};
    data.muS   = muS;
    data.a     = pars{2};
    data.e     = pars{3};
    data.i     = pars{4}*(pi/180);
    data.Omega = pars{5}*(pi/180);
    data.omega = pars{6}*(pi/180);
    data.M0    = pars{7}*(pi/180);
    kkk        = kkk+1;
    lst{kkk}   = Astro();
    lst{kkk}.setup_Keplerian( data );
    %lst{kkk}.print();
  end
end
