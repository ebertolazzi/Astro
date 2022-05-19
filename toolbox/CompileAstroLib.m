clc;
clear functions;

old_dir = cd(fileparts(which(mfilename)));

NAMES = { 'Astro', 'Lambert', 'Check_EQ_for_consistency' };
 
LIB_NAMES = { ...
  'Astro.cc', ...
  'CPUinfo.cc', ...
  'Console.cc', ...
  'GenericContainer/GenericContainer.cc', ...
  'GenericContainer/GenericContainerCinterface.cc', ...
  'GenericContainer/GenericContainerSerialize.cc', ...
  'GenericContainer/GenericContainerSupport.cc', ...
  'GenericContainer/GenericContainerTables.cc', ...
  'Kepler.cc', ...
  'Malloc.cc', ...
  'Numbers.cc', ...
  'PolynomialRoots-1-Quadratic.cc', ...
  'PolynomialRoots-2-Cubic.cc', ...
  'PolynomialRoots-3-Quartic.cc', ...
  'PolynomialRoots-Jenkins-Traub.cc', ...
  'PolynomialRoots-Utils.cc', ...
  'Rocket.cc', ...
  'Table.cc', ...
  'TicToc.cc', ...
  'Trace.cc', ...
  'Units.cc', ...
  'Utils.cc', ...
  'fmt.cc', ...
  'lambert.cc', ...
  'rang.cc', ...
};

MROOT = matlabroot;
%DEFS  = '-DG2LIB_DEBUG ';
DEFS  = ' ';

CMDBASE = [ 'mex -c -largeArrayDims -Isrc -Isrc/Utils ', DEFS];
if ismac
  CMDBASE = [CMDBASE, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g -F/Library/Frameworks " '];
elseif isunix
  CMDBASE = [CMDBASE, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g" '];
elseif ispc
  CMDBASE = [CMDBASE, 'COMPFLAGS="\$COMPFLAGS -O2" '];
end

LIB_OBJS = '';
for k=1:length(LIB_NAMES)
  [filepath,bname,ext] = fileparts(LIB_NAMES{k});
  NAME = [' src/', filepath, '/', bname, ext ];
  if isunix
    LIB_OBJS = [ LIB_OBJS, bname, '.o ' ];
  elseif ispc
    LIB_OBJS = [ LIB_OBJS, bname, '.obj ' ];
  end
  CMD = [CMDBASE ' -c ' NAME];
  disp('---------------------------------------------------------');
  disp(CMD);
  eval(CMD);
end

for k=1:length(NAMES)
  N=NAMES{k};
  disp('---------------------------------------------------------');
  fprintf(1,'Compiling: %s\n',N);

  CMD = [ 'while mislocked(''' N '''); munlock(''' N '''); end;'];
  eval(CMD);

  CMD = [ 'mex ', DEFS, ' -Isrc -output bin/', N, 'MexWrapper' ];
  CMD = [ CMD, ' -largeArrayDims src_mex/mex_', N ];
  CMD = [ CMD, '.cc src_mex/GenericContainerMatlabInterface.cc ', LIB_OBJS ];

  if ismac
    CMD = [CMD, ' CXXFLAGS="\$CXXFLAGS -Wall -O2 -g -F/Library/Frameworks "'];
    %CMD = [CMD, ' LDFLAGS="\$LDFLAGS -F/Library/Frameworks -framework MechatronixAstro  -framework MechatronixCore "'];
  elseif isunix
    % Workaround for MATLAB 2020 that force dynamic link with old libstdc++
    % solution: link with static libstdc++
    % ARCH  = computer('arch');
    % PATH1 = [MROOT, '/bin/', ARCH];
    % PATH2 = [MROOT, '/extern/bin/', ARCH];
    CMD = [ CMD, ...
      ' CXXFLAGS="\$CXXFLAGS -Wall -O2 -g"' ...
      ' LDFLAGS="\$LDFLAGS -static-libgfortran -static-libstdc++ -static-libgcc"' ...
      ' LINKLIBS="-ldl -L\$MATLABROOT/bin/\$ARCH -L\$MATLABROOT/extern/bin/\$ARCH -lMatlabDataArray -lmx -lmex -lmat -lm "' ...
    ];
  elseif ispc
    CMD = [CMD, 'COMPFLAGS="\$COMPFLAGS -O2" '];
  end

  disp(CMD);
  eval(CMD);
end

if isunix
  delete *.o
else
  delete *.obj
end

cd(old_dir);

disp('----------------------- DONE ----------------------------');
