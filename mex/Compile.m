clc;

NAMES = { 'LambertProblem' } ;

SRCS_BASE = '../src/' ;
SRCS = [ 'common_Mex.cc ' ...
         SRCS_BASE 'lambert.cc ' ...
       ] ;

LIBS = ['-largeArrayDims '] ;
INC  = [ '-I' SRCS_BASE ] ;

disp('---------------------------------------------------------');
for k=1:1
  N=NAMES{k} ;
  fprintf(1,'Compiling: %s\n',N) ;

  if isunix
    CMD = [ 'mex ' INC ' -output ' N ' mex' N '.cc ' SRCS ' ' LIBS ] ;
    if ismac
      CMD = [CMD, ...
             ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0 ' ...
             ' -iframework /Library/Frameworks/" ' ...
             ' LDFLAGS="\$LDFLAGS -framework MechatronixCore -framework MechatronixSolver  -framework Accelerate "' ] ;
    else
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0"'] ;
    end
  elseif ispc
    % -llegacy_stdio_definitions
    CMD  = [ 'mex ' INC ' -output ' N ' mex' N '.cc ' SRCS ' ' LIBS ] ;
  end
  disp(CMD);
  eval(CMD);
end
disp('----------------------- DONE ----------------------------');
