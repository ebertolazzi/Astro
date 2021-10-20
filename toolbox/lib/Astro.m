classdef Astro < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = protected, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    function self = Astro( varargin )
      self.objectHandle = AstroMexWrapper( 'new', varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
      %% Destroy the C++ class instance
      if self.objectHandle ~= 0
        AstroMexWrapper( 'delete', self.objectHandle );
      end
      self.objectHandle = 0; % avoid double destruction of object
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function change_handle( self, new_objectHandle )
      AstroMexWrapper( 'new_handle', self.objectHandle , new_objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Return the `pointer` of the interbal stored c++ object
    %>
    %> **Usage**
    %>
    %> \rst
    %> .. code-block:: matlab
    %>
    %>   obj = ref.obj_handle();
    %> 
    %> \endrst
    %>
    function obj = obj_handle( self )
      obj = self.objectHandle;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Make of copy of a curve object
    %>
    %> **Usage**
    %>
    %> \rst
    %> .. code-block:: matlab
    %>
    %>   ref.copy( C );
    %> 
    %> \endrst
    %>
    %> where `C` id the curve object to be copied.
    %>
    function copy( self, C )
      AstroMexWrapper( 'copy', self.objectHandle, C.obj_handle() );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = name( self )
      res = AstroMexWrapper( 'name', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = position( self, t )
      if nargout == 1
        varargout{1} = AstroMexWrapper( 'position', self.objectHandle, t );
      else
        [varargout{1},varargout{2},varargout{3}] = ...
          AstroMexWrapper( 'position', self.objectHandle, t );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = velocity( self, t )
      if nargout == 1
        varargout{1} = AstroMexWrapper( 'velocity', self.objectHandle, t );
      else
        [varargout{1},varargout{2},varargout{3}] = ...
          AstroMexWrapper( 'velocity', self.objectHandle, t );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = acceleration( self, t )
      if nargout == 1
        varargout{1} = AstroMexWrapper( 'acceleration', self.objectHandle, t );
      else
        [varargout{1},varargout{2},varargout{3}] = ...
          AstroMexWrapper( 'acceleration', self.objectHandle, t );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = jerk( self, t )
      if nargout == 1
        varargout{1} = AstroMexWrapper( 'jerk', self.objectHandle, t );
      else
        [varargout{1},varargout{2},varargout{3}] = ...
          AstroMexWrapper( 'jerk', self.objectHandle, t );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % name, t0, a, e, Omega, omega, iorb, M0, muS
    function setup_Keplerian( self, varargin )
      AstroMexWrapper( 'setup_Keplerian', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % name, t0, p, f, g, h, k, retrograde, L, muS
    function setup_Equinoctial( self, varargin )
      AstroMexWrapper( 'setup_Equinoctial', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setup_PV( self, name, P, V, t0, muS )
      AstroMexWrapper( 'setup_PV', self.objectHandle, ...
        name, P, V, t0, muS ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ma = mean_anomaly( self, t )
      ma = AstroMexWrapper( 'mean_anomaly', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ta = true_anomaly( self, t )
      ta = AstroMexWrapper( 'true_anomaly', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = p_orbital( self )
      [varargout{1:nargout}] = AstroMexWrapper( 'p_orbital', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = f_orbital( self )
      [varargout{1:nargout}] = AstroMexWrapper( 'f_orbital', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = g_orbital( self )
      [varargout{1:nargout}] = AstroMexWrapper( 'g_orbital', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = h_orbital( self )
      [varargout{1:nargout}] = AstroMexWrapper( 'h_orbital', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = k_orbital( self )
      [varargout{1:nargout}] = AstroMexWrapper( 'k_orbital', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = L_orbital( self, varargin )
      [varargout{1:nargout}] = AstroMexWrapper( 'L_orbital', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function L_D = L_orbital_D( self, t )
      L_D = AstroMexWrapper( 'L_orbital_D', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function L_DD = L_orbital_DD( self, t )
      L_DD = AstroMexWrapper( 'L_orbital_DD', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function lp = latitude_of_periapsis( self )
      lp = AstroMexWrapper( 'latitude_of_periapsis', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function la = latitude_of_apoapsis( self )
      la = AstroMexWrapper( 'latitude_of_apoapsis', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = info( self )
      res = AstroMexWrapper( 'info', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function print_info( self )
      AstroMexWrapper( 'print_info', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function n = number_of_revolution( self, t )
      n = AstroMexWrapper( 'number_of_revolution', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function tm = time_from_L_angle( self, t, L )
      tm = AstroMexWrapper( 'time_from_L_angle', self.objectHandle, t, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function P = absolute_position( self, t )
      P = AstroMexWrapper( 'absolute_position', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function V = absolute_velocity( self, t )
      V = AstroMexWrapper( 'absolute_velocity', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function p = period( self )
      p = AstroMexWrapper( 'period', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function a = apoapsis( self )
      a = AstroMexWrapper( 'apoapsis', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function p = periapsis( self )
      p = AstroMexWrapper( 'periapsis', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = muS( self )
      res = AstroMexWrapper( 'muS', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = orbit_energy( self )
      res = AstroMexWrapper( 'orbit_energy', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = radius_by_L( self, L )
      res = AstroMexWrapper( 'radius_by_L', self.objectHandle, L  );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = radius_by_L_D( self, L )
      res = AstroMexWrapper( 'radius_by_L_D', self.objectHandle, L  );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = radius_by_L_DD( self, L )
      res = AstroMexWrapper( 'radius_by_L_DD', self.objectHandle, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = M0( self )
      res = AstroMexWrapper( 'M0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = M0_EQ_gradient( self )
      res = AstroMexWrapper( 'M0_EQ_gradient', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = E0_angle( self )
      res = AstroMexWrapper( 'E0_angle', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = E0_EQ_gradient( self )
      res = AstroMexWrapper( 'E0_EQ_gradient', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = E_angle( self, t )
      res = AstroMexWrapper( 'E_angle', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = E_EQ_gradient( self, t )
      res = AstroMexWrapper( 'E_EQ_gradient', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = H0_angle( self )
      res = AstroMexWrapper( 'H0_angle', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = H0_EQ_gradient( self )
      res = AstroMexWrapper( 'H0_EQ_gradient', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = H_angle( self, t )
      res = AstroMexWrapper( 'H_angle', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = H_EQ_gradient( self, t )
      res = AstroMexWrapper( 'H_EQ_gradient', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = theta0( self )
      res = AstroMexWrapper( 'theta0', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = theta0_EQ_gradient( self )
      res = AstroMexWrapper( 'theta0_EQ_gradient', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = radius_EQ_gradient( self, t )
      res = AstroMexWrapper( 'radius_EQ_gradient', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = position0_EQ_jacobian( self )
      res = AstroMexWrapper( 'position0_EQ_jacobian', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = position_EQ_jacobian( self, t )
      res = AstroMexWrapper( 'position_EQ_jacobian', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = absolute_velocity_by_angle( self, L )
      res = AstroMexWrapper( 'absolute_velocity_by_angle', self.objectHandle, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = absolute_velocity_EQ_gradient( self, t )
      res = AstroMexWrapper( 'absolute_velocity_EQ_gradient', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = velocity0_EQ_jacobian( self )
      res = AstroMexWrapper( 'velocity0_EQ_jacobian', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = velocity_EQ_jacobian( self, t )
      res = AstroMexWrapper( 'velocity_EQ_jacobian', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = L_from_true_anomaly( self, nu )
      res = AstroMexWrapper( 'L_from_true_anomaly', self.objectHandle, nu );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_retrograde( self )
      AstroMexWrapper( 'make_retrograde', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_not_retrograde( self )
      AstroMexWrapper( 'make_not_retrograde', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = retrograde( self )
      res = AstroMexWrapper( 'retrograde', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function N = normal( self )
      N = AstroMexWrapper( 'normal', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = local_frame( self, t )
      F = AstroMexWrapper( 'local_frame', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = local_frame_by_L( self, L )
      F = AstroMexWrapper( 'local_frame_by_L', self.objectHandle, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = ellipse_frame( self )
      F = AstroMexWrapper( 'ellipse_frame', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_E( self, t )
      [varargout{1:nargout}] = AstroMexWrapper( 'eval_E', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_L( self, t )
      [varargout{1:nargout}] = AstroMexWrapper( 'eval_L', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function print( self )
      res = AstroMexWrapper( 'info', self.objectHandle );
      fprintf('\n''%s''\n',self.name());
      fprintf('Equinoctial\n');
      fprintf('  p = %g\n',res.Equinoctial.p);
      fprintf('  f = %g\n',res.Equinoctial.f);
      fprintf('  g = %g\n',res.Equinoctial.g);
      fprintf('  h = %g\n',res.Equinoctial.h);
      fprintf('  k = %g\n',res.Equinoctial.k);
      fprintf('  retrograde = %g\n',res.Equinoctial.retrograde);
      fprintf('\nKepler\n');
      fprintf('  a     = %g\n',res.Kepler.a);
      fprintf('  e     = %g\n',res.Kepler.e);
      fprintf('  i     = %g [%g degree]\n',res.Kepler.i,res.Kepler.i*(180/pi));
      fprintf('  omega = %g [%g degree]\n',res.Kepler.omega,res.Kepler.omega*(180/pi));
      fprintf('  Omega = %g [%g degree]\n',res.Kepler.Omega,res.Kepler.Omega*(180/pi));
      fprintf('\nINFO\n');
      fprintf('  M0 = %g [%g degree]\n',res.M0,res.M0*(180/pi));
      fprintf('  t0 = %g\n',res.t0);
      fprintf('  period       = %g\n',self.period());
      fprintf('  apoapsis     = %g\n',self.apoapsis());
      fprintf('  periapsis    = %g\n',self.periapsis());
      fprintf('  muS          = %g\n',self.muS());
      fprintf('  orbit_energy = %g\n',self.orbit_energy());
      fprintf('\n\n');
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, varargin )
      year = self.period();
      t    = 0:year/1000:year;
      [ X, Y, Z ] = self.position( t );
      plot3( X, Y, Z, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
