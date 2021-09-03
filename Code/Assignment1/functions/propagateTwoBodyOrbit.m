function [rr,vv] = propagateTwoBodyOrbit(mu,y0,tspan)
% 
% Function for the ODE propagation of the two body orbit.
% 
% PROTOTYPE:
%  [rr,vv] = propagateTwoBodyOrbit(mu,y0,tspan)
% 
% INPUT:
%  mu [1]        gravitational constant                               [km^3/s^2]
%  y0 [6,1]      initial conditions (position and velocity)           [km/s]
%  tspan []      time vector                                          [km]
% 
% OUTPUT:
%  rr [3,1]      position vector                                      [km]
%  vv [3,1]      velocity vector                                      [km/s]

% 
% CONTRIBUTORS:
%  Marco Adorno
%  Giuseppe Esposito 
%  Davide Gravina 
%  David Reina
% 
% VERSIONS:
%  20-01-2021: First version
% 

options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[~, Y] = ode113( @(t,y) TwoBodyOrbit(tspan,y, mu), tspan, y0, options );

rr = Y(:,1:3);
vv = Y(:,4:6);
