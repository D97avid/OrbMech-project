function dx = kepl_orbit(x,mu)
% 
% Function to compute the derivative of orbital state in Cartesian coordinates 
% to propagate the keplerian orbit in Restricted 2BP assumptions.
% 
% PROTOTYPE:
%  dx = kepl_orbit(x,mu)
% 
% INPUT:
%  x [6,1]       orbital state in cartesian coordinates x=[rx,ry,rz,vx,vy,vz]  [km],[km/s]
%  mu [1]        gravitational parameter of primary body    [km^3/s^2]
% 
% OUTPUT:
%  dx [6,1]         derivative of orbital state in cartesian coordinates dx=[vx,vy,vz,ax,ay,az]  [km/s],[km/s^2]
%
% CONTRIBUTORS:
%  Marco Adorno
%  Giuseppe Esposito 
%  Davide Gravina 
%  David Reina
% 
% VERSIONS:
%  20-01-2021: First version

r = norm([x(1:3)]); %modulus of position vector [km]
dx = [x(4);x(5);x(6);-mu*x(1)/(r^3);-mu*x(2)/(r^3);-mu*x(3)/(r^3)]; %derivative of orbital state
