function dx = kepl_orbit_J2(x,mu,J2,R)
% 
% Function to compute the derivative of orbital state in Cartesian coordinates 
% to propagate the keplerian orbit in Restricted 2BP assumptions taking into account J2 effect.
% 
% PROTOTYPE:
%  dx = kepl_orbit_J2(x,mu,J2,R)
% 
% INPUT:
%  x [6,1]       orbital state in cartesian coordinates x=[rx,ry,rz,vx,vy,vz]  [km],[km/s]
%  mu [1]        gravitational parameter of primary body    [km^3/s^2]
%  J2 [1]        second zonal harmonic of primary body   
%  R [1]         mean radius of primary body    [km]
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

r = norm([x(1:3)]); % modulus of position vector [km]
a_J2 = (3/2)*(J2/(r^4))*mu*(R^2)*[(x(1)/r)*(5*((x(3)^2)/r^2)-1);(x(2)/r)*(5*((x(3)^2)/r^2)-1);(x(3)/r)*(5*((x(3)^2)/r^2)-3)]; %acceleration vector due to J2 [km/s^2]
dx = [x(4);x(5);x(6);((-mu*x(1))/(r^3))+a_J2(1); ... 
                     ((-mu*x(2))/(r^3))+a_J2(2); ... 
                     ((-mu*x(3))/(r^3))+a_J2(3)];