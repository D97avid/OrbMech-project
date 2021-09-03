function dx = kepl_orbit_J2_moon(x,t,mu_E,mu_M,J2,R_E)
% 
% Function to compute the derivative of orbital state in Cartesian coordinates 
% to propagate the keplerian orbit in Restricted 2BP assumptions taking into account J2 and Moon perturbing effects.
% 
% PROTOTYPE:
%  dx = kepl_orbit_J2_moon(x,mu,J2,R)
% 
% INPUT:
%  x [6,1]       orbital state in cartesian coordinates x=[rx,ry,rz,vx,vy,vz]  [km],[km/s]
%  t [1]         time (MJD2000) expressed in seconds [s]
%  mu_E [1]      gravitational parameter of the Earth   [km^3/s^2]
%  mu_M [1]      gravitational parameter of the Moon   [km^3/s^2]
%  J2 [1]        second zonal harmonic of the Earth   
%  R_E [1]       mean radius of the Earth    [km]
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

r = norm(x(1:3));               % modulus of position vector [km]
a_J2 = (3/2)*(J2/(r^4))*mu_E*(R_E^2)*[(x(1)/r)*(5*((x(3)^2)/r^2)-1);(x(2)/r)*(5*((x(3)^2)/r^2)-1);(x(3)/r)*(5*((x(3)^2)/r^2)-3)]; %perturbing acceleration due to J2 [km/s^2]

t_MJD2000=t/86400;              % time MJD2000 [days]
[rm, ~] = ephMoon(t_MJD2000);   % extract the Moon position vector from the ephemerides [km]
nrm=norm(rm);                   % modulus of the Moon position vector [km]
rr_sc = x(1:3);                 % position vector of the s/c [km]
rms=rm'-rr_sc;                  % position vector of the Moon wrt the s/c [km]
nrms=norm(rms);                 % modulus of Moon position vector wrt the s/c [km]

% Perturbation acceleration due to Moon gravity:
a_moon = mu_M .* ( (rms./(nrms^3) ) - (rm'./(nrm^3)) ); % perturbing acceleration due to the Moon [km/s^2]

% Computation of the derivative of the orbital state
dx = [x(4);x(5);x(6);((-mu_E*x(1))/(r^3))+a_J2(1)+a_moon(1); ... 
                     ((-mu_E*x(2))/(r^3))+a_J2(2)+a_moon(2); ... 
                     ((-mu_E*x(3))/(r^3))+a_J2(3)+a_moon(3)];
