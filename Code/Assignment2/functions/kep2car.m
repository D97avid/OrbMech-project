function [rr,vv] = kep2car(kep,mu)
% 
% Function to compute the position and velocity vectors from the keplerian
% elements vector.
% 
% PROTOTYPE:
%  [rr,vv] = kep2car(kep,mu)
% 
% INPUT:
%  kep [6]       keplerian elements array: kep = [a e i OM om theta][km, rad]
%  mu [1]        gravitational constant                             [km^3/s^2]
% 
% OUTPUT:
%  rr [3,1]      position vector                                    [km]
%  vv [3,1]      velocity vector                                    [km/s]
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

p = kep(1)*(1-kep(2)^2);
r = p/(1+kep(2)*cos(kep(6)));

rr = r*[cos(kep(6));sin(kep(6));0];
vv = sqrt(mu/p)*[-sin(kep(6));kep(2)+cos(kep(6));0];

ROM = [ cos(kep(4)) sin(kep(4)) 0; -sin(kep(4)) cos(kep(4)) 0; 0 0 1];
Ri = [ 1 0 0; 0 cos(kep(3)) sin(kep(3)); 0 -sin(kep(3)) cos(kep(3))];
Rom = [ cos(kep(5)) sin(kep(5)) 0; -sin(kep(5)) cos(kep(5)) 0; 0 0 1];

T = Rom*Ri*ROM;

rr = T'*rr;
vv = T'*vv;
