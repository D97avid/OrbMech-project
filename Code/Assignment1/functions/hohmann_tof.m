function [dt] = hohmann_tof(a1,a2,mu)
%
% hohmann_tof time of flight of a Hohmann transfer
% 
% Function to compute the time of flight of a Hohmann transfer between two
% circular and coplanar orbits.
% 
% PROTOTYPE:
%  dt = hohmann_tof(a1,a2,mu)
% 
% INPUT:
%  a1 [1]        semi-major axis of departure orbit                 [km]
%  a2 [1]        semi-major axis of arrival orbit                   [km]
%  mu [a]        planetary constants of the centre of attraction    [km^3/s^2]
% 
% OUTPUT:
%  dt [1]        time of flight                                     [s]
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

at = (a1+a2)/2;             % semi-major axis of the transfer orbit

dt = pi*sqrt(at^3/mu);      % time of flight computation (half of the period)


