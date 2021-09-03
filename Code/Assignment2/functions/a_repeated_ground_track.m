function a = a_repeated_ground_track(k,m,mu,omega_E)
% 
% Function that computes the semi-major axis required to obtain a repeating ground track with the given values of k,m
% 
% PROTOTYPE:
%  a = a_repeated_ground_track(k,m,mu,omega_E)
%
% INPUT:
%  k [1]                  number of revolutions of the s/c to obtain the ground track repetition
%  m [1]                  number of revolutions of the planet to obtain the ground track repetition
%  mu [1]                 gravitational parameter of primary body [km^3/s^2]
%  omega_E [1]            Earth's spin rate [rad/s]
%
% OUTPUT: 
%  a [1]                  semi-major axis for the ground track repetition [km]
%
% CONTRIBUTORS:
%  Marco Adorno
%  Giuseppe Esposito 
%  Davide Gravina 
%  David Reina
% 
% VERSIONS:
%  20-01-2021: First version

a = nthroot(((mu*m^2)/(omega_E*k)^2),3); %computes the semi-major axis 

