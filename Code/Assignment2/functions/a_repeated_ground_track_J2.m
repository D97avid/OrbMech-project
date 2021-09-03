function a = a_repeated_ground_track_J2(k,m,mu,omega_E,e,J2,R_e,i,a_repeated)
% 
% Function that computes the semi-major axis required to obtain a repeating
% ground track with the given values of k,m, taking into account the
% secular effect of J2
% 
% PROTOTYPE:
%  a = a_repeated_ground_track(k,m,mu,omega_E)
%
% INPUT:
%  k [1]                  number of revolutions of the s/c to obtain the ground track repetition
%  m [1]                  number of revolutions of the planet to obtain the ground track repetition
%  mu [1]                 gravitational parameter of primary body [km^3/s^2]
%  omega_E [1]            Earth's spin rate [rad/s]
%  e [1]                  s/c orbit's eccentricity
%  i [1]                  s/c orbit's inclination [rad]
%  J2 [1]                 second zonal harmonic
%  R_e [1]                Earth's mean radius [km]
%  a_repeated [1]         semi-major axis for the repetition of the unperturbed ground track [km]
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

c = (-3/2)*((sqrt(mu)*J2*(R_e^2))/((1-e^2)^2)); %constant coefficient

% a_repeated is used as initial guess for fzero
a = fzero(@(a) (m/k) - ((omega_E - (c/a^(7/2))*cos(i))/((sqrt(mu/a^3))+(c/a^(7/2))*((5/2)*(sin(i))^2-2)-((c*sqrt(1-e^2))/a^(7/2))*(1-(3/2)*(sin(i))^2))),a_repeated);
