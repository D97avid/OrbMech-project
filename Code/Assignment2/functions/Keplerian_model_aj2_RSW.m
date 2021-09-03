function aj2_RSW = Keplerian_model_aj2_RSW (kep,muE,RE,J2)
% 
% Function to compute the pertubation accelaration, due to J2 effect, in
% the RSW (Radial-Transversal_Out of plane) frame. This function is used
% into an ODE SOLVER (ode113,ode45,..).
%
% PROTOTYPE:
%  aj2_RSW = Keplerian_model_aj2_RSW (kep,muE,RE,J2)
% 
% INPUT:
% kep [1x6]             kepler parameters (a,e,i,OM,om,th) [km,-,rad,rad,rad,rad]
% muE [1]               Earth's gravitational parameter [km^3/s^2]
% RE  [1]               Earth's radius [km]
% J2  [1]               Coefficient for the second zonal harmonic [-] for Earth
% 
% OUTPUT:
% aj2_RSW [1x3]         Perturbed acceleration in RSW frame
% 
% CONTRIBUTORS:
%  Marco Adorno
%  Giuseppe Esposito 
%  Davide Gravina 
%  David Reina
% 
% VERSIONS:
%  20-01-2021: First version

% Set of kepler parameters:
a=kep(1);   % semi-major axis ;
e=kep(2);   % eccentricity ;
i=kep(3);   % inclination ;
OM=kep(4);  % right ascension of the ascending node ;
om=kep(5);  % arguument of periapsis ;
th=kep(6);  % true anomaly ;

% Set cartesian initial data:
[rr,~] = kep2car([a,e,i,OM,om,th],muE); 

% Compute useful parameter:
norm_r=norm(rr);    % modulus of position vector
u=th+om;            % argoument of latitude


% Formula of Perturbed acceleration (aj2), due to J2 (oblatness of Earth), as
% function of Keplerian elements, in RSW frame (Radial - Transversal - Out of plane):
aj2_coeff = -(3/2)*(J2*muE*(RE^2))/(norm_r^4);

aj2_matrix = [ 1-(3*(sin(i)^2)*(sin(u)^2));...
             (sin(i)^2)*sin(2*u);...
             sin(2*i)*sin(u)];
               
aj2_RSW = aj2_coeff * aj2_matrix; 
