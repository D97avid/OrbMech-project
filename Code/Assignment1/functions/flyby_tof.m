function dt_flyby = flyby_tof(a_hyp_m,e_hyp_m,a_hyp_p,e_hyp_p,rrE_flyby)
% 
% Function to compute the duration of the flyby around Earth considering
% a finite sphere of influence
% 
% PROTOTYPE:
% dt_flyby = flyby_tof(a_hyp_m,e_hyp_m,a_hyp_p,e_hyp_p,rrE_flyby)
%
% INPUT:
%  a_hyp_m [1]          semi-major axis of incoming hyperbola           [km]
%  e_hyp_m [1]          eccentricity of incoming hyperbola              [-]
%  a_hyp_p [1]          semi-major axis of outgoing hyperbola           [km]
%  e_hyp_p [1]          eccentricity of outgoing hyperbola              [-]
%  rrE_flyby [3,1]      Earth position vector at flyby                  [km]
% 
% OUTPUT:
%  dt_flyby [1]        duration of the flyby                            [s]
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

% Constants:
muS = astroConstants(4);                   % Sun gravitational constant
muE = astroConstants(13);                  % Earth gravitational constant

% Radius of Earth's SOI
r_SOI = norm(rrE_flyby)*(muE/muS)^(2/5); 

% Apply hyperbola time law:
% Incoming hyperbola:
th_m = acos((a_hyp_m*(1-e_hyp_m^2)-r_SOI)/e_hyp_m/r_SOI);
F_m = 2*atanh(tan(th_m/2)/sqrt((1+e_hyp_m)/(e_hyp_m-1)));
dt_m = sqrt((abs(a_hyp_m))^3/muE) * (e_hyp_m*sinh(F_m)-F_m);

% Outgoing hyperbola:
th_p = acos((a_hyp_p*(1-e_hyp_p^2)-r_SOI)/e_hyp_p/r_SOI);
F_p = 2*atanh(tan(th_p/2)/sqrt((1+e_hyp_p)/(e_hyp_p-1)));
dt_p = sqrt((abs(a_hyp_p))^3/muE) * (e_hyp_p*sinh(F_p)-F_p);

dt_flyby = dt_m + dt_p;