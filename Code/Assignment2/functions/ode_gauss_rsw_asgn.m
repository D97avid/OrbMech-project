function dkep = ode_gauss_rsw_asgn (t,kep,muE,aj2_RSW,a_moon_RSW)
% 
% Function to compute Planetary Gauss equations, due to Moon's gravity and J2 perturbation, 
% in the RSW (Radial-Transversal_Out of plane) frame. This function is used
% into an ODE SOLVER (ode113,ode45,..).
% 
% PROTOTYPE:
%  dkep = ode_gauss_rsw_asgn (t,kep,muE,aj2_RSW,a_moon_RSW)
% 
% INPUT:
% t [1]                      time for ephemerides of Moon and for ODE Solver [s]  
% kep [1x6]                  kepler parameters (a,e,i,OM,om,th)
% muE [1]                    Earth's gravitational parameter [km^3/s^2]
% aj2_RSW [1x3]              Perturbed acceleration (J2) in RSW frame [km/s^2]
% a_moon_RSW [1x3]           Perturbed acceleration (Moon) in RSW frame [km/s^2]
% 
% OUTPUT:
% dkep [6x1]                  set of derivatives of kepler elements: [da,de,di,dOM,dom,dth].
% 
% CONTRIBUTORS:
%  Marco Adorno
%  Giuseppe Esposito 
%  Davide Gravina 
%  David Reina
% 
% VERSIONS:
%  20-01-2021: First version


%SPACECRAFT COMPUTATION:
% Set of kepler parameters:
as=kep(1);      % semi-major axis ;
es=kep(2);      % eccentricity ;
is=kep(3);      % inclination ;
OMs=kep(4);     % right ascension of the ascending node ;
oms=kep(5);     % arguument of periapsis ;
ths=kep(6);     % true anomaly ;

% Computation of usefull realations :
ps=as*(1-es^2);                     % semi-lactus rectum
rs=ps/(1+es*cos(ths));              % distance of position (modulus radius)
bs=as*sqrt(1-es^2);                 % semi-minor axis
ns=sqrt(muE/as^3);                  % mean motion
hs=ns*as*bs;                        % angular momentum

norm_r=rs;

%Component of AJ2 in RSW frame:
aj2_R= aj2_RSW(1);
aj2_S= aj2_RSW(2);
aj2_W= aj2_RSW(3);

%Component of A_MOON in RSW frame:
a_moon_R= a_moon_RSW(1);
a_moon_S= a_moon_RSW(2);
a_moon_W= a_moon_RSW(3);

%Component of A_PERTURBED in RSW frame:
a_pert_R= aj2_R + a_moon_R;
a_pert_S= aj2_S + a_moon_S;
a_pert_W= aj2_W + a_moon_W;

%GAUSS PLANETARY EQUATION in RSW frame (set of derivates):

da = (2*(as^2)/hs) * (es*sin(ths)*a_pert_R + (ps/norm_r)*a_pert_S);
de = (1/hs) * ( ps*sin(ths)*a_pert_R + ( (ps+norm_r)*cos(ths) + norm_r*es )*a_pert_S );
di = ( norm_r*cos(ths+oms)/hs ) * a_pert_W;
dOM = ( norm_r*sin(ths+oms)/(hs*sin(is)) ) * a_pert_W;
dom = (1/(hs*es)) * ( -ps*cos(ths)*a_pert_R + (ps+norm_r)*sin(ths)*a_pert_S ) - ( (norm_r*sin(ths+oms)*cos(is)) / (hs*sin(is)) )*a_pert_W ;
dth =  (hs/(norm_r^2)) + (1/(es*hs)) * ( ps*cos(ths)*a_pert_R - (ps+norm_r)*sin(ths)*a_pert_S );

dkep=[da;de;di;dOM;dom;dth];
