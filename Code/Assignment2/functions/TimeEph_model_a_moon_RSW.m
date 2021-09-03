function a_moon_RSW = TimeEph_model_a_moon_RSW (t,kep,muMoon,muE)
% 
% Function to compute the pertubation accelaration, due to Moon's gravity, 
% in the RSW (Radial-Transversal_Out of plane) frame. This function is used
% into an ODE SOLVER (ode113,ode45,..).
% 
% PROTOTYPE:
% a_moon_RSW = TimeEph_model_a_moon_RSW (t,kep,muMoon,muE)
%
% INPUT:
% t                                            time array [s]
% kep [1x6]                                    Kepler parameters (a,e,i,OM,om,th) [km,-,rad,rad,rad,rad]
% muMoon [1]                                   Moon's gravitational parameter [km^3/s^2]
% mu_E [1]                                     Earth's gravitational parameter [km^3/s^2] 
%
% OUTPUT:
% a_moon_RSW                                   Perturbed acceleration (Moon) in RSW frame [km/s^2]
% 
% CONTRIBUTORS:
%  Marco Adorno
%  Giuseppe Esposito 
%  Davide Gravina 
%  David Reina
% 
% VERSIONS:
%  20-01-2021: First version

%Position vector of spacecraft, as function of Keplerian elements:
kep_sc=kep;
[rr_sc,~] = kep2car(kep_sc,muE);
              
% Set of Cartesian position of MOON (XYZ) from EPHEMERIDES OF MOON:
t_MJD2000=t/86400; %[day] time from second to day of MJD2000 time
[rm, ~] = ephMoon(t_MJD2000); %ATTENTION !!! DEFINE 't' as MJD2000 time.
nrm=norm(rm); %[km] modulus of Moon's position wrt the Earth
%where: rm,vm : are respectively position and velocity vector of Moon wrt the Earth.

% Set of Cartesian position of MOON wrt the S/C (XYZ):
rms=rm'-rr_sc; %[km] Position vector of Moon wrt the s/c
nrms=norm(rms); %[km] modulus of Moon's position wrt the s/c

%Perturbation acceleration due to Moon gravity:
a_moon_XYZ = muMoon * ( (rms/(nrms^3) ) - (rm'/(nrm^3)) );

%Keplerian elements OF SPACECRAFT:
    a_sc = kep_sc(1);
    e_sc = kep_sc(2);
    i_sc = kep_sc(3);
    OM_sc = kep_sc(4);
    om_sc = kep_sc(5);
    th_sc = kep_sc(6);

%Rotation matrix (direction cosines) from XYZ (cartesian) frame to RSW frame, as function of Keplerian elements of SPACECRAFT:
u_sc=om_sc+th_sc; %[rad] argoument of latitude OF SPACECRAFT;
ROT_xyz_rsw=[-sin(OM_sc)*cos(i_sc)*sin(u_sc)+cos(OM_sc)*cos(u_sc), cos(OM_sc)*cos(i_sc)*sin(u_sc)+sin(OM_sc)*cos(u_sc), sin(i_sc)*sin(u_sc);
             -sin(OM_sc)*cos(i_sc)*cos(u_sc)-cos(OM_sc)*sin(u_sc), cos(OM_sc)*cos(i_sc)*cos(u_sc)-sin(OM_sc)*sin(u_sc), sin(i_sc)*cos(u_sc);
                          sin(OM_sc)*sin(i_sc)                   ,              -cos(OM_sc)*sin(i_sc)                 ,      cos(i_sc)    ];
%Computation of aj2_XYZ in RSW frame:
a_moon_RSW=ROT_xyz_rsw*a_moon_XYZ;
                

end