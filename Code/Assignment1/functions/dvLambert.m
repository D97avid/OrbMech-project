function [dv_tot,dv1,dv2] = dvLambert(t1,t2,dep_planet,arr_planet)
% 
% Function to compute the delta v cost of a Lambert transfer between two
% planets of the Solar System.
% 
% PROTOTYPE:
%  [dv_tot,dv1,dv2] = dvLambert(t1,t2,dep_planet,arr_planet)
% 
% INPUT:
%  t1 [1]        departure time from starting orbit in MJD2000      [day]
%  t2 [1]        arrival time to final orbit  in MJD2000            [day]
%  dep_planet[1] departure planet ID (for LambertMR function)       [-]
%  arr_planet[1] arrival planet ID (for LambertMR function)         [-]
% 
% OUTPUT:
%  dv_tot [1]    total delta v cost of the Lambert transfer         [km/s]
%  dv1 [1]       delta v cost of the first manoeuvre                [km/s]
%  dv2 [1]       delta v cost of the second manoeuvre               [km/s]
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

% Initialization of the values in case the constraints aren't respected:
dv_tot = NaN;
dv1 = NaN;
dv2 = NaN;

% First constraint: negative time transfer window

if t2 - t1 >= 0
    muS = astroConstants(4);    % Sun gravitational constant [km^3/s^2]

    % departure planet 
    kep1 = uplanet(t1, dep_planet); % keplerian elements
    [rr1,vv1] = kep2car(kep1,muS);  % position and velocity

    % arrival planet
    kep2 = uplanet(t2, arr_planet); % keplerian elements
    [rr2,vv2] = kep2car(kep2,muS);  % position and velocity

    % conversion of time in seconds
    deltat =(t2-t1)*3600*24;    

    % transfer computation
    [~,~,~,ERROR,vv1_t,vv2_t,~,~] = lambertMR( rr1, rr2, deltat, muS, 0, 0, 0 );
    
    % delta v computation
    if ERROR == 0
        dv1 = norm(vv1_t' - vv1);
        dv2 = norm(vv2 - vv2_t');
        dv_tot = dv1 + dv2;
    end
end
