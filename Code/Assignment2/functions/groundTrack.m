function [alpha_deg,delta_deg,lon_wrapped,lat] = groundTrack(input_type,orbit_initial_state,thetaG_0,time_span,omega_E,mu,t0)
% 
% Function to plot the Ground Track of an orbit during a given time period.
% 
% PROTOTYPE:
% [alpha_deg,delta_deg,lon_wrapped,lat] = groundTrack(input_type,orbit_initial_state,thetaG_0,time_span,omega_E,mu,t0)
%
% INPUT:
%  input_type [char]                     describes how the orbit initial state is given ('Keplerian' = Keplerian Elements; 'Cartesian' = Cartesian Coordinates)
%  orbit_initial_state [6,1]             orbit initial state vector which with Keplerian Elements [km],[rad] in case of Keplerian or initial position and velocity in case of 'Cartesian' [km],[km/s]
%  thetaG_0 [1]                          right ascension of Greenwich meridian at t0    [rad]
%  time_span [1,N]                       time span vector (N is the number of points used for the linspace) [s]
%  omega_E [1]                           Earth's spin rate [deg/h]
%  mu [1]                                gravitational parameter of primary body [km^3/s^2]
%  t0 [1]                                initial time [s]
%
% OUTPUT: 
%  alpha_deg [N,1]                       vector of right-ascension value at each time-step [deg]
%  delta_deg [N,1]                       vector of declination value at each time-step [deg]
%  lon_wrapped [N,1]                     vector of longitude value at each time-step [deg]
%  lat [N,1]                             vector of latitude value at each time-step [deg]
%
% CONTRIBUTORS:
%  Marco Adorno
%  Giuseppe Esposito 
%  Davide Gravina 
%  David Reina
% 
% VERSIONS:
%  20-01-2021: First version

%% Orbit Initial State Analysis

switch input_type %analize the type of input, and eventually converts it into a Cartesian initial state
    case 'Keplerian'
        initial_state = ones(6,1);
        [initial_state(1:3),initial_state(4:6)] = kep2car(orbit_initial_state,mu);
    case 'Cartesian'
        initial_state = orbit_initial_state;
    otherwise
        error('Input type not valid');
end

%% Orbit Propagation

options = odeset('RelTol',1e-13,'AbsTol',1e-14);
[~ , x] = ode113(@(t,x) kepl_orbit(x,mu),time_span,initial_state,options); %propagate the orbit during the time-span

%% Conversion to angles for Ground Track

delta = NaN(length(time_span),1); %initialization of delta vector
for i=1:length(time_span)
    delta(i) = asin(x(i,3)/norm(x(i,1:3))); %computes the declination
end
delta_deg = rad2deg(delta); %converts the declination in degrees

alpha = NaN(length(time_span),1); %initialization of alpha vector
for i=1:length(time_span)
    if (x(i,2)/norm(x(i,1:3)))>0
        alpha(i) = acos((x(i,1)/norm(x(i,1:3)))/cos(delta(i))); %computes the right ascension
    else
        alpha(i) = (2*pi) - acos((x(i,1)/norm(x(i,1:3)))/cos(delta(i)));
    end
end
alpha_deg = rad2deg(alpha); %converts the right ascension in degrees

lat = delta_deg; %computes the latitude vector in degrees

thetaG = ones(length(time_span),1)*thetaG_0 +( deg2rad(omega_E/3600)) * (time_span' - t0*ones(length(time_span),1)); %computes the right ascension of Greenwich at each time-step
thetaG_deg = rad2deg(thetaG); %converts the thetaG vector in degrees
lon = alpha_deg - thetaG_deg; %computes the longitude vector in degrees
lon_wrapped = wrapTo180(lon); %wraps the longitude vector between [-180 180]

for i=2:length(lon_wrapped)
    if abs(lon_wrapped(i)-lon_wrapped(i-1))>180 %solving of the ambiguity when crossing the anti-meridian for plot purposes
        lon_wrapped(i) = NaN;
    end
end
