function kep = car2kep(rr,vv,mu)
% 
% Function to compute the keplerian elements vector from the position and
% velocity vector and the gravitational constant of the centre of attraction.
% 
% PROTOTYPE:
%  kep = car2kep(rr,vv,mu)
% 
% INPUT:
%  rr [3,1]      position vector                                    [km]
%  vv [3,1]      velocity vector                                    [km/s]
%  mu [1]        gravitational constant                             [km^3/s^2]
% 
% OUTPUT:
%  kep [6]       keplerian elements array: kep = [a e i OM om theta][km, rad]
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

i = [1 0 0];
j = [0 1 0];
k = [0 0 1];

r = norm(rr);
v = norm(vv);

a = 1/(2/r - (v^2)/mu);

hh = cross(rr,vv);
h = norm(hh);

ee = cross(vv,hh)/mu - rr/r;
e = norm(ee);

i = acos(hh(3)/h);

N = cross(k,hh);

if norm(N) ~= 0
    N = N/norm(N);
end

if N(2) >= 0 
    OM = acos(N(1));
else
    OM = 2*pi - acos(N(1));
end

if ee(3) >= 0
    om = acos(dot(N,ee)/e);
else
    om = 2*pi - acos(dot(N,ee)/e);
end

vr = dot(vv,rr)/r;

if vr > 0 
    th = acos(dot(rr,ee)/r/e);
elseif vr < 0
    th = 2*pi - acos(dot(rr,ee)/r/e);
else
    th = 0;
end

kep = [a e i OM om th];


