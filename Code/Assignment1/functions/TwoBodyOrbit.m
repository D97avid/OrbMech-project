function dy = TwoBodyOrbit(~,y,mu)
% 
% ODE system for the two body orbit propagation
% 
% PROTOTYPE:
%  dy = TwoBodyOrbit(t,y,mu)
% 
% INPUT:
%  t [1]     time (can be omitted because, as the system is autonomus)[s]
%  y [6,1]   state of the system y = [position; velocity]             [km,km/s]
%  mu [1]    gravitational constant                                   [km^3/s^2]
% 
% OUTPUT:
%  dy [1]    derivative of the state                                  [km/s^2 km/s^3]
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

rr = y(1:3);
vv = y(4:6);

r = norm(rr);

dy = [ vv(1);vv(2);vv(3);
       -mu/r^3 * rr(1); -mu/r^3 * rr(2); -mu/r^3 * rr(3) ];

end