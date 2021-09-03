function [dv,rp,vp_m,vp_p,v_inf_m_n,v_inf_p_n,a_hyp_m,a_hyp_p,e_hyp_m,e_hyp_p] = poweredGravityAssist(v_inf_m,v_inf_p)
% 
% Function to compute the powered gravity assist manoeuvre around Earth.
% 
% PROTOTYPE:
%  [dv,rp,vp_m,vp_p,v_inf_m_n,v_inf_p_n,a_hyp_m,a_hyp_p,e_hyp_m,e_hyp_p] = poweredGravityAssist(v_inf_m,v_inf_p)
% 
% INPUT:
%  v_inf_m [3,1]    incoming excess velocity vector                    [km/s]
%  v_inf_p [3,1]    outgoing excess velocity vector                    [km/s]
% 
% OUTPUT:
%  dv [1]           delta v cost of the powered gravity assist         [km/s]
%  rp [1]           radius of perigee of the flyby trajectory          [km]
%  vp_m [1]         incoming velocity magnitude at perigee             [km/s]
%  vp_p [1]         outgoing velocity magnitude at perigee             [km/s]
%  v_inf_m_n [1]    magnitude of incoming excess velocity              [km/s]
%  v_inf_p_n [1]    magnitude of outgoing excess velocity              [km/s]
%  a_hyp_m [1]      incoming hyperbola semi-major axis                 [km]
%  a_hyp_p [1]      outgoing hyperbola semi-major axis                 [km]
%  e_hyp_m [1]      incoming hyperbola eccentricity                    [-]
%  e_hyp_p [1]      outgoing hyperbola eccentricity                    [-]
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

% Constants
muE = astroConstants(13);

% Initial computations
turn_angle = acos((dot(v_inf_p',v_inf_m'))/(norm(v_inf_m)*norm(v_inf_p)));
v_inf_m_n = norm(v_inf_m);
v_inf_p_n = norm(v_inf_p);

% fsolve
rp_guess = 6378;
options = optimset('Display','off');
rp = fzero( @(rp) -turn_angle + asin(1/(1 + rp/muE*norm(v_inf_m)^2)) + asin(1/(1 + rp/muE*norm(v_inf_p)^2)),rp_guess,options);

rp_critical = 6378.137 + 200;   % critical radius of perigee (Earth radius plus atmosphere height)

dv = NaN;
vp_m = NaN;
vp_p = NaN;
a_hyp_m = NaN;
a_hyp_p = NaN;
e_hyp_m = NaN;
e_hyp_p = NaN;

if rp > rp_critical
    a_hyp_m = -muE/norm(v_inf_m)^2;
    a_hyp_p = -muE/norm(v_inf_p)^2;
    %turn_angle_m = asin(1/(1 + rp/muE*norm(v_inf_m)^2));
    %turn_angle_p = asin(1/(1 + rp/muE*norm(v_inf_p)^2));
    e_hyp_m = 1+rp*norm(v_inf_m)^2/muE;
    e_hyp_p = 1+rp*norm(v_inf_p)^2/muE;
    vp_m = sqrt(2*muE*(1./rp + 1/2/abs(a_hyp_m)));
    vp_p = sqrt(2*muE*(1./rp + 1/2/abs(a_hyp_p)));
    dv = abs(vp_p - vp_m);
else
    rp = NaN;
end


