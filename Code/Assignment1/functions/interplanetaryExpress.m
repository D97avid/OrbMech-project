function [dv_tot,dv_dep,dv_flyby,dv_arr,rp_flyby,vp_m,vp_p,inc_flyby,vv1_t1,vv1_t2,vv2_t1,vv2_t2,rrE_flyby,vvE_flyby,kep1_t1,kep2_t1,kep1_t2,kep2_t2] = interplanetaryExpress(t1,t2,t3)
% 
% Function to compute entire interplanetary transfer cost and geometry.
% The transfer computed is: Jupiter - Earth (powered gravity assist) - Mercury.
% However this configuration can be modified if the dep_planet, arr_planet
% and flyby_planet variables are changed accordingly.
% 
% PROTOTYPE:
% [dv_tot,dv_dep,dv_flyby,dv_arr,rp_flyby,vp_m,vp_p,inc_flyby,vv1_t1,vv1_t2,vv2_t1,vv2_t2,rrE_flyby,kep1_t1,kep2_t1,kep1_t2,kep2_t2] = interplanetaryExpress(t1,t2,t3)
% 
% INPUT:
%  t1 [1]          departure time in MJD2000                           [days]
%  t2 [1]          flyby time in MJD2000                               [days]
%  t3 [1]          arrival time in MJD2000                             [days]
% 
% OUTPUT:
%  dv_tot [1]       total delta v cost of the interplanetary transfer  [km/s]
%  dv_dep [1]       delta v cost of the departure manoeuvre            [km/s]
%  dv_flyby [1]     delta v cost of the powered gravity assist         [km/s]
%  dv_arr [1]       delta v cost of the arrival manoeuvre              [km/s]
%  rp_flyby [1]     modulus of the radius of perigee of the flyby      [km]
%  vp_m [1]         incoming velocity magnitude at perigee             [km/s]
%  vp_p [1]         outgoing velocity magnitude at perigee             [km/s]
%  inc_flyby [1]    flyby hyperbola plane inclination                  [rad]
%  vv1_t1 [3,1]     starting velocity on first Lambert arc             [km/s]
%  vv1_t2 [3,1]     final velocity on first Lambert arc                [km/s]
%  vv2_t1 [3,1]     starting velocity on second Lambert arc            [km/s]
%  vv2_t2 [3,1]     final velocity on second Lambert arc               [km/s]
%  rrE_flyby [3,1]  Earth flyby position vector in heliocentric frame  [km]
%  vvE_flyby [3,1]  Earth flyby velociy vector in heliocentric frame   [km]
%  kep1_t1[1,6]     keplerian elements array of the first transfer leg
%                   initial point                                      [km rad]
%  kep2_t1[1,6]     keplerian elements array of the first transfer leg
%                   final point                                        [km rad]
%  kep1_t2[1,6]     keplerian elements array of the second transfer leg
%                   initial point                                      [km rad]
%  kep2_t2[1,6]     keplerian elements array of the second transfer leg
%                   final point                                        [km rad]
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

dv_tot = NaN;
dv_dep = NaN;
dv_flyby = NaN;
dv_arr = NaN;
rp_flyby = NaN;
inc_flyby = NaN;
vv1_t1 = NaN;
vv1_t2 = NaN;
vv2_t1 = NaN;
vv2_t2 = NaN;
vp_m = NaN;
vp_p = NaN;
kep1_t1 = NaN;
kep2_t1 = NaN;
kep1_t2 = NaN;
kep2_t2 = NaN;


if t2-t1 >= 0 && t3-t2 >= 0
    dep_planet = 5;         % Jupiter
    flyby_planet = 3;       % Earth
    arr_planet = 1;         % Mercury
    muS = astroConstants(4);

    % departure planet
    kep1 = uplanet(t1, dep_planet);
    [rr1,vv1] = kep2car(kep1,muS);

    % fly-by planet
    kep2 = uplanet(t2, flyby_planet);
    [rrE_flyby,vvE_flyby] = kep2car(kep2,muS);

    % arrival planet
    kep3 = uplanet(t3, arr_planet);
    [rr3,vv3] = kep2car(kep3,muS);

    % transfer leg 1
    [~,~,~,ERROR1,vv1_t1,vv2_t1,~,~] = lambertMR( rr1, rrE_flyby, (t2-t1)*3600*24, muS, 0, 0, 0 );
    vv1_t1 = vv1_t1';
    vv2_t1 = vv2_t1';
    kep1_t1 = car2kep(rr1,vv1_t1,muS);
    kep2_t1 = car2kep(rrE_flyby,vv2_t1,muS);
    % transfer leg 2
    [~,~,~,ERROR2,vv1_t2,vv2_t2,~,~] = lambertMR( rrE_flyby, rr3, (t3-t2)*3600*24, muS, 0, 0, 0 );
    vv1_t2 = vv1_t2';
    vv2_t2 = vv2_t2';
    kep1_t2 = car2kep(rrE_flyby,vv1_t2,muS);
    kep2_t2 = car2kep(rr3,vv2_t2,muS);
    
    if ERROR1 == 0 && ERROR2 == 0
        v_inf_m = vv2_t1 - vvE_flyby;
        v_inf_p = vv1_t2 - vvE_flyby;
        [dv_flyby,rp_flyby,vp_m,vp_p,inc_flyby] = poweredGravityAssist(v_inf_m,v_inf_p);
        if not(isnan(rp_flyby))
            dv_dep = norm(vv1_t1 - vv1);
            dv_arr = norm(vv2_t2 - vv3);
            dv_tot = dv_dep + dv_flyby + dv_arr;
        end
    end
end

