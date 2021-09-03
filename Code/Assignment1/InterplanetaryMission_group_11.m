%% Main script for Interplanetary Mission assignment %%
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
%% Initialization:
clear,clc,close all;

% Path adjustments:
addpath(strcat(pwd,'/functions'));
addpath(strcat(pwd,'/functions/textures'));
addpath('../Provided-scripts');
addpath('../Provided-scripts/time');

% Loading of constants and input data:
muS = astroConstants(4);                % Sun gravitational constant
muE = astroConstants(13);               % Earth gravitational constant

% Departure, flyby and arrival planets identification numbers:
Departure_planet = 5;                  
Flyby_planet = 3;
Arrival_planet = 1;

% Conversion of the earliest departure and latest arrival dates 
% in Modified Julian day 2000:
Earliest_departure = date2mjd2000([2030 6 1 0 0 0]);
Latest_arrival = date2mjd2000([2070 6 1 0 0 0]);

%% Preliminary estimations
% Orbital periods computation
kep_dep = uplanet(0, Departure_planet);             % Keplerian elements vector of departure planet  
T_dep = 1/(86400) * 2*pi*sqrt(kep_dep(1)^3/muS);    % Orbital period of departure planet [days]

kep_FB  = uplanet(0, Flyby_planet);                 % Keplerian elements vector of flyby planet
T_FB = 1/(86400)* 2*pi*sqrt(kep_FB(1)^3/muS);       % Orbital period of flyby planet [days]

kep_arr = uplanet(0, Arrival_planet);               % Keplerian elements vector of arrival planet
T_arr = 1/(86400) * 2*pi*sqrt(kep_arr(1)^3/muS);    % Orbital period of arrival planet [days]

% Synodic periods computation
T_syn_DepwrtFB = T_dep*T_FB / abs(T_dep-T_FB);      % Synodic period of departure planet w.r.t flyby planet
T_syn_FBwrtArr = T_FB*T_arr / abs(T_FB-T_arr);      % Synodic period of flyby planet w.r.t arrival planet
T_syn_DepwrtArr = T_dep*T_arr / abs(T_dep-T_arr);   % Synodic period of departure planet w.r.t arrival planet

% Narrowing of the departure, flyby and arrival windows
t1 = [];                                            % departure window
t2 = [];                                            % flyby window
t3 = [];                                            % arrival window

% Due to the short synodic periods, the departure window has been
% restricted to a single orbital period of the departure planet (Jupiter).

t1 = [Earliest_departure Earliest_departure+T_dep];

% Upper and lower constraints have been imposed to the flyby and arrival 
% windows, considering 70% and 130% of a generic Hohmann transfer time of flight 
% between the planets assuming circular and coplanar orbits:

A = hohmann_tof(kep_dep(1),kep_FB(1),muS)/3600/24;   % Hohmann Jupiter-Earth
B = hohmann_tof(kep_FB(1),kep_arr(1),muS)/3600/24;   % Hohmann Earth - Mercury

t2 = [t1(1)+0.7*A, t1(end)+1.3*A];        
t3 = [t2(1)+0.7*B, t2(end)+1.3*B];   
%%
fprintf('----------------------------------------------\n');
fprintf('     INTERPLANETARY MISSION ASSIGNMENT        \n');
fprintf('----------------------------------------------\n');
%% Genetic Algorithm (ga)
rng default                               % For riproducibility
N_runs = 5;                               % Number of runs
lb_ga = [t1(1) t2(1) t3(1)];              % Lower boundary   
ub_ga = [t1(end) t2(end) t3(end)];        % Upper boundary

% Initialization of output array
dv_ga_vect = ones(N_runs,1);              % minimum dv array obtained via Genetic Algorithm
x_ga = zeros(N_runs,3);                   % minimum dv time array obtained via Genetic Algorithm
                               
% Options definition for the Genetic Algorithm:
options = optimoptions('ga','MaxGeneration',100,'PopulationSize',1000,...
    'FunctionTolerance',0.01,...
    'PlotFcn',"gaplotbestf",'Display','off');

fprintf('Genetic Algorithm running...\n');
fprintf('Progress:   0%%');
tic
for i = 1:N_runs
    [x_ga(i,:),dv_ga_vect(i)] = ga(@(x) interplanetaryExpress(x(1),x(2),x(3)),3,[],[],[],[],lb_ga,ub_ga,[],options);
    if round(i/N_runs*100) < 10
        fprintf('\b\b');
        fprintf('%d%%',round(i/N_runs*100));
    else
        fprintf('\b\b\b');
        fprintf('%d%%',round(i/N_runs*100));
    end
end
fprintf('\nFinished. \n');
toc

dt_ga_vect = x_ga(:,3)-x_ga(:,1);               % dt array obtained from all the runs of Genetic Algorithm
[dv_ga,ii] = min(dv_ga_vect);
t_ga = x_ga(ii,:);                              % time array of the minimum dv obtained via Genetic Algorithm
dt_ga = t_ga(3) - t_ga(1);
fprintf('Results:\n');
fprintf('dv_ga = %.4f km/s\n', dv_ga)
fprintf('dt_ga = %.4f days\n', dt_ga)
fprintf('----------------------------------------------\n');

%% Refining the solution obtained using other algorithms: %%
% Definition of a shorter time window in the neighborhood of the output
% solution of the Genetic Algorithm:
delta1 = T_FB/2;
delta2 = T_FB*delta1/T_dep;
delta3 = T_arr*delta1/T_dep;

% Refined departure, flyby and arrival time windows:
t1_new = linspace(t_ga(1)-delta1,t_ga(1)+delta1,2);
t2_new = linspace(t_ga(2)-delta2,t_ga(2)+delta2,2);
t3_new = linspace(t_ga(3)-delta3,t_ga(3)+delta3,2);

%% Gradient based optimization (fmincon)
lb_fmincon = [t1_new(1) t2_new(1) t3_new(1)];           % lower boundary
ub_fmincon = [t1_new(end) t2_new(end) t3_new(end)];     % upper boundary
x0 = t_ga;                                              % starting point

options = optimoptions('fmincon','Algorithm','sqp','PlotFcn',"optimplotfval",...
    'FunctionTolerance',0.01,'StepTolerance',1e-10,'ConstraintTolerance',1e-10,...
    'display','off');

fprintf('Gradient-based optimization running...\n');

% Initialization of output
dv_fmincon = [];                           % minimum dv obtained via Gradient based optimization
t_fmincon = [];                            % minimum dv time array obtained via Gradient based optimization

tic
[t_fmincon,dv_fmincon] = fmincon(@(x) interplanetaryExpress(x(1),x(2),x(3)),x0,[],[],[],[],lb_fmincon,ub_fmincon,[],options);
fprintf('Finished. \n')
toc

dt_fmincon = t_fmincon(3) - t_fmincon(1);

fprintf('Results:\n');
fprintf('dv_fmincon = %.4f km/s\n', dv_fmincon)
fprintf('dt_fmincon = %.4f days\n', dt_fmincon)
fprintf('----------------------------------------------\n');
%% Grid search
% Definition of the new time windows for the grid search algorithm
t1_gs = linspace(t1_new(1),t1_new(end),100);        
t2_gs = linspace(t2_new(1),t2_new(end),100);
t3_gs = linspace(t3_new(1),t3_new(end),100);

% Initialization of output array
dv_gs_mat = zeros(length(t1_gs),length(t2_gs),length(t3_gs));

fprintf('Grid search algorithm running...\n');
fprintf('Progress:   0%%');
tic
for i = 1:length(t1_gs)
    for j = 1:length(t2_gs)
        for k = 1:length(t3_gs)
            [dv_gs_mat(i,j,k)] = interplanetaryExpress(t1_gs(i),t2_gs(j),t3_gs(k));
        end
    end
    if round(i/length(t1_gs)*100) < 10
        fprintf('\b\b');
        fprintf('%d%%',round(i/length(t1_gs)*100));
    else
        fprintf('\b\b\b');
        fprintf('%d%%',round(i/length(t1_gs)*100));
    end
end
fprintf('\nFinished. \n')
toc

% Find minimum in dv_gs_mat
[dv_gs,loc] = min(dv_gs_mat(:));
[iii,jjj,kkk] = ind2sub(size(dv_gs_mat),loc);
t_gs = [t1_gs(iii) t2_gs(jjj) t3_gs(kkk)];      % time array of minimum dv obtained via Grid Search
dt_gs = t_gs(3)-t_gs(1);

fprintf('Results:\n')
fprintf('dv_gs = %.4f km/s\n', dv_gs)
fprintf('dt_gs = %.4f days\n', dt_gs)
fprintf('----------------------------------------------\n');
%% Methods comparison
% Comparison between the two methods and choosing of the one with a lower
% delta v cost.
if dv_gs <= dv_fmincon
    fprintf('The Grid search algorithm provided a more \nconvenient solution with respect to the \nGradient based optimization algorithm:\n');
    fprintf('Grid search:   %.4f km/s   %.4f years\n',dv_gs,dt_gs/365.25);
    fprintf('Fmincon:       %.4f km/s   %.4f years\n',dv_fmincon,dt_fmincon/365.25);
    t_chosen = t_gs;
    DELTAV = dv_gs;
    DELTAT = dt_gs;
else
    fprintf('The Gradient based optimization algorithm \nprovided a more convenient solution with \nrespect to the Grid search algorithm:\n');
    fprintf('Grid Search:    %.4f km/s    %.4f years\n',dv_gs,dt_gs/365.25);
    fprintf('fmincon:        %.4f km/s    %.4f years\n',dv_fmincon,dt_fmincon/365.25);
    t_chosen = t_fmincon;
    DELTAV = dv_fmincon;
    DELTAT = dt_fmincon;
end
t_chosen_dates = zeros(3,6);
t_chosen_dates(1,:) = mjd20002date(t_chosen(1));
t_chosen_dates(2,:) = mjd20002date(t_chosen(2));
t_chosen_dates(3,:) = mjd20002date(t_chosen(3));

%% Final results %%
% Computation of all the parameters:
[~,dv_dep,dv_flyby,dv_arr,rp_flyby,~,~,inc_flyby,vv1_t1,vv1_t2,vv2_t1,vv2_t2,rrE_flyby,vvE_flyby,kep1_t1,kep2_t1,kep1_t2,kep2_t2] = interplanetaryExpress(t_chosen(1),t_chosen(2),t_chosen(3));
[~,rp,vp_m,vp_p,v_inf_m_n,v_inf_p_n,a_hyp_m,a_hyp_p,e_hyp_m,e_hyp_p] = poweredGravityAssist(vv2_t1 - vvE_flyby,vv1_t2 - vvE_flyby);
dvtot_flyby = norm(vv1_t2 - vv2_t1);
dt_flyby = flyby_tof(a_hyp_m,e_hyp_m,a_hyp_p,e_hyp_p,rrE_flyby);

%% Contour plots calculations
t1_plot = linspace(t_chosen(1)-T_syn_DepwrtFB,t_chosen(1)+T_syn_DepwrtFB,100);
t2_plot = linspace(t_chosen(2)-T_syn_DepwrtFB,t_chosen(2)+T_syn_DepwrtFB,100);
t3_plot = linspace(t_chosen(3)-T_syn_DepwrtFB,t_chosen(3)+T_syn_DepwrtFB,100);

dv_plot1 = zeros(length(t1_plot),length(t2_plot));
dv_plot2 = zeros(length(t2_plot),length(t3_plot));

tic
fprintf('----------------------------------------------\n');
fprintf('Porkchop plot calculations running...\n');
fprintf('Progress:   0%%')
for i = 1:length(t1_plot)
    for j = 1:length(t2_plot)
        for k = 1:length(t3_plot)
            [~,dv_plot1(i,j),~] = dvLambert(t1_plot(i),t2_plot(j),Departure_planet,Flyby_planet);
            [~,~,dv_plot2(j,k)] = dvLambert(t2_plot(j),t3_plot(k),Flyby_planet,Arrival_planet);
        end
    end
    if round(i/length(t1_plot)*100) < 10
        fprintf('\b\b');
        fprintf('%d%%',round(i/length(t1_plot)*100));
    else
        fprintf('\b\b\b');
        fprintf('%d%%',round(i/length(t1_plot)*100));
    end
end
fprintf('\nFinished. \n')

t1_contour = zeros(1,length(t1_plot));
t2_contour = zeros(1,length(t1_plot));
t3_contour = zeros(1,length(t1_plot));

for i = 1:length(t1_plot)
    t1_contour(i) = datenum(mjd20002date(t1_plot(i)));
    t2_contour(i) = datenum(mjd20002date(t2_plot(i)));
    t3_contour(i) = datenum(mjd20002date(t3_plot(i)));
end
toc

%% Final results display %%
fprintf('----------------------------------------------\n');
fprintf('                FINAL RESULTS                 \n');
fprintf('----------------------------------------------\n');
fprintf('The chosen transfer window consists of:\n');
fprintf('Departure date:         %d/%d/%d %d:%d:%.2f\n',t_chosen_dates(1,3),t_chosen_dates(1,2),t_chosen_dates(1,1),t_chosen_dates(1,4),t_chosen_dates(1,5),t_chosen_dates(1,6));
fprintf('Flyby date:             %d/%d/%d %d:%d:%.2f\n',t_chosen_dates(2,3),t_chosen_dates(2,2),t_chosen_dates(2,1),t_chosen_dates(2,4),t_chosen_dates(2,5),t_chosen_dates(2,6));
fprintf('Arrival date:           %d/%d/%d %d:%d:%.2f\n',t_chosen_dates(3,3),t_chosen_dates(3,2),t_chosen_dates(3,1),t_chosen_dates(3,4),t_chosen_dates(3,5),t_chosen_dates(3,6));
fprintf('Total transfer time:%10.4f years\n',DELTAT/365.25);
fprintf('Total \x394v cost:       %10.4f km/s\n',DELTAV);
fprintf('----------------------------------------------\n');
fprintf('----------------------------------------------\n');

%% Plots %%
% Jupiter orbit
kep_Jup = uplanet (t_chosen(1), Departure_planet);
[rr_dep,vv_dep] = kep2car(kep_Jup,muS);
[rr_Jup,vv_Jup] = propagateTwoBodyOrbit(muS,[rr_dep;vv_dep],linspace(0,2*pi*sqrt(kep_Jup(1)^3/muS),1000));

% Earth orbit
kep_E = uplanet (t_chosen(2), Flyby_planet);
[rr_flyby,vv_flyby] = kep2car(kep_E,muS);
[rr_E,vv_E] = propagateTwoBodyOrbit(muS,[rr_flyby;vv_flyby],linspace(0,2*pi*sqrt(kep_E(1)^3/muS),1000));

% Mercury orbit
kep_Mer = uplanet (t_chosen(3), Arrival_planet);
[rr_arr,vv_arr] = kep2car(kep_Mer,muS);
[rr_Mer,vv_Mer] = propagateTwoBodyOrbit(muS,[rr_arr;vv_arr],linspace(0,2*pi*sqrt(kep_Mer(1)^3/muS),1000));

% Transfer orbit 1
deltat1 = (t_chosen(2)-t_chosen(1));
[rr_t1,vv_t1] = propagateTwoBodyOrbit(muS,[rr_dep;vv1_t1],linspace(0,deltat1*3600*24,1000));

% Transfer orbit 2
deltat2 = (t_chosen(3)-t_chosen(2));
[rr_t2,vv_t2] = propagateTwoBodyOrbit(muS,[rr_flyby;vv1_t2],linspace(0,deltat2*3600*24,1000));

% Astronomical units conversion
AU = astroConstants(2);     % Astronomical Unit
rr_Jup = rr_Jup./AU;
rr_E = rr_E./AU;
rr_Mer = rr_Mer./AU;
rr_dep = rr_dep./AU;
rr_flyby = rr_flyby./AU;
rr_arr = rr_arr./AU;
rr_t1 = rr_t1./AU;
rr_t2 = rr_t2./AU;

% Plot transfer
dth = 0.01;
figure()
plot3(rr_Jup(:,1),rr_Jup(:,2),rr_Jup(:,3),'--r','linewidth',2)
hold on
plot3(rr_E(:,1),rr_E(:,2),rr_E(:,3),'--k','linewidth',2)
hold on
plot3(rr_Mer(:,1),rr_Mer(:,2),rr_Mer(:,3),'--g','linewidth',2)
hold on
Plot_Planet(5,1.1443e-06,rr_dep(1),rr_dep(2),rr_dep(3));
hold on
Plot_Planet(3,1.2557e-05,rr_flyby(1),rr_flyby(2),rr_flyby(3));
hold on
Plot_Planet(1,3.2791e-05,rr_arr(1),rr_arr(2),rr_arr(3));
hold on
plot3(rr_t1(:,1),rr_t1(:,2),rr_t1(:,3),'c','linewidth',2)
hold on
plot3(rr_t2(:,1),rr_t2(:,2),rr_t2(:,3),'b','linewidth',2)
hold on
Plot_Planet(10,1.4378e-07);

axis equal
grid on

title('Interplanetary transfer strategy');
xlabel('X [AU]');
ylabel('Y [AU]');
zlabel('Z [AU]');
legend('Jupiter orbit','Earth orbit','Mercury orbit','Jupiter at departure',...
    'Earth at fly-by','Mercury at arrival','Transfer orbit 1','Transfer orbit 2');

% Plot flyby
rp_vect = rp_flyby*[1 0 0]';
vp_m_vect = vp_m*[0 1 0]';
vp_p_vect = vp_p*[0 1 0]'; 
[rr_m,vv_m] = propagateTwoBodyOrbit(muE,[rp_vect;vp_m_vect],0:-1:-1500);
[rr_p,vv_p] = propagateTwoBodyOrbit(muE,[rp_vect;vp_p_vect],0:1:1500);

figure()
plot3(rr_m(:,1),rr_m(:,2),rr_m(:,3),'c-','LineWidth',2)
hold on
plot3(rr_p(:,1),rr_p(:,2),rr_p(:,3),'b-','LineWidth',2)
hold on
Plot_Planet(3);
axis equal
grid on
legend('Incoming hyperbola','Outgoing hyperbola');
title('Fly-by');
xlabel('x');
ylabel('y');
zlabel('z');

% Contour plot Jupiter - Earth
figure()
[C1,h1] = contour(t1_contour,t2_contour,dv_plot1',min(min(dv_plot1))+(1:1:20),'HandleVisibility','off');
hold on
plot(datenum(mjd20002date(t_chosen(1))),datenum(mjd20002date(t_chosen(2))),'cp','LineWidth',2)
xlabel('Departure time');
ylabel('Arrival time');
xtickangle(45);
ytickangle(45);
datetick('x','yyyy mmm dd', 'keeplimits')
datetick('y','yyyy mmm dd', 'keeplimits')
grid on
axis equal
title('Jupiter - Earth')
hcb = colorbar;
hcb.Title.String = '$\Delta v$ [km/s]';
hcb.Title.Interpreter = 'latex';
hcb.Title.FontSize = 15;
ha = gca;
ha.FontSize = 13;
legend('Chosen time window');
% Contour plot Earth - Mercury
figure()
[C2,h2] = contour(t2_contour,t3_contour,dv_plot2',min(min(dv_plot2))+(1:1:20),'HandleVisibility','off');
hold on
plot(datenum(mjd20002date(t_chosen(2))),datenum(mjd20002date(t_chosen(3))),'bp','LineWidth',2)
xlabel('Departure time');
ylabel('Arrival time');
xtickangle(45);
ytickangle(45);
datetick('x','yyyy mmm dd', 'keeplimits')
datetick('y','yyyy mmm dd', 'keeplimits')
grid on
axis equal
title('Earth - Mercury')
hcb = colorbar;
hcb.Title.String = '$\Delta v$ [km/s]';
hcb.Title.Interpreter = 'latex';
hcb.Title.FontSize = 15;
ha = gca;
ha.FontSize = 13;
legend('Chosen time window');
