%% Main script for Planetary Explorer Mission assignment %%
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
addpath('../Provided-scripts');

% Loading of constants and input data:
% Nominal orbit:
a = 39899;      % semi-major axis [km]
e = 0.2510;     % eccentricity
i = 56.6144;    % inclination [deg]
OMEGA = 0;      % right ascension of ascending node [deg]
omega = 0;      % argument of perigee [deg]
f0 = 0;         % initial true anomaly [deg]

% Constants:
mu_E = astroConstants(13);  % gravitational parameter of the Earth [km^3/s^2]
mu_M = astroConstants(20);  % gravitational parameter of the Moon [km^3/s^2]
J2 = astroConstants(9);     % second zonal harmonic of the Earth
R_E = astroConstants(23);   % mean equatorial radius of the Earth [km]

T = 2*pi*sqrt((a^3)/mu_E);  % orbital period [s]

%% Task selection
task = 2; 
% 1 : Nominal and Repeating Ground Tracks (Perturbed and Unperturbed)
% 2 : Orbit Propagation with Perturbations (Gauss and Cartesian methods)
% 3 : Orbit Evolution 3D Representation
% 4 : Filtering of Keplerian Elements
% 5 : Comparison with real Satellite Data

switch task
    case 1
        %% Nominal Ground Tracks (Perturbed and Unperturbed)
        % Ground Track Parameters
        N = 1000;                               % number of points for the timespan discretization
        time_span1 = linspace(0,T,N);           % timespan of 1 orbit [s]
        time_span2 = linspace(0,3600*24,N);     % timespan of 1 day [s]
        time_span3 = linspace(0,3600*10*24,N);  % timespan of 10 days [s]
        
        t0 = 0;                                 % initial time [s]
        thetaG_0 = 0;                           % right ascension of Greenwich meridian at t0 [rad]
        omega_E = 15.04;                        % Earth angular rate [deg/h]
        
        % Unperturbed Ground Tracks
        [~,~,lon1,lat1] = groundTrack('Keplerian',[a,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0),mu_E],thetaG_0,time_span1,omega_E,mu_E,t0);
        [~,~,lon2,lat2] = groundTrack('Keplerian',[a,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0),mu_E],thetaG_0,time_span2,omega_E,mu_E,t0);
        [~,~,lon3,lat3] = groundTrack('Keplerian',[a,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0),mu_E],thetaG_0,time_span3,omega_E,mu_E,t0);      
        % Perturbed Ground Tracks with J2
        [~,~,lon1J2,lat1J2] = groundTrack_J2('Keplerian',[a,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0),mu_E],thetaG_0,time_span1,omega_E,mu_E,t0,J2,R_E);
        [~,~,lon2J2,lat2J2] = groundTrack_J2('Keplerian',[a,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0),mu_E],thetaG_0,time_span2,omega_E,mu_E,t0,J2,R_E);
        [~,~,lon3J2,lat3J2] = groundTrack_J2('Keplerian',[a,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0),mu_E],thetaG_0,time_span3,omega_E,mu_E,t0,J2,R_E);
        
        % Plot of unperturbed and perturbed ground tracks
        figure()
        plot_Earth_2D();
        hold on
        plot(lon1,lat1,'g');
        plot(lon1(1),lat1(1),'go','MarkerSize',8,'LineWidth',2);
        plot(lon1(end),lat1(end),'gs','MarkerSize',8,'LineWidth',2);
        plot(lon1J2,lat1J2,'r--');
        plot(lon1J2(1),lat1J2(1),'ro','MarkerSize',8,'LineWidth',2);
        plot(lon1J2(end),lat1J2(end),'rs','MarkerSize',8,'LineWidth',2);
        title('Ground Track for 1 Orbit');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Unperturbed','Start','End','Perturbed','Start','End');
        
        figure()
        plot_Earth_2D();
        hold on
        plot(lon2,lat2,'g');
        plot(lon2(1),lat2(1),'go','MarkerSize',8,'LineWidth',2);
        plot(lon2(end),lat2(end),'gs','MarkerSize',8,'LineWidth',2);
        plot(lon2J2,lat2J2,'r--');
        plot(lon2J2(1),lat2J2(1),'ro','MarkerSize',8,'LineWidth',2);
        plot(lon2J2(end),lat2J2(end),'rs','MarkerSize',8,'LineWidth',2);
        title('Ground Track for 1 Day');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Unperturbed','Start','End','Perturbed','Start','End');
        
        figure()
        plot_Earth_2D();
        hold on
        plot(lon3,lat3,'g');
        plot(lon3(1),lat3(1),'go','MarkerSize',8,'LineWidth',2);
        plot(lon3(end),lat3(end),'gs','MarkerSize',8,'LineWidth',2);
        plot(lon3J2,lat3J2,'r--');
        plot(lon3J2(1),lat3J2(1),'ro','MarkerSize',8,'LineWidth',2);
        plot(lon3J2(end),lat3J2(end),'rs','MarkerSize',8,'LineWidth',2);
        title('Ground Track for 10 Days');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Unperturbed','Start','End','Perturbed','Start','End');
        
        %% Repeating Ground Tracks (Perturbed and Unperturbed)
        k = 1; % number of revolutions of the s/c to obtain the ground track repetition
        m = 1; % number of revolutions of the planet to obtain the ground track repetition
        
        % Semi-major axis for unperturbed repeating ground track
        a_repeated = a_repeated_ground_track(k,m,mu_E,deg2rad(omega_E/3600));   % a for repeating ground track [km]
        T_rep = 2*pi*sqrt((a_repeated^3)/mu_E);                                 % modified orbital period of the repeating orbit [s]
        time_span1_rep = linspace(0,T_rep,N);                                   % timespan for 1 repeating orbit [s]
         
        [~,~,lon1_rep,lat1_rep] = groundTrack('Keplerian',[a_repeated,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0)],thetaG_0,time_span1_rep,omega_E,mu_E,t0);
        [~,~,lon2_rep,lat2_rep] = groundTrack('Keplerian',[a_repeated,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0)],thetaG_0,time_span2,omega_E,mu_E,t0);
        [~,~,lon3_rep,lat3_rep] = groundTrack('Keplerian',[a_repeated,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0)],thetaG_0,time_span3,omega_E,mu_E,t0);
        
        % Semi-major axis for perturbed repeating ground track
        a_repeated_pert = a_repeated_ground_track_J2(k,m,mu_E,deg2rad(omega_E/3600),e,J2,R_E,deg2rad(i),a_repeated); % a for repeating and perturbed ground track [km]
        T_rep_pert = 2*pi*sqrt((a_repeated_pert^3)/mu_E);           % modified orbital period of the repeating and perturbed orbit [s]
        time_span1_rep_pert = linspace(0,T_rep_pert,N);             % timespan for 1 repeating and perturbed orbit [s]       
        
        [~,~,lon1_repJ2,lat1_repJ2] = groundTrack_J2('Keplerian',[a_repeated_pert,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0)],thetaG_0,time_span1_rep_pert,omega_E,mu_E,t0,J2,R_E);
        [~,~,lon2_repJ2,lat2_repJ2] = groundTrack_J2('Keplerian',[a_repeated_pert,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0)],thetaG_0,time_span2,omega_E,mu_E,t0,J2,R_E);
        [~,~,lon3_repJ2,lat3_repJ2] = groundTrack_J2('Keplerian',[a_repeated_pert,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0)],thetaG_0,time_span3,omega_E,mu_E,t0,J2,R_E);

        % Plot of unperturbed and perturbed repeating ground tracks:
        figure()
        plot_Earth_2D();
        hold on
        plot(lon1_rep,lat1_rep,'g');
        plot(lon1_rep(1),lat1_rep(1),'go','MarkerSize',8,'LineWidth',2);
        plot(lon1_rep(end),lat1_rep(end),'gs','MarkerSize',8,'LineWidth',2);
        plot(lon1_repJ2,lat1_repJ2,'r--');
        plot(lon1_repJ2(1),lat1_repJ2(1),'ro','MarkerSize',8,'LineWidth',2);
        plot(lon1_repJ2(end),lat1_repJ2(end),'rs','MarkerSize',8,'LineWidth',2);
        title('Repeating Ground Track for 1 Orbit');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Unperturbed','Start','End','Perturbed','Start','End');
        
        figure()
        plot_Earth_2D();
        hold on
        plot(lon2_rep,lat2_rep,'g');
        plot(lon2_rep(1),lat2_rep(1),'ro','MarkerSize',8,'LineWidth',2);
        plot(lon2_rep(end),lat2_rep(end),'gs','MarkerSize',8,'LineWidth',2);
        plot(lon2_repJ2,lat2_repJ2,'r--');
        plot(lon2_repJ2(1),lat2_repJ2(1),'ro','MarkerSize',8,'LineWidth',2);
        plot(lon2_repJ2(end),lat2_repJ2(end),'rs','MarkerSize',8,'LineWidth',2);
        title('Repeating Ground Track for 1 Day');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Unperturbed','Start','End','Perturbed','Start','End');
        
        figure()
        plot_Earth_2D();
        hold on
        plot(lon3_rep,lat3_rep,'g');
        plot(lon3_rep(1),lat3_rep(1),'go','MarkerSize',8,'LineWidth',2);
        plot(lon3_rep(end),lat3_rep(end),'gs','MarkerSize',8,'LineWidth',2);
        plot(lon3_repJ2,lat3_repJ2,'r');
        plot(lon3_repJ2(1),lat3_repJ2(1),'ro','MarkerSize',8,'LineWidth',2);
        plot(lon3_repJ2(end),lat3_repJ2(end),'rs','MarkerSize',8,'LineWidth',2);
        title('Repeating Ground Track for 10 Days');
        xlabel('Longitude [deg]');
        ylabel('Latitude [deg]');
        legend('Unperturbed','Start','End','Perturbed','Start','End');
        
    case 2
        %% Orbit Propagation with Perturbations (Gauss and Cartesian methods)
        kep_0 = [a,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0)];     % initial keplerian elements
        x_car_0 = ones(6,1);
        [x_car_0(1:3),x_car_0(4:6)] = kep2car(kep_0,mu_E);        % initial cartesian orbital state
        date_in = [2020,10,31,12,00,00];                                        % initial date for propagation
        t0 = date2mjd2000(date_in);                                             % initial date in MJD2000 [days]
        t00=t0*86400;                                                           % inidial date in MJD2000 converted in seconds [s]
        k_prop = 500;                                                           % number of orbits for the propagation
        N_prop = 10000;                                                         % number of points for the discretization of the timespan
        t_span = linspace(t00,t00 + k_prop*T ,N_prop);                          % timespan for propagation [s]
        options = odeset('RelTol',1e-13,'AbsTol',1e-14); 
        
        % Gauss propagation:
        GAUSS_PLANETARY_EQs = @(t,kep) ode_gauss_rsw_asgn (t,kep,mu_E,Keplerian_model_aj2_RSW(kep,mu_E,R_E,J2),TimeEph_model_a_moon_RSW(t,kep,mu_M,mu_E));
        [T_Gauss,kep_gauss] = ode113(GAUSS_PLANETARY_EQs,t_span,kep_0,options);
        
        % Cartesian propagation:
        [T_Car,x_car] = ode113(@(t,x) kepl_orbit_J2_moon(x,t,mu_E,mu_M,J2,R_E),t_span,x_car_0,options);
        kep_car = zeros(N_prop,6); 
        for i=1:N_prop                                                          % Conversion from cartesian orbital states to keplerian elements
            kep_car(i,:) = car2kep(x_car(i,1:3),x_car(i,4:6),mu_E);
        end
        kep_car(:,4:6) = unwrap(kep_car(:,4:6));                                % Unwrapping of angular keplerian elements
        
        % Errors between Gauss and Cartesian propagation:
        err_a_plot = abs(kep_car(:,1)-kep_gauss(:,1))./abs(kep_0(1));           % relative error of semi-major axis
        err_e_plot = abs(kep_car(:,2)-kep_gauss(:,2));                          % absolute error of eccentricity
        err_i_plot = abs(kep_car(:,3)-kep_gauss(:,3))/(2*pi);                   % relative error of inclination
        err_OM_plot = abs(kep_car(:,4)-kep_gauss(:,4))/(2*pi);                  % relative error of right ascension of ascending node
        err_om_plot = abs(kep_car(:,5)-kep_gauss(:,5))/(2*pi);                  % relative error of argoument of periapsis
        err_th_plot = abs(kep_car(:,6)-kep_gauss(:,6))./abs(kep_gauss(:,6));    % relative error of true anomaly

        % Plot of comparison between keplerian elements with Gauss and
        % Cartesian propagation
        figure()
        subplot(2,3,1)
        plot((T_Car-t00)/T,kep_car(:,1),'LineWidth',2)
        title('evolution of a _ Cartesian method [km]');
        xlabel('n° orbit period');
        ylabel('a_{car} [km]');
        grid on
        
        subplot(2,3,2)
        plot((T_Gauss-t00)/T,kep_gauss(:,1),'LineWidth',2)
        title('evolution of a _ Gauss method [km]');
        xlabel('n° orbit period');
        ylabel('a_{gauss} [km]');
        grid on
        
        subplot(2,3,3)
        semilogy((T_Gauss-t00)/T,err_a_plot,'LineWidth',2)
        title('rel. err. of a ');
        xlabel('n° orbit period');
        ylabel('|a_{car} - a_{gauss}|/|a_0| ');
        grid on
        
        subplot(2,3,4)
        plot((T_Car-t00)/T,kep_car(:,2),'LineWidth',2)
        title('evolution of e _ Cartesian method [-]');
        xlabel('n° orbit period');
        ylabel('e_{car} [km]');
        grid on
        
        subplot(2,3,5)
        plot((T_Gauss-t00)/T,kep_gauss(:,2),'LineWidth',2)
        title('evolution of e _ Gauss method [-]');
        xlabel('n° orbit period');
        ylabel('e_{gauss} [-]');
        grid on
        
        subplot(2,3,6)
        semilogy((T_Gauss-t00)/T,err_e_plot,'LineWidth',2)
        title('abs. err. of e [-]');
        xlabel('n° orbit period');
        ylabel('|e_{car} - e_{gauss}|  [-]');
        grid on
       
        figure()
        subplot(2,3,1)
        plot((T_Car-t00)/T,rad2deg(kep_car(:,3)),'LineWidth',2)
        title('evolution of i _ Cartesian method [deg]');
        xlabel('n° orbit period');
        ylabel('i_{car} [deg]');
        grid on
        
        subplot(2,3,2)
        plot((T_Gauss-t00)/T,rad2deg(kep_gauss(:,3)),'LineWidth',2)
        title('evolution of i _ Gauss method [deg]');
        xlabel('n° orbit period');
        ylabel('i_{gauss} [deg]');
        grid on
        
        subplot(2,3,3)
        semilogy((T_Gauss-t00)/T,rad2deg(err_i_plot),'LineWidth',2)
        title('rel. err. of i');
        xlabel('n° orbit period');
        ylabel('|i_{car} - i_{gauss}|/360° ');
        grid on
        
        subplot(2,3,4)
        plot((T_Car-t00)/T,rad2deg(kep_car(:,4)),'LineWidth',2)
        title('evolution of \Omega [deg] _ Cartesian method');
        xlabel('n° orbit period');
        ylabel('\Omega_{car} [deg]');
        grid on
        
        subplot(2,3,5)
        plot((T_Gauss-t00)/T,rad2deg(kep_gauss(:,4)),'LineWidth',2)
        title('evolution of \Omega [deg] _ Gauss method');
        xlabel('n° orbit period');
        ylabel('\Omega_{gauss} [deg]');
        grid on
        
        subplot(2,3,6)
        semilogy((T_Gauss-t00)/T,rad2deg(err_OM_plot),'LineWidth',2)
        title('rel. err. of \Omega ');
        xlabel('n° orbit period');
        ylabel('|\Omega_{car} - \Omega_{gauss}|/360°');
        grid on
                
        figure()
        subplot(2,3,1)
        plot((T_Car-t00)/T,rad2deg(kep_car(:,5)),'LineWidth',2)
        title('evolution of \omega _ Cartesian method [deg]');
        xlabel('n° orbit period');
        ylabel('\omega_{car} [deg]');
        grid on
        
        subplot(2,3,2)
        plot((T_Gauss-t00)/T,rad2deg(kep_gauss(:,5)),'LineWidth',2)
        title('evolution of \omega _ Gauss method [deg]');
        xlabel('n° orbit period');
        ylabel('\omega_{gauss} [deg]');
        grid on
        
        subplot(2,3,3)
        semilogy((T_Gauss-t00)/T,rad2deg(err_om_plot),'LineWidth',2)
        title('rel. err. of \omega');
        xlabel('n° orbit period');
        ylabel('|\omega_{car} - \omega_{gauss}|/360°');
        grid on
        
        subplot(2,3,4)
        plot((T_Car-t00)/T,rad2deg(kep_car(:,6)),'LineWidth',2)
        title('evolution of \theta _ Cartesian method [deg]');
        xlabel('n° orbit period');
        ylabel('\theta_{car} [deg]');
        grid on
        
        subplot(2,3,5)
        plot((T_Gauss-t00)/T,rad2deg(kep_gauss(:,6)),'LineWidth',2)
        title('evolution of \theta _ Gauss method [deg]');
        xlabel('n° orbit period');
        ylabel('\theta_{gauss} [deg]');
        grid on
        
        subplot(2,3,6)
        semilogy((T_Gauss-t00)/T,rad2deg(err_th_plot),'LineWidth',2)
        title('rel. err. of \theta');
        xlabel('n° orbit period');
        ylabel('|\theta_{car} - \theta_{gauss}|/360°');
        grid on
    case 3 
        %% Orbit Evolution 3D Representation
        kep_0 = [a,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0)];     % initial keplerian elements for propagation
        date_in = [2020,10,31,12,00,00];                                        % initial date for propagation
        t0 = date2mjd2000(date_in);                                             % initial date in mjd2000
        t00=t0*86400;                                                           % conversion of mjd2000 in seconds [s]
        k_prop = 10000;                                                         % number of orbits for propagation
        N_prop = 50*k_prop;                                                     % number of points for the timespan discretization
        fprintf('Orbit evolution representation calculations in progress...\n');
        tic
        t_span = linspace(t00,t00 + k_prop*T ,N_prop);                          % timespan for propagation [s]
        options = odeset('RelTol',1e-13,'AbsTol',1e-14);
        GAUSS_PLANETARY_EQs = @(t,kep) ode_gauss_rsw_asgn (t,kep,mu_E,Keplerian_model_aj2_RSW(kep,mu_E,R_E,J2),TimeEph_model_a_moon_RSW(t,kep,mu_M,mu_E));
        [T_Gauss,kep_gauss] = ode113(GAUSS_PLANETARY_EQs,t_span,kep_0,options); % propagation of the orbit through Gauss equations
        periods_change_color=1;                                                 % number which indicates after how many period the change of color must occur
        plot_perturbed_orbit(T_Gauss,kep_gauss,T,mu_E,N_prop,k_prop,periods_change_color); % plots the 3d evolution of the orbit
        toc
        fprintf('Finished.\n');
    case 4
        %% Filtering of Keplerian Elements
        kep_0 = [a,e,deg2rad(i),deg2rad(OMEGA),deg2rad(omega),deg2rad(f0)]; % initial keplerian elements
        date_in = [2020,10,31,12,00,00];                                    % initial date for propagation
        t0 = date2mjd2000(date_in);                                         % initial date in mjd2000
        t00=t0*86400;                                                       % conversion of mjd2000 in seconds [s]
        k_prop = 500;                                                       % number of orbits for propagation
        N_prop = 100000;                                                    % number of points for the timespan discretization
        t_span = linspace(t00,t00 + k_prop*T ,N_prop);                      % timespan for propagation [s]
        options = odeset('RelTol',1e-13,'AbsTol',1e-14);
        
        % Gauss propagation:
        GAUSS_PLANETARY_EQs = @(t,kep) ode_gauss_rsw_asgn (t,kep,mu_E,Keplerian_model_aj2_RSW(kep,mu_E,R_E,J2),TimeEph_model_a_moon_RSW(t,kep,mu_M,mu_E));
        [T_Gauss,kep_gauss] = ode113(GAUSS_PLANETARY_EQs,t_span,kep_0,options);
        T_short = T;                                                        % period of short-periodic oscillations [s]
        T_long = 15*T;                                                      % period of long-periodic oscillations [s]
        nwindow_short = nearest(T_short/(sum(diff(T_Gauss)/(numel(T_Gauss)-1))));       % number of points of the time window of the length of 1 T_short [s]
        nwindow_long = nearest(T_long/(sum(diff(T_Gauss)/(numel(T_Gauss)-1))));         % number of points of the time window of the length of 1 T_long [s]
        kep_filtered_long = movmean(kep_gauss,nwindow_short,1);             % filtering of short-periodic oscillations (only secular and long periodic components remaining)
        kep_filtered_secular = movmean(kep_gauss,nwindow_long,1);           % filtering of long-periodic oscillations (only secular components remaining)
        
        % Plot of filtered keplerian elements:
        figure()
        subplot(3,2,1)
        plot((T_Gauss-t00)/T,kep_gauss(:,1),'LineWidth',2);
        hold on
        plot((T_Gauss-t00)/T,kep_filtered_secular(:,1),'g','LineWidth',2);
        xlabel('n° orbit period');
        ylabel(' a [km]');
        legend('Complete','Secular');
        title('Semi-major axis');
        grid on
        
        subplot(3,2,2)
        plot((T_Gauss-t00)/T,kep_gauss(:,2),'LineWidth',2);
        hold on
        plot((T_Gauss-t00)/T,kep_filtered_long(:,2),'r','LineWidth',2);
        plot((T_Gauss-t00)/T,kep_filtered_secular(:,2),'g','LineWidth',2);
        xlabel('n° orbit period');
        ylabel(' e [-] ');
        legend('Complete','Long Term','Secular');
        title('Eccentricity');
        grid on
        
        subplot(3,2,3)
        plot((T_Gauss-t00)/T,rad2deg(kep_gauss(:,3)),'LineWidth',2);
        hold on
        plot((T_Gauss-t00)/T,rad2deg(kep_filtered_long(:,3)),'r','LineWidth',2);
        plot((T_Gauss-t00)/T,rad2deg(kep_filtered_secular(:,3)),'g','LineWidth',2);
        xlabel('n° orbit period');
        ylabel('i [deg]');
        legend('Complete','Long Term','Secular');
        title('Inclination');
        grid on
        
        subplot(3,2,4)
        plot((T_Gauss-t00)/T,rad2deg(kep_gauss(:,4)),'LineWidth',2);
        hold on
        plot((T_Gauss-t00)/T,rad2deg(kep_filtered_long(:,4)),'r','LineWidth',2);
        plot((T_Gauss-t00)/T,rad2deg(kep_filtered_secular(:,4)),'g','LineWidth',2);
        xlabel('n° orbit period');
        ylabel('\Omega [deg]');
        legend('Complete','Long Term','Secular');
        title('Right Ascension of Ascending Node');
        grid on
        
        subplot(3,2,5)
        plot((T_Gauss-t00)/T,rad2deg(kep_gauss(:,5)),'LineWidth',2);
        hold on
        plot((T_Gauss-t00)/T,rad2deg(kep_filtered_long(:,5)),'r','LineWidth',2);
        plot((T_Gauss-t00)/T,rad2deg(kep_filtered_secular(:,5)),'g','LineWidth',2);
        xlabel('n° orbit period');
        ylabel('\omega [deg]');
        legend('Complete','Long Term','Secular');
        title('Argument of Perigee');
        grid on
        
        subplot(3,2,6)
        plot((T_Gauss-t00)/T,rad2deg(kep_gauss(:,6)),'LineWidth',2);
        hold on
        plot((T_Gauss-t00)/T,rad2deg(kep_filtered_long(:,6)),'r','LineWidth',2);
        plot((T_Gauss-t00)/T,rad2deg(kep_filtered_secular(:,6)),'g','LineWidth',2);
        xlabel('n° orbit period');
        ylabel('\theta [deg]');
        legend('Complete','Long Term','Secular');
        title('True Anomaly');
        grid on

    case 5
        %% Real Satellite Orbit Comparison
        % Satellite used:    (NORAD Catalog number: 24282 ; Name: ZHONGXING-7 )
              
        % Initial keplerian elements of the satellite and other inputs:
        kep0_sat = [4.046496437446377E+04, 2.917034269288031E-01, deg2rad(2.183757731769000E+01), deg2rad(2.616733747605937E+02), deg2rad(1.500679599533694E+02), deg2rad(1.372205091124491E+02)] ;
        
        initial_date = [2019,01,01,00,00,00];           % initial date for timespan
        final_date = [2020,01,01,00,00,00];             % final date for timespan
        t_in = date2mjd2000(initial_date);              % initial date in MJD2000
        t_fin = date2mjd2000(final_date);               % final date in MJD2000
        t_in_sec = t_in*86400;                          % initial MJD2000 in seconds (for propagation purposes)
        t_fin_sec = t_fin*86400;                        % initial MJD2000 in seconds (for propagation purposes)
        Step = 86400;                                   % length in seconds of the time-step for propagation, chosen equal to the one of the ephemerides (1 day = 86400s)
        t_span_sat = [t_in_sec:Step:t_fin_sec];         % defintion of the timespan

        options = odeset('RelTol',1e-13,'AbsTol',1e-14);
        
        % Propagation of the Ideal Prturbed Model 
        GAUSS_PLANETARY_EQs = @(t,kep) ode_gauss_rsw_asgn (t,kep,mu_E,Keplerian_model_aj2_RSW(kep,mu_E,R_E,J2),TimeEph_model_a_moon_RSW(t,kep,mu_M,mu_E));
        [Time_sat,kep_sat] = ode113(GAUSS_PLANETARY_EQs,t_span_sat,kep0_sat,options);
        
        % Output of Real Ephemerides
        namefileexcel= 'Real_ephemerides_ICRF-CUTforMATLAB.xlsx';       % excel file containg the ephemerides data             
        [time_real_eph,kep_real_eph]= Real_ephemerides(namefileexcel);  % extracts the time and keplerian elements from the ephemerides file
        T_sat= 2*pi*sqrt(kep0_sat(1)^3/mu_E);                           % orbital period of the satellite [s]
        
        kep_real_eph(:,3:6) = deg2rad(kep_real_eph(:,3:6));             %converts keplerian elements in [rad]
        kep_real_eph(:,4:6) = unwrap(kep_real_eph(:,4:6));              %unwraps keplerian elements obtained from the ephemerides
        kep_sat(:,4:6) = unwrap(kep_sat(:,4:6));                        %unwraps keplerian elements obtained from the propagation
        
        % Errors between Gauss and Cartesian propagation
        err_a_plot1 = abs(kep_sat(:,1)-kep_real_eph(:,1))./abs(kep0_sat(1));         % relative error of semi-major axis
        err_e_plot1 = abs(kep_sat(:,2)-kep_real_eph(:,2));                           % absolute error of eccentricity
        err_i_plot1 = abs(kep_sat(:,3)-kep_real_eph(:,3))/(2*pi);                    % relative error of inclination
        err_OM_plot1 = abs(kep_sat(:,4)-kep_real_eph(:,4))/(2*pi);                   % relative error of right ascension of ascending node
        err_om_plot1 = abs(kep_sat(:,5)-kep_real_eph(:,5))/(2*pi);                   % relative error of argoument of periapsis
        err_th_plot1 = abs(kep_sat(:,6)-kep_real_eph(:,6))./abs(kep_real_eph(:,6));  % relative error of true anomaly
        
        % Maximum errors (absolute and relative): 
        % relative errors:
        max_err_a = max(err_a_plot1);
        max_err_e = max(err_e_plot1);
        max_err_i = max(err_i_plot1);
        max_err_OM = max(err_OM_plot1);
        max_err_om = max(err_om_plot1);
        max_err_th = max(err_th_plot1);
        
        % absolute errors:
        max_errass_a = max(abs(kep_sat(:,1)-kep_real_eph(:,1)));
        max_errass_e = max(abs(kep_sat(:,2)-kep_real_eph(:,2)));
        max_errass_i = rad2deg(max(abs(kep_sat(:,3)-kep_real_eph(:,3))));
        max_errass_OM = rad2deg(max(abs(kep_sat(:,4)-kep_real_eph(:,4))));
        max_errass_om = rad2deg(max(abs(kep_sat(:,5)-kep_real_eph(:,5))));
        max_errass_th = rad2deg(max(abs(kep_sat(:,6)-kep_real_eph(:,6))));
        
        % plot of comparison between keplerian elements obtained with the numerical propagation
        % and the ephemerides data:
        
        figure()
        subplot(3,2,1)
        plot((Time_sat-t_in_sec)/T_sat,kep_sat(:,1),'LineWidth',2)
        hold on
        plot((time_real_eph-t_in_sec)/T_sat,kep_real_eph(:,1),'r','LineWidth',2)
        xlabel('n° orbit period');
        ylabel('a [km]');
        legend('Gauss propagation','Real Data');
        title('Semi-major Axis');
        grid on
        
        subplot(3,2,2)
        plot((Time_sat-t_in_sec)/T_sat,kep_sat(:,2),'LineWidth',2)
        hold on
        plot((time_real_eph-t_in_sec)/T_sat,kep_real_eph(:,2),'r','LineWidth',2)
        xlabel('n° orbit period');
        ylabel('e [-]');
        legend('Gauss propagation','Real Data');        
        title('Eccentricity');
        grid on
        
        subplot(3,2,3)
        plot((Time_sat-t_in_sec)/T_sat,rad2deg(kep_sat(:,3)),'LineWidth',2)
        hold on
        plot((time_real_eph-t_in_sec)/T_sat,rad2deg(kep_real_eph(:,3)),'r','LineWidth',2)
        xlabel('n° orbit period');
        ylabel('i [deg]');
        legend('Gauss propagation','Real Data');
        title('Inclination');
        grid on
        
        subplot(3,2,4)
        plot((Time_sat-t_in_sec)/T_sat,rad2deg(kep_sat(:,4)),'LineWidth',2)
        hold on
        plot((time_real_eph-t_in_sec)/T_sat,rad2deg(kep_real_eph(:,4)),'r','LineWidth',2)
        xlabel('n° orbit period');
        ylabel('\Omega [deg]');
        legend('Gauss propagation','Real Data');
        title('Right Ascension of Ascending Node');
        grid on
        
        subplot(3,2,5)
        plot((Time_sat-t_in_sec)/T_sat,rad2deg(kep_sat(:,5)),'LineWidth',2)
        hold on
        plot((time_real_eph-t_in_sec)/T_sat,rad2deg(kep_real_eph(:,5)),'r','LineWidth',2)
        xlabel('n° orbit period');
        ylabel('\omega [deg]');
        legend('Gauss propagation','Real Data');
        title('Argument of Pericenter');
        grid on
        
        subplot(3,2,6)
        plot((Time_sat-t_in_sec)/T_sat,rad2deg(kep_sat(:,6)),'LineWidth',2)
        hold on
        plot((time_real_eph-t_in_sec)/T_sat,rad2deg(kep_real_eph(:,6)),'r','LineWidth',2)
        xlabel('n° orbit period');
        ylabel('\theta [deg]');
        legend('Gauss propagation','Real Data');
        title('True Anomaly');
        grid on
        
end

