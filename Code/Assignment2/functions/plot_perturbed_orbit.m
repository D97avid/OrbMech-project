function [] = plot_perturbed_orbit(T_Gauss,kep_gauss,T,mu_E,N_prop,k_prop,periods_change_color)
% 
% Function that computes the plot of propagation orbits in all time period
% considerated, in Cartesian Coordinates. The figure give also a scaling
% change of color between the initial orbit and the final one.
% 
% PROTOTYPE:
%  plot_perturbed_orbit(T_Gauss,kep_gauss,T,mu_E,N_prop,k_prop,periods_change_color)
%
% INPUT: 
% T_Gauss [Nx1]                 Time array of integration of ODE Solver
% kep_gauss [Nx6]               Kepler parameters of N points (a,e,i,OM,om,th) [km,-,rad,rad,rad,rad]
% T [1]                         Orbital Period of orbit considerated [s]
% mu_E [1]                      Earth's gravitational parameter [km^3/s^2] 
% N_prop [1]                    Number of point considerated, inside tSpan, to propagate the orbit.
% k_prop [1]                    Number of revolution considerated for the propagation
% periods_change_color [1]      number of periods considered beforge change color in the plot.
%
% OUTPUT: 
%            
%
% CONTRIBUTORS:
%  Marco Adorno
%  Giuseppe Esposito 
%  Davide Gravina 
%  David Reina
% 
% VERSIONS:
%  20-01-2021: First version


%% GAUSS PLOT of Pertubated orbit:
% Transformation from keplerian elements to Cartesian for plot: (PAY ATTENTION about which 'kep2car' are using)
x_gauss = zeros(N_prop,3);

for i=1:N_prop
    [x_gauss(i,1:3),~] = kep2car (kep_gauss(i,:),mu_E);
end

hold on
title('Perturbed orbit for J2 and Moon');
grid on

% Planet and orbit drawing:
%Set of step to change the color of plot (for each period of orbit):
num_orbit=k_prop;
period_max=periods_change_color; %orbit period considereted before change the color
%initialization of index:
index1=1;
index2=zeros(ceil(num_orbit/period_max),1);
%computation of indexSteps for changing color plot:
for ii=1:1:length(x_gauss(:,1))
    if T_Gauss(ii)>(index1*T*period_max)+T_Gauss(1,1)
        index2(index1)=ii;
        index1=index1+1; 
    end      
end
%plot with different color for each period of orbit:
% initial plot:
plot3(x_gauss(1:index2(1),1),x_gauss(1:index2(1),2),x_gauss(1:index2(1),3),'--b','LineWidth',2.5); %plot of s/c orbit
% middle plots:
for ii=1:1:index1-2
pip=plot3(x_gauss(index2(ii):index2(ii+1),1),x_gauss(index2(ii):index2(ii+1),2),x_gauss(index2(ii):index2(ii+1),3),'color',[(ii/num_orbit),0.275,1-(ii/num_orbit)]); %plot of s/c orbit
pip.HandleVisibility = 'off'; %delete these plot from the legend
end
% final plot: 
plot3(x_gauss(index2(end-1):end,1),x_gauss(index2(end-1):end,2),x_gauss(index2(end-1):end,3),'--r','LineWidth',2.5); %plot of s/c orbit

axis equal
% Starting and Ending Points of s/c:
plot3(x_gauss(1,1),x_gauss(1,2),x_gauss(1,3),'ob','LineWidth',2)         %Pointer for starting position of s/c
plot3(x_gauss(end,1),x_gauss(end,2),x_gauss(end,3),'sr','LineWidth',2)   %Pointer for ending position of s/c
hold on
Plot_Planet(3)
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('Initial Orbit','Final Orbit','Starting Point','Ending Point')


end
