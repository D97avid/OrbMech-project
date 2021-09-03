function plot_Earth_2D()
% 
% Function to plot the 2D representation of Earth's surface for Ground Track plotting.
% 
% PROTOTYPE:
% plot_Earth_2D()
%
% INPUT:
%  [-]
% 
% OUTPUT:
%  [-]
%
% CONTRIBUTORS:
%  Marco Adorno
%  Giuseppe Esposito 
%  Davide Gravina 
%  David Reina
% 
% VERSIONS:
%  20-01-2021: First version

I = (imread('Earth.jpg'));          %loads the Earth map
image([-180 180],[90 -90], I);      %puts the axes limits as the latitude and longitude limits (in deg)
axis xy
