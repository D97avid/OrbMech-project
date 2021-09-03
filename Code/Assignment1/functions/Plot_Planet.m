function [] = Plot_Planet(ID,A,xx,yy,zz)
% 
% function for 3D plot of the solar system planets and Sun.
% 
% PROTOTYPE:
%  [] = Plot_Planet(ID,A,xx,yy,zz)
% 
% INPUT:
%  ID [1]    planet identification number                             [-]
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%                   10:  Sun
%  A [1]     size modifier for plotting Earth                         [-]
%  xx [1]    x-axis coordinate in inertial frame                      [km]
%  yy [1]    y-axis coordinate in inertial frame                      [km]
%  zz [1]    z-axis coordinate in inertial frame                      [km]
% 
% OUTPUT:
% n.a.
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

if nargin == 1
    xx = 0;
    yy = 0;
    zz = 0;
    A = 1;
elseif nargin == 2
    xx = 0;
    yy = 0;
    zz = 0;
end

switch ID
    case 1
        C = imread('mercury.jpg');
        [x, y, z] = ellipsoid(xx,yy,zz, astroConstants(21)*A, astroConstants(21)*A,astroConstants(21)*A,1E2);
        surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2))]), 'FaceColor', 'texturemap','EdgeColor','none');
    case 2
        C = imread('venus.jpg');
        [x, y, z] = ellipsoid(xx,yy,zz, astroConstants(22)*A, astroConstants(22)*A,astroConstants(22)*A,1E2);
        surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2))]), 'FaceColor', 'texturemap','EdgeColor','none');
    case 3
        C = imread('earth.jpg');
        [x, y, z] = ellipsoid(xx,yy,zz, astroConstants(23)*A, astroConstants(23)*A,astroConstants(23)*A,1E2);
        surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2))]), 'FaceColor', 'texturemap','EdgeColor','none');
    case 4
        C = imread('mars.jpg');
        [x, y, z] = ellipsoid(xx,yy,zz, astroConstants(24)*A, astroConstants(24)*A,astroConstants(24)*A,1E2);
        surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2))]), 'FaceColor', 'texturemap','EdgeColor','none');
    case 5
        C = imread('jupiter.jpg');
        [x, y, z] = ellipsoid(xx,yy,zz, astroConstants(25)*A, astroConstants(25)*A,astroConstants(25)*A,1E2);
        surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2))]), 'FaceColor', 'texturemap','EdgeColor','none');
    case 6
        C = imread('saturn.jpg');
        [x, y, z] = ellipsoid(xx,yy,zz, astroConstants(26)*A, astroConstants(26)*A,astroConstants(26)*A,1E2);
        surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2))]), 'FaceColor', 'texturemap','EdgeColor','none');
    case 7
        C = imread('uranus.jpg');
        [x, y, z] = ellipsoid(xx,yy,zz, astroConstants(27)*A, astroConstants(27)*A,astroConstants(27)*A,1E2);
        surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2))]), 'FaceColor', 'texturemap','EdgeColor','none');
    case 8
        C = imread('neptune.jpg');
        [x, y, z] = ellipsoid(xx,yy,zz, astroConstants(28)*A, astroConstants(28)*A,astroConstants(28)*A,1E2);
        surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2))]), 'FaceColor', 'texturemap','EdgeColor','none');
    case 9
        C = imread('pluto.jpg');
        [x, y, z] = ellipsoid(xx,yy,zz, astroConstants(29)*A, astroConstants(29)*A,astroConstants(29)*A,1E2);
        surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2))]), 'FaceColor', 'texturemap','EdgeColor','none');
    case 10
        C = imread('sun.jpg');
        [x, y, z] = ellipsoid(xx,yy,zz, astroConstants(3)*A, astroConstants(3)*A,astroConstants(3)*A,1E2);
        surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2))]), 'FaceColor', 'texturemap','EdgeColor','none');
end
