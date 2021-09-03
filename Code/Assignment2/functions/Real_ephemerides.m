function [time_real_eph,kep_real_eph]= Real_ephemerides(namefileexcel)
% 
% This function convert the data of real ephemerides into the excel sheet
% into a Matlab matrix.
%
% PROTOTYPE:
% [time_real_eph,kep_real_eph]= Real_ephemerides(namefileexcel)
%
% INPUT:
% namefileexcel [string]      name of file excel to convert in matlab matrix
% 
% OUTPUT: 
% time_real_eph [1]             time of Modified Julian day 2000 number [s]
% kep_real_eph [1x6]            [a,e,i,OM,om,th] Keplerian parameters [km,deg]
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




% HOW TO SET THE DATA INTO THE FILE EXCEL :
% In this case the column vector into the matrix correspond at this
% keplerian parameter, into ICRF (International Celestial Reference Frame) frame centred into Earth 
% 
%   M_excel[:,1] =   JDTDB    Julian Day Number, Barycentric Dynamical Time
%   M_excel[:,2] =   EC     Eccentricity, e
%   M_excel[:,3] =   QR     Periapsis distance, q (km)
%   M_excel[:,4] =   IN     Inclination w.r.t X-Y plane, i (degrees)
%   M_excel[:,5] =   OM     Longitude of Ascending Node, OMEGA, (degrees)
%   M_excel[:,6] =   W      Argument of Perifocus, w (degrees)
%   M_excel[:,7] =   Tp     Time of periapsis (Julian Day Number)
%   M_excel[:,8] =   N      Mean motion, n (degrees/sec)
%   M_excel[:,9] =   MA     Mean anomaly, M (degrees)
%   M_excel[:,10] =   TA     True anomaly, nu (degrees)
%   M_excel[:,11] =   A      Semi-major axis, a (km)
%   M_excel[:,12] =   AD     Apoapsis distance (km)
%   M_excel[:,13] =   PR     Sidereal orbit period (sec)
% 
% 
%REMEMBER_1: the file excel must be into the path folder of functions and
%must be correct and caracterize particularly to fit into the funcion used
%to read it.
% (CHECK THE 'Real_ephemerides_ICRF-CUTforMATLAB.xlsx' AS EXAMPLE)
%REMEMBER_2: to create properly the file excel of ephemerides use this procedure:
% 1 - Download the TLEs line for all the time step that you need using the
%     Website 'Spacetrack'.
% 2 - Then, use this TLEs to set the 'NASA horizon' website, correctly.
% 3 - Set, properly, the other parameter of 'NASA horizon' website, and run
%     the generator of ephemerides.
% 4 - Download the result of ephemerides ( in KEPLERIAN ELEMENTS ).
% 5 - Modify the .txt file downloaded, changing the ',' into ';', and then
%     the '.' into ','.
% 6 - Then, change the extention of file.txt into '.csv' or '.xlsx' (excel)
% 7 - Check the number and eliminate the useless column, then set correctly
%     the page (look the example below).
% 8 - Use this beatiful function to convert this file excel(.xlsx) into a
%     matlab matrix.


% Convert the number into the excel sheet into a matlab matrix: 
M_excel = readmatrix(namefileexcel);

for i=1:1:length(M_excel(:,1))
% Conversion from Julian Day number to Modified Julian day 2000 number:
mjd2000_date(i) = jd2mjd2000(M_excel(i,1)); %[days]

% Conversion from modified Julian day 2000 number to Gregorian calenda date (our calendar) :
% date(i,:) = mjd20002date(mjd2000_date(i)); 
end
% Conversion from mjd2000_date in [days] to mjd2000_date in [sec]:
time_real_eph = (mjd2000_date')*86400; %[sec]

% Set of usefull Keplerian parameter from the M_excel:
kep_real_eph(:,1) = M_excel(:,11); % [km] - semi-major axis
kep_real_eph(:,2) = M_excel(:,2); % [-] - eccentricity
kep_real_eph(:,3) = M_excel(:,4); % [deg] - inclination
kep_real_eph(:,4) = M_excel(:,5); % [deg] - longitude of ascending node
kep_real_eph(:,5) = M_excel(:,6); % [deg] - argoument of perifocus
kep_real_eph(:,6) = M_excel(:,10); % [deg] - true anomaly

end