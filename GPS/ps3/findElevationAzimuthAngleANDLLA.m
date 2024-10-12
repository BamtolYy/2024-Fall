function [E,A,LLA]=findElevationAzimuthAngleANDLLA(recPos,satPos)
%
% INPUTS
%
% recPos --- 3x1 position vector of receiver in ECEF
% 
% satPos --- 3x1 position vector of satellite in ECEF
%
% OUTPUTS
% 
% E ---- Elevation Angle from receiver seeing satellite
% 
% A ---- Azimuth Angle from receiver seeing satellite
% 
% r_LLA -- [Latitude, Longitude, Alitude] Geodetic of receiver
% 
% s_LLA -- [Latitude, Longitude, Alitude] Geodetic of satellite
%+------------------------------------------------------------------------------+
% Reference: https://www.ngs.noaa.gov/CORS/Articles/SolerEisemannJSE.pdf
% https://www.youtube.com/watch?v=dwccsh_aRiM&t=466s
%+==============================================================================+
%% Calculate it using LLA and Considering WGS84 Model
r_lla    = ecef2lla(recPos,'WGS84')/180*pi;  % Also can be calculated with:
                                                % z = r*sin(lat); y = r*cos(lat)*sin(long)
s_lla    = ecef2lla(satPos,'WGS84')/180*pi;

r_r = norm(receiverPos); % Length of the receiver postition vector or
                              % typically radius of earth for receivers on
                              % the ground check iif it is around 6378 km

t_r      = norm(GPSPos);      % Length of the GPS postition vector
r_lat    = r_lla(1); 
s_lat    = s_lla(1);
r_long   = r_lla(2);
s_long   = s_lla(2);
% Calculate subtanding angle between reciver vector and satellite subpoint
% vector. (ASD05 Look Angles on Youtube)
gamma  = acos(sin(t_lat)*sin(r_lat)+cos(t_lat)*cos(r_lat)...
    *cos(t_long-r_long));

% Elevation Angle (ASD05 Look Angles on Youtube)
d1 = sqrt(1+(r_r/t_r)^2-2*(r_r/t_r)*cos(gamma));
E = acos(sin(gamma)/(d1))/(2*pi)*360;

% Azimuth Angle (ASD05 Look Angles on Youtube)
alpha = atan(tan(abs(t_long-r_long))/sin(r_lat));
A = asin(sin(abs(r_lla(2)-t_lla(2)))*cos(t_lat)/sin(gamma));


% Try 2
d2 =norm (t_ecef-r_ecef); 
E2 = acos(t_r*sin(gamma)/d2)/(2*pi)*360;

% Try 3

d3 = t_r*sqrt(1+(r_r/t_r)^2-2*(r_r/t_r)*cos(gamma));
E3 = acos(t_r*sin(gamma)/d3)/(2*pi)*360;


%--------------------------------------------------------------------------



% Calcualte Elevation angle just using vectors:
% Reference https://www.youtube.com/watch?v=dwccsh_aRiM&t=466s for Geometry
r_r    = norm(receiverPos); % Length of the receiver postition vector or
                            % typically radius of earth for receivers on
                            % the ground check iif it is around 6378 km
t_r    = norm(GPSPos);      % Length of the GPS postition vector                            
d = norm (t_ecef-r_ecef); 
psi = acos((t_r^2-d^2-r_r^2)/(-2*d*r_r));
E = (psi-pi/2)/(2*pi)*360;


