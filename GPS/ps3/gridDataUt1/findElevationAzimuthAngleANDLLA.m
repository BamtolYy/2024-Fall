function [E, A, r_lla, s_lla]=findElevationAzimuthAngleANDLLA(recPos,satPos)
%%  Calculate Elevation and Azimuth angle using LLA and Considering WGS84 Model
% INPUTS
%
% recPos --- 3x1 position vector of receiver in ECEF
%
% satPos --- 3x1 position vector of satellite in ECEF
%
% OUTPUTS
%
% E ---- Elevation Angle from receiver seeing satellite in rad
%
% A ---- Azimuth Angle from receiver seeing satellite in rad
%
% r_LLA -- [Latitude, Longitude, Alitude] Geodetic of receiver in rad
%
% s_LLA -- [Latitude, Longitude, Alitude] Geodetic of satellite in rad
%+------------------------------------------------------------------------------+
% Reference: https://www.ngs.noaa.gov/CORS/Articles/SolerEisemannJSE.pdf
% https://www.youtube.com/watch?v=dwccsh_aRiM&t=466s
%+==============================================================================+

%% Get LLA
r_lla    = ecef2lla(recPos','WGS84')/180*pi; % Also can be calculated with:
                                            % z = r*sin(lat); y = r*cos(lat)*sin(long)
s_lla    = ecef2lla(satPos','WGS84')/180*pi;
r_lat    = r_lla(1);
s_lat    = s_lla(1);
r_long   = r_lla(2);
s_long   = s_lla(2);

%% Calculate subtanding angle between receiver vector and satellite
%  subpoint vector.(ASD05 Look Angles on Youtube)

r_r   = norm(recPos); % Length of the receiver postition vector or
% typically radius of earth for receivers on
% the ground check iif it is around 6378 km
s_r   = norm(satPos); % Length of the GPS postition vector

gamma = acos(sin(s_lat)*sin(r_lat)+cos(s_lat)*cos(r_lat)...
    *cos(s_long-r_long));

%% Elevation Angle (ASD05 Look Angles on Youtube)
d = sqrt(1+(r_r/s_r)^2-2*(r_r/s_r)*cos(gamma));
E = acos(sin(gamma)/(d));

%% Azimuth Angle (ASD05 Look Angles on Youtube)


alpha = asin(sin(abs(r_long-s_long))*cos(s_lat)/sin(gamma));
if s_long <  r_long && tan(r_lat)*cos(s_long-r_long) < tan(s_lat)
    % NorthWest
    A = 2*pi - alpha;
elseif s_long < r_long && tan(r_lat)*cos(s_long-r_long) > tan(s_lat)
    % SouthWest
    A = pi + alpha;
elseif s_long > r_long && tan(r_lat)*cos(s_long-r_long) < tan(s_lat)
    % NorthEast
    A = alpha;
elseif s_long > r_long && tan(r_lat)*cos(s_long-r_long) > tan(s_lat)
    % SouthEast
    A = pi - alpha;
end








%% --------------------------------------------------------------------------
% IGNORE BELOW WORKS IT WAS FOR TRYING DIFFERENT METHODS
% Try 2
% d2 =norm (satPos-recPos);
% E2 = acos(s_r*sin(gamma)/d2)/(pi)*180;
% % Try 3
% d3 = s_r*sqrt(1+(r_r/s_r)^2-2*(r_r/s_r)*cos(gamma));
% E3 = acos(s_r*sin(gamma)/d3)/(1*pi)*180;

% Calcualte Elevation angle just using vectors:
% Reference https://www.youtube.com/watch?v=dwccsh_aRiM&t=466s for Geometry
% d5 = norm (satPos-recPos);
% psi = acos((s_r^2-d5^2-r_r^2)/(-2*d5*r_r));
% E5 = (psi-pi/2)/(pi)*180;

% Another form of alpha equation
% alpha2 = atan(tan(abs(s_long-r_long))/sin(r_lat));
% if s_long <  r_long && tan(r_lat)*cos(s_long-r_long) < tan(s_lat)
%     % NorthWest
%     A2 = 2*pi - alpha2;
% elseif s_long < r_long && tan(r_lat)*cos(s_long-r_long) > tan(s_lat)
%     % SouthWest
%     A2 = pi + alpha2;
% elseif s_long > r_long && tan(r_lat)*cos(s_long-r_long) < tan(s_lat)
%     % NorthEast
%     A2 = alpha2;
% elseif s_long > r_long && tan(r_lat)*cos(s_long-r_long) > tan(s_lat)
%     % SouthEast
%     A2 = pi - alpha2;
% end


