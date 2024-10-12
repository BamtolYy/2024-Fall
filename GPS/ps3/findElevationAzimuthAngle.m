function [E,A]=findElevationAzimuthAngle(receiverPos,GPSPos)

r_ecef = receiverPos;       % Receiver postion in ECEF
t_ecef = GPSPos;            % GPS or Transmitter position in ECEF
r_lla  = ecef2lla(r_ecef,'WGS84')/360*2*pi;  % Also can be calculated with:
                            % z = r*sin(lat); y = r*cos(lat)*sin(long)
t_lla  = ecef2lla(t_ecef,'WGS84')/360*2*pi;
r_r    = norm(receiverPos); % Length of the receiver postition vector or
                            % typically radius of earth for receivers on
                            % the ground check iif it is around 6378 km

t_r    = norm(GPSPos);      % Length of the GPS postition vector
t_lat  = t_lla(1); 
r_lat  = r_lla(1);
t_long = t_lla(2);
r_long = r_lla(2);

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
% Reference: https://www.ngs.noaa.gov/CORS/Articles/SolerEisemannJSE.pdf
% https://www.youtube.com/watch?v=dwccsh_aRiM&t=466s


% Try 4

d3 = norm (t_ecef-r_ecef); 

psi = acos((t_r^2-d3^2-r_r^2)/(-2*d*r_r));
E4 = (psi-pi/2)/(2*pi)*360;


