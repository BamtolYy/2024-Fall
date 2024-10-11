function [E,A]=findElevationAzimuthAngle(receiverPos,GPSPos)

r_ecef = receiverPos;       % Receiver postion in ECEF
t_ecef = GPSPos;            % GPS or Transmitter position in ECEF
r_lla  = ecef2lla(r_ecef);  % Also can be calculated with:
                            % z = r*sin(lat); y = r*cos(lat)*sin(long)
t_lla  = ecef2lla(t_ecef);
r_r    = norm(receiverPos); % Length of the receiver postition vector or
                            % typically radius of earth for receivers on
                            % the ground check iif it is around 6378 km

t_r    = norm(GPSPos);      % Length of the GPS postition vector
 

% Calculate subtanding angle between reciver vector and satellite subpoint
% vector. (ASD05 Look Angles on Youtube)
gamma  = acos(sin(t_lla(1))*sin(r_lla(1))+cos(t_lla(1))*cos(r_lla(1))...
    *cos(t_lla(2)-r_lla(2)));

% Elevation Angle (ASD05 Look Angles on Youtube)
E = acos(sin(gamma)/sqrt(1+(r_r/t_r)^2-2*(r_r/t_r)*cos(gamma)))/(2*pi)*360;

% Azimuth Angle (ASD05 Look Angles on Youtube)
A = asin(sin(abs(r_lla(2)-t_lla(2)))*cos(t_lla(1))/sin(gamma))/(2*pi)*360;
    


