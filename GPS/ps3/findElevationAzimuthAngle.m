function [E] = findElevationAzimuthAngle(receiverPos,GPSPos)

rs        = -receiverPos+GPSPos;% Receiver to Satellite vector
rs_norm   = rs/norm(rs); % Normalized rs to find unit vecotor of rs or its direction
rs_norm_p = sqrt(rs_norm(1)^2 + rs_norm(2)^2);

E = atan(rs(3)/norm(rs))/(2*pi)*360;