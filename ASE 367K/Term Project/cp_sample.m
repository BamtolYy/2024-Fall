function [cp, C_N_alpha, S_ref, cp_moment] = cp_sample()
% This functiion determines the center of pressure of TREL rocket, which
% dimensions are provided in the ASE 367K Flight Dynamics Term Project
%
%
% ---- Assumptions:
%       1. Dimensions start from the nose
%       2. At sea level

%sample errors
errors = 0.002*(1-2*rand([1,9]));
%aerodynamic geometry
body_diameter = 0.4 + errors(1);
nosecone_length = 1.2 + errors(2);
fin_diameter = (0.8 + errors(3))*sqrt(2);
fin_tip_chord = 0.7 + errors(4);
fin_root_chord = 0.8 + errors(5);
fin_station = 1.2 + errors(6);
body_length = 6 + errors(7);
nozzle_length = 0.8 + errors(8);
nozzle_diameter = 0.25 + errors(9);

%nosecone properties
C_N_nosecone_alpha = 2;
x_nosecone = nozzle_length+body_length+(1-0.466)*nosecone_length;
C_N_body_alpha = 0;%low angle of attack

%fin properties
N_fin = 4;
f = 1;%Interference coefficient for 3,4 fins
radius = body_diameter / 2;
span = (fin_diameter-body_diameter)/2;
C_N_fin_alpha = (1+f*(radius)./(span+radius))*((4*N_fin*(span/body_diameter).^2)...
/(1+sqrt(1+(2*span/(fin_root_chord+fin_tip_chord)).^2) ) );
x_R = (fin_root_chord-fin_tip_chord)/2;
x_fin = x_R/3*(fin_root_chord+2*fin_tip_chord)/(fin_root_chord+fin_tip_chord) +...
1/6*((fin_root_chord+fin_tip_chord)-(fin_root_chord*fin_tip_chord)/...
((fin_root_chord+fin_tip_chord)));
x_fin = fin_station+nozzle_length-x_fin;

% Bonical Boattail
% CP of Conical transition XT
Xp = 7.2;               % m; Nose to d1
LT = 0.8;               % m; d1 to d2
XT = Xp + LT/3*(1+(1-body_diameter/nozzle_diameter)/(1-(body_diameter/nozzle_diameter)^2));
rocket_length = 1.2+6+0.8;
XT = rocket_length - XT;
S1 = pi*body_diameter^2/4;
S2 = pi*nozzle_diameter^2/4;
CNaCB = 8/(pi*body_diameter^2)*(S2-S1);

%combined properties
C_N_alpha = C_N_fin_alpha + C_N_nosecone_alpha +CNaCB;
cp = (x_fin*C_N_fin_alpha + x_nosecone*C_N_nosecone_alpha+CNaCB*(XT))/C_N_alpha;
S_ref = 0.25*pi*body_diameter.^2;
cp_moment = (x_fin.^2*C_N_fin_alpha + x_nosecone.^2*C_N_nosecone_alpha+(XT).^2*CNaCB)/C_N_alpha;

