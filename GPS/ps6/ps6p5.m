Bn = 10;    % Hz

% 1st Order
k1 = 4*Bn;
H1 = tf(k1,[1 k1]);
% 2nd Order
k2 = 8/3*Bn;
a2 = k2/2;
H2 = tf([k2 a2*k2],[1 k2 k2*a2]);
% 3rd Order
a3 = 1.2*Bn;
b3 = a3^2/2;
k3 = 2*a3;
H3 = tf([k3 k3*a3 k3*b3],[1 k3 k3*a3 k3*b3]);
