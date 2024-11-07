clear all;close all;clc;


F_N_static = 216000;









const_D = 0.5*rho*S*C_D_0;
const_L = 0.5*rho*S*C_L_0;

W0 = 790100;






for ii=1:Num-1
    m_Vdot = const_F*(k_prime_F0 + k_prime_F1*V_Mat(ii) + k_prime_F2*V_Mat(ii)^2)...
        -const_D*V_Mat(ii)^2 - mu*((W0-Wf_dot*t_k)-const_L*V_Mat(ii)^2);
    Vdot   = m_Vdot*g/(W0-Wf_dot*t_k);
    V_Mat(ii+1) = Vdot*dt + V_Mat(ii);
    tk = t_k+dt;
end