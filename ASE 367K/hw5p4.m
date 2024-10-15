syms hx hy hz p q r pd qd rd Ixx Iyy Izz Ixz L M N


IB         = [Ixx 0 -Ixz; 0 Iyy 0; -Ixz 0 Izz];
omegaBdot  = [pd;qd;rd];
omegatilde = [0 -r q; r 0 -p; -q p 0];
omegaB     = [p;q;r];
hB        = [hx; hy; hz];

MB = IB*omegaBdot + omegatilde*(IB*omegaB+hB)