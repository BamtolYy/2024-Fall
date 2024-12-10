function [vTotal] = updateDll(s)
% updateDll : Perform a single update step of a carrier-aided first-order code
%             tracking loop.
%
% INPUTS
%
% s ------------- A structure with the following fields:
%   
%       Bn_target -- Target bandwidth of the closed-loop code tracking loop, in
%                    Hz.
%
%       IsqQsqAvg -- The average value of |Sk|^2 = Ik^2 + Qk^2; also equal to
%                    (Nk*Abark/2)^2 + 2*sigmaIQ^2.
%
%       sigmaIQ ---- The standard deviation of the noise in the prompt in-phase
%                    and quadrature accumulations.
%
%       Ip --------- The in-phase prompt accumulation over the interval from
%                    tkm1 to tk.
%
%       Qp --------- The quadrature prompt accumulation over the interval from
%                    tkm1 to tk.
%
%       Ie --------- The in-phase early accumulation over the interval from
%                    tkm1 to tk.
%
%       Qe --------- The quadrature early accumulation over the interval from
%                    tkm1 to tk.
%
%       Il --------- The in-phase late accumulation over the interval from
%                    tkm1 to tk.
%
%       Ql --------- The quadrature late accumulation over the interval from tkm1
%                    to tk.
%
%       vp --------- The aiding signal from the phase tracking loop, in seconds
%                    per second. This is equal to the Doppler frequency shift
%                    that will be used to drive the receiver s carrier-tracking
%                    numerically controlled oscillator during the time interval
%                    from tk to tkp1, divided by the carrier frequency and
%                    multiplied by -1 (high-side mixing) or 1 (low-side mixing)
%                    depending on whether the RF front-end peforms high- or
%                    low-side mixing. Thus, vp = sMix*fDk/fc, with sMix = +/- 1.
%
%       Tc --------- Spreading code chip interval, in seconds.
%
% OUTPUTS
%       vTotal ------ The code tracking loop’s estimate of the code phase rate at
%                     sample time tk, in sec/sec. vTotal is equal to the code
%                     tracking loop s correction term v plus the carrier aiding
%                     term vp.
%
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+

% Configure Loop
loopOrder = 1;
[Ad,Bd,Cd,Dd,Bn_act] = configureLoopFilter(0.1,s.Ta,loopOrder);

% Determine Error
C = s.Tc/2*1/(s.IsqQsqAvg-2*s.sigmaIQ^2);
e = C*((s.Ie-s.Il)*s.Ip+(s.Qe-s.Ql)*s.Qp);
if isempty(Cd) 
    Cd =[0 0];
end

if isempty(Dd) 
    Dd =0;
end

vk     = Cd*s.xk + Dd*e;
vTotal = vk+s.vp;
