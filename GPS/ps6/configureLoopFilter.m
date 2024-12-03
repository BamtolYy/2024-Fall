function [Ad,Bd,Cd,Dd,Bn_act] = configureLoopFilter(Bn_target,Ta,loopOrder)
% configureLoopFilter : Configure a discrete-time loop filter for a feedback
%                       tracking loop.
%
%
% INPUTS
%
% Bn_target ----- Target loop noise bandwidth of the closed-loop system, in
%                 Hz.
%
% Ta ------------ Accumulation interval, in seconds. This is also the loop
%                 update (discretization) interval.
% loopOrder ----- The order of the closed-loop system. Possible choices
%                 are 1, 2, or 3.
%
%
% OUTPUTS
%
% Ad,Bd,Cd,Dd --- Discrete-time state-space model of the loop filter.
%
% Bn_act -------- The actual loop noise bandwidth (in Hz) of the closed-loop
%                 tracking loop as determined by taking into account the
%                 discretized loop filter, the implicit integration of the
%                 carrier phase estimate, and the length of the accumulation
%                 interval.
%
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+

Az   = tf([1 1],[2 0],Ta);
NCO  = tf(1,[1 0]);
NCOd = c2d(NCO,Ta);
switch loopOrder
    case 1
        k = 4*Bn_target;
        Ds = tf(k,1);
        Dz = c2d(Ds,Ta);
        Hz = feedback(Az*NCOd*Dz,1);
    case 2
        k = 8/3*Bn_target;
        a= k/2;
        Ds = tf(k*[1 a],[1 0]);
        Dz = c2d(Ds,Ta);
        Hz = feedback(Az*NCOd*Dz,1);
    case 3
        a = 1.2*Bn_target;
        b = a^2/2;
        k = 2*a;
        Ds = tf(k*[1 a b],[1 0 0]);
        Dz = c2d(Ds,Ta);
        Hz = feedback(Az*NCOd*Dz,1);
end


% Calculate Bn_act
walias = pi/Ta;
wvec = [0:10000]'*(walias/10000);
[magvec,~] = bode(Hz,wvec);
magvec = magvec(:);
Bn_act = sum(magvec.^2)*mean(diff(wvec))/(2*pi*(magvec(1,1)^2));

% Convert the loop filter to a discrete-time state-space model
[Ad,Bd,Cd,Dd]=ssdata(Dz);
