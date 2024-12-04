close all; clear all; clc;
%% a)
Ta = 10/1000;       % seconds
Bn = 10;            % Hz
%% b)
%---- Define Loop Filter
loopOrder = 3;
[s.Ad,s.Bd,s.Cd,s.Dd,Bn_act] = configureLoopFilter(Bn,Ta,loopOrder);
Bn_act
%---- Generate Ficticious Phase Time History
fs = 10000;                 % Hz; Ficticious Sampling rate of signal
T  = 1/fs;                  % seconds; Sampling Interval
t = (0:T:60)';
PhaseHist = 0.001*t;         % Raw Phase History
% Average over subinterval Ta
N_Ta    = floor(Ta/T);        % Number of samples in an accumulation
N_sub   = floor(length(t)/N_Ta);     % Number of subintervals in the signal
xk0     = zeros(loopOrder-1,1);
s.xk    = xk0; 
thetahat = zeros(N_sub,1);
for ii = 1:N_sub
    start = (ii-1)*N_Ta+1;
    endpt = ii*N_Ta;
    phase_segment = PhaseHist(start:endpt);
    avg = mean(phase_segment);
    s.Ip =cos(avg);
    s.Qp =sin(avg);
    % Update PLL
    [xkp1,vk] = updatePll(s);
    s.xk = xkp1;
    vkk(ii)=vk;
    thetahat(ii+1) = thetahat(ii)+vk*Ta;
end
t2 = (0:Ta:60)';
plot(t,PhaseHist,t2,thetahat)
legend('True','Estimate')





