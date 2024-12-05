close all; clear all; clc;
%% a)
Ta = 10/1000;       % seconds
Bn = 10;            % Hz
%% b)
%---- Define Loop Filter
loopOrder = 3;
[s.Ad,s.Bd,s.Cd,s.Dd,Bn_act] = configureLoopFilter(Bn,Ta,loopOrder);
%---- Generate Ficticious Phase Time History
fs = 10000;                 % Hz; Ficticious Sampling rate of signal
T  = 1/fs;                  % seconds; Sampling Interval
t = (0:T:5)';
fd = 1;
PhaseHist = 2*pi*fd*t;         % Raw Phase History

% Average over subinterval Ta
Nk    = floor(Ta/T);        % Number of samples in an accumulation
N_sub   = floor(length(t)/Nk);     % Number of subintervals in the signal
xk0     = zeros(loopOrder-1,1);
s.xk    = xk0;
thetahat = zeros(length(t),1);
deltaThetak = 0;
s.Ip = cos(deltaThetak);
s.Qp = sin(deltaThetak);
for ii = 1:N_sub
    start = (ii-1)*Nk+1;
    endpt = ii*Nk;
    [xkp1,vk] = updatePll(s);
    s.xk=xkp1;
    if ii == 1
         thetahat(start:endpt)=thetahat(start)+vk*[0:T:Ta-T];
    else
    thetahat(start:endpt)=thetahat(start-1)+vk*[0:T:Ta-T];
    deltaTheta = PhaseHist(start:endpt)-thetahat(start:endpt);
    deltaThetak = mean(deltaTheta);
    s.Ip = cos(deltaThetak);
    s.Qp = sin(deltaThetak);
    end
end
figure,
plot(t,PhaseHist,t(1:end-1),thetahat(1:end-1),'--')
legend('True','Estimate')
xlabel('Time (seconds')
ylabel('\theta')
title('Phase Detection without Noise')



%% c)
clear all; close all; clc;
PhaseHist = 2*pi*fd*t;         % Raw Phase History


