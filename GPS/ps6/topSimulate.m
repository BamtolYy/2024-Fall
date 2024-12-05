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
t = (0:T:20)';
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

PhaseHist2 = 2*pi*fd*t;         % Raw Phase History
a= round(length(PhaseHist2)/2);
PhaseHist2(a:end) = PhaseHist2(a:end)+pi;
s.xk    = xk0;
thetahat2 = zeros(length(t),1);
deltaThetak2 = 0;
s.Ip = cos(deltaThetak2);
s.Qp = sin(deltaThetak2);
for ii = 1:N_sub
    start = (ii-1)*Nk+1;
    endpt = ii*Nk;
    [xkp1,vk] = updatePll(s);
    s.xk=xkp1;
    if ii == 1
         thetahat2(start:endpt)=thetahat2(start)+vk*[0:T:Ta-T];
    else
    thetahat2(start:endpt)=thetahat2(start-1)+vk*[0:T:Ta-T];
    deltaTheta2 = PhaseHist2(start:endpt)-thetahat2(start:endpt);
    deltaThetak2 = mean(deltaTheta2);
    s.Ip = cos(deltaThetak2);
    s.Qp = sin(deltaThetak2);
    end
end
figure,
plot(t,PhaseHist2,t(1:end-1),thetahat2(1:end-1),'--')
legend('True','Estimate')
xlabel('Time (seconds')
ylabel('\theta')
title('Phase Detection with 180 Phase Transition without Noise')

fprintf(['When there is a 180 degrees phase transition, the Loop\n' ...
    'overestimates the phase for a small amount of time. Just like\n' ...
    'any control system there is the there is a delay in error calculation,\n' ...
    'so the loop overestimates the error for an instance of time and causes\n' ...
    'the overshoot'])

%% d)
sigman=0.9;
PhaseHist3 = 2*pi*fd*t;         % Raw Phase History
% Average over subinterval Ta
Nk    = floor(Ta/T);        % Number of samples in an accumulation
N_sub   = floor(length(t)/Nk);     % Number of subintervals in the signal
xk0     = zeros(loopOrder-1,1);
s.xk    = xk0;
thetahat3 = zeros(length(t),1);
deltaThetak3 = 0;
s.Ip = cos(deltaThetak3);
s.Qp = sin(deltaThetak3);
for ii = 1:N_sub
    start = (ii-1)*Nk+1;
    endpt = ii*Nk;
    [xkp1,vk] = updatePll(s);
    s.xk=xkp1;
    if ii == 1
         thetahat3(start:endpt)=thetahat3(start)+vk*[0:T:Ta-T];
    else
    thetahat3(start:endpt)=thetahat3(start-1)+vk*[0:T:Ta-T];
    deltaTheta3 = PhaseHist3(start:endpt)-thetahat3(start:endpt);
    deltaThetak3 = mean(deltaTheta3);
    s.Ip = cos(deltaThetak3)+sigman*randn;
    s.Qp = sin(deltaThetak3)+sigman*randn;
    end
end
figure,
plot(t,PhaseHist3,t(1:end-1),thetahat3(1:end-1),'--')
legend('True','Estimate')
xlabel('Time (seconds)')
ylabel('\theta')
title('Phase Detection with Noise')
CN0 = 
