% topSimulateTrainDoppler.m
% Top-level script for train Doppler simulation

clear; clc; clf;
%----- Setup
fc = 440;
vTrain = 20;
t0 = 0;
x0 = 0;
delt = 0.01;
N = 1000;
vs = 343;
xObs = 73.6;
dObs = 7.1;

%----- Simulate
[fDVec,tVec] = simulateTrainDoppler(fc,vTrain,t0,x0,xObs,dObs,delt,N,vs);
fApparentVec = fDVec + fc;

%----- Plot
% figure(1)
% plot(tVec,fDVec + fc, 'r');
% xlabel('Time (seconds)');
% ylabel('Apparent horn frequency (Hz)');
% grid on;
% shg;

%----- Generate a sound vector
T = delt*N;                    % simulation time (sec)
fs = 22050;                    % sample frequency (Hz)
deltSamp = 1/fs;               % sampling interval (sec)
Ns = floor(T/deltSamp);        % number of samples
tsamphist = [0:Ns-1]'*deltSamp;
Phihist = zeros(Ns,1);
fApparentVecInterp = interp1(tVec,fApparentVec,tsamphist,'spline');
for ii=2:Ns
  fii = fApparentVecInterp(ii);
  Phihist(ii) = Phihist(ii-1) + 2*pi*fii*deltSamp;
end
soundVec = sin(Phihist);
%  

%----- Write to audio file
audiowrite('trainout.wav',soundVec,fs);

%----- Write frequency time history to output file
save trainData fApparentVec tVec

%% Import Original File
%----- Read provided audio file
[y, fs] = audioread('originalSoundWave.wav');

%----- Frequency Domain

OPhihist = asin(y);
for ii=2:Ns
  Ofii (ii)=(OPhihist(ii)-OPhihist(ii-1))/(2*pi*deltSamp);
  original(ii) = abs(Ofii(ii));
  
end

OtVec = linspace(t0,tVec(end),length(y));
plotoriginal = interp1(OtVec,original,tVec,'spline');
figure(2)
plot(tVec,plotoriginal,'o',tVec,fDVec + fc, 'r')
xlabel('Time [s]')
ylabel('Doppler Frequency(Hz)')
ylim([410 470]);
title('Doppler Frequency Change over Time with observer postion estimate at [73.6,7.1]')
legend('trainout.wav','Estimated')

% ----- Play the sound vector
sound([soundVec,y], fs);   

% spectrogram(y,2^10,[],[],100,'yaxis')
% hold on,
% % ylim([0.5 0.4])
% % plot(tVec,fDVec + fc, 'r');



