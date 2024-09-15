% topSimulateTrainDoppler.m
%
% Top-level script for train Doppler simulation


%clear; clc; clf;
%----- Setup
fc = 440;
vTrain = 20;
t0 = 0;
x0 = 0;
delt = 0.01;
N = 1000;
vs = 343;
xObs = 50;
dObs = -60;

%----- Simulate
[fDVec,tVec] = simulateTrainDoppler(fc,vTrain,t0,x0,xObs,dObs,delt,N,vs);
fApparentVec = fDVec + fc;

%----- Plot
plot(tVec,fDVec + fc, 'r');
xlabel('Time (seconds)');
ylabel('Apparent horn frequency (Hz)');
grid on;
shg;

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


%----- Read provided audio file
[y, fs] = audioread('H.wav');

%----- Frequency Domain

f    = linspace(0,fs,fs);
Y    = abs(fft(soundVec));
Z    = abs(fft(y));



figure(1);
plot(f(1:fs/2),Y(1:fs/2), f(1:fs/2),Z(1:fs/2)),
legend('Guess','Original')


% ----- Play the sound vector
sound([soundVec,y], fs);   