% clear; clc;
% %% Load Data from dfDataHead.bin
% % Use the document fftAcqTheory.pdf found on Canvas as your guide.
% % Recall that you studied the GP2015 front end in Problem Set 4. The GP2015
% % produces digitized data with an intermediate frequency
% % fIF = 1.405396825396879 MHz and a sampling rate Ns = 40e6/7 samples per
% % second. In the absence of Doppler, there would be Ns/1000 = 40000/7 ≈ 5714 samples per GPS L1 C/A code.
% %----- Setup
% Tfull = 4;                % Time interval of data to load
% fs = 40e6/7;                % Sampling frequency (Hz)
% T = 1/fs;
% N = fs*Tfull;
% N = floor(N/16)*16;         % Number of data samples to load
% nfft = 2^10;                % Size of FFT used in power spectrum estimation
% fIF  =  1.405396825396879e6; % Hz
% %----- Load data
% fid = fopen(["C:\Users\gsh04\Desktop\2024-Fall\GPS\ps5\dfDataHead.bin"], 'r','l');
% [Y,count] = binloadSamples(fid,N,'dual');
% 
% Y = Y(floor(fs*3/16)*16:end,1);
% if(count ~= N)
%     error('Insufficient data');
% end
% %% Coarse Search
% disp('----------------------------------------------------------')
% fprintf('                    Coarse Search\n')
% disp('----------------------------------------------------------')
% 
% % Coarse Search Parameter
% fDRange = [-5000:100:5000];
% prn = [1:37];
% NC = 10;                    % Noncoherent sum number
% Ta = 0.001;                 % Accumulation time in seconds
% coarsets = zeros(1, length(prn));
% coarsefD = zeros(1, length(prn));
% % Estimate
% for m = 1:length(prn)
%     [ts, fD,~] = acquisition(Y,prn(m),fDRange,NC,Ta,fs,fIF);
%     coarsets(m) = ts;
%     coarsefD(m) = fD;
% end
% 
% disp('----------------------------------------------------------')
% 
% %% Fine Search
% disp('----------------------------------------------------------')
% fprintf('                    Fine Search\n')
% disp('----------------------------------------------------------')
% 
% prnFine =  find(~isnan(coarsefD));
% Ta =0.01;
% NC = 1;
% tsFine = zeros(length(prnFine));
% fDFine = zeros(length(prnFine));
% for h = 1
%     fDmaxFine = coarsefD(prnFine(h))+20;
%     fDminFine = coarsefD(prnFine(h))-20;
%     fDRangeFine = [fDminFine:5:fDmaxFine];
%     [ts, fD, Sk2,noiseVariance] = acquisition(Y,prnFine(h),fDRangeFine,NC,Ta,fs,fIF);
%     tsFine(h) = ts;
%     fDFine(h) = fD;
%     peakSk2(h) = Sk2;
%     sigmaIQ2(h) = noiseVariance;
% end

close all;
%% Open Loop Values
fc = 1575.42*1e6;
fDOpenloop = -2208;
Tcode = 0.001/(1+fDOpenloop/fc);

s.IsqQsqAvg = peakSk2;
s.sigmaIQ = sqrt(sigmaIQ2(1));

thetaHat  = 0 ;

loopOrder = 3;
Bn_target = 10;
[s.Ad,s.Bd,s.Cd,s.Dd,Bn_act] = configureLoopFilter(Bn_target,Ta,loopOrder);
vTheta=2*pi*-fDOpenloop;
[V,D] = eig(s.Ad);
q    = vTheta./(s.Cd*V(:,1));
s.xk = q*V(:,1);
teml = 0.5;

vTheta_history(1) = vTheta;
Nk = round(Ta/T);

tshat(1) = tsFine(1);
tshat(2:round(N/Nk))= tshat(1)+[1:round(N/Nk)-1].*Tcode;

for k = 2:length(tshat)
    [Ske,Skp,Skl] = correlate(Y,fIF,tshat(k-1),vTheta,thetaHat,teml,14,T,Ta);
    s.Ip = real(Skp);
    s.Qp = imag(Skp);
    s.Ie = real(Ske);
    s.Qe = imag(Ske);
    s.Il = real(Skl);
    s.Ql = real(Skl);
    s.IsqQsqAvg = abs(Skp)^2;
    [xkp1,vTheta] = updatePll(s);
    s.vp = -vTheta/(2*pi*fc);
    s.xk = xkp1;
    vTheta_history(k) = vTheta;
    thetaHat = thetaHat+vTheta*Ta;
    s.Ta = Ta;
    s.Tc = 1e-3 / 1023; % Seconds per chip
    [vTotal] = updateDll(s);

    Sk2_history(k) = abs(Skp^2);
    Sk_history(k) = Skp;
    Skl_history(k) = Skl;
end
figure,
plot(Sk2_history)
figure,
plot(-vTheta_history/(2*pi))
figure,
plot(Sk_history)

% Plot IQ plane
figure;
plot(real(Sk_history(1:400)), imag(Sk_history(1:400)), 'o', 'MarkerSize', 5, 'LineWidth', 1.5);
hold on;
plot([-max(abs(real(Sk_history))) max(abs(real(Sk_history)))], [0 0], 'r--'); % Real axis
grid on;
axis equal;
xlabel('Re(S_{p,k})');
ylabel('Im(S_{p,k})');
title('Phasor Plot of S_{p,k}');
hold off;

