% clear; clc;
% %% Load Data from dfDataHead.bin
% % Use the document fftAcqTheory.pdf found on Canvas as your guide.
% % Recall that you studied the GP2015 front end in Problem Set 4. The GP2015
% % produces digitized data with an intermediate frequency
% % fIF = 1.405396825396879 MHz and a sampling rate Ns = 40e6/7 samples per
% % second. In the absence of Doppler, there would be Ns/1000 = 40000/7 ≈ 5714 samples per GPS L1 C/A code.
% %----- Setup
% Tfull = 70;                % Time interval of data to load
% fs = 40e6/7;                % Sampling frequency (Hz)
% T = 1/fs;
% N = floor(fs*Tfull);
% N = floor(N/16)*16;         % Number of data samples to load
% nfft = 2^10;                % Size of FFT used in power spectrum estimation
% fIF  =  1.405396825396879e6; % Hz
% %----- Load data
% % fid = fopen(["C:\Users\gsh04\Desktop\2024-Fall\GPS\ps5\dfDataHead.bin"], 'r','l');
% fid = fopen(["C:\Users\gsh04\Desktop\2024-Fall\GPS\ps7\dfDataHead.bin"], 'r','l');
% [Z,count] = binloadSamples(fid,N,'dual');
% 
% timeTotal = [0:T:N*T-T]';
% acquisitionStartTime = 3;
% % Y = Z(round(fs*acquisitionStartTime)+1:end,1);
% Y = Z(:,1);
% 
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
% NC = 2;                    % Noncoherent sum number
% Ta = 0.01;                 % Accumulation time in seconds
% coarsets = zeros(1, length(prn));
% coarsefD = zeros(1, length(prn));
% % Estimate
% for m = 1:length(prn)
%     [ts, fD,~,~,~] = acquisition(Y,prn(m),fDRange,NC,Ta,fs,fIF);
%     coarsets(m) = ts;
%     coarsefD(m) = fD;
% end
% 
% disp('----------------------------------------------------------')

%% Fine Search
disp('----------------------------------------------------------')
fprintf('                    Fine Search\n')
disp('----------------------------------------------------------')

prnFine =  find(~isnan(coarsefD));
Ta =0.01;
NC = 6;
tsFine = zeros(length(prnFine));
fDFine = zeros(length(prnFine));
for h = 4
    fDmaxFine = coarsefD(prnFine(h))+100;
    fDminFine = coarsefD(prnFine(h))-100;
    fDRangeFine = [fDminFine:1 :fDmaxFine];
    [ts, fD, Sk2,noiseVariance,sk] = acquisition(Y,prnFine(h),fDRangeFine,NC,Ta,fs,fIF);
    tsFine(h) = ts;
    fDFine(h) = fD;
    peakSk2(h) = Sk2;
    sigmaIQ2(h) = noiseVariance;
end

%% (c) Initialize the beat carrier phase estimate
thetaHat = 0;
thetaHat_history(1) = thetaHat;

%% (d) Initialize Moving Window Average
g = h;                      % PRN index
s.IsqQsqAvg = peakSk2(g);
movingAverage = peakSk2(g);
s.Ip = real(sk);
s.Qp = imag(sk);
s.sigmaIQ = sqrt(sigmaIQ2(g));
%% (e)
loopOrder = 3;
Bn_target = 10;
s.Ta = Ta;
s.Tc = 1e-3/1023;             % Chip interval in seconds
[s.Ad,s.Bd,s.Cd,s.Dd,Bn_act] = configureLoopFilter(Bn_target,s.Ta,loopOrder);


%% (f) x_k=0 calculation
vTheta=2*pi*-fDFine(g);
% vTheta=2*pi*2202;

% vTheta=2*pi*3341.3;


[V,D] = eig(s.Ad);
q    = vTheta/(s.Cd*V(:,1));
s.xk = q*V(:,1);

%% (g) Correlate
teml = 0.5;                 % Chips
fc = 1575.42*1e6;
tstart = tsFine(g);
% s.sigmaIQ = sqrt(1.05e5);
Nk = floor(Ta/T);
NumberofAccumulation = floor(length(Y(floor(tsFine(h)/T):end))/Nk);
% vTheta_history = zeros(NumberofAccumulation-1,1);
vTheta_history(1) = vTheta;
Sk2_history    = zeros(NumberofAccumulation-1,1);

for k = 1 : NumberofAccumulation-1
    [Se_k, Sp_k, Sl_k] = performCorrelations(Y, fs, fIF, tstart, vTheta, thetaHat, teml, prnFine(g), Ta);
    %% (h) Update Moving Window AVerage
    s.Ip = real(Sp_k);
    s.Qp = imag(Sp_k);
    s.Ie = real(Se_k);
    s.Qe = imag(Se_k);
    s.Il = real(Sl_k);
    s.Ql = imag(Sl_k);
    if k <100
        movingAverage(k+1) = abs(Sp_k)^2;
        s.IsqQsqAvg = mean(movingAverage);
    else
        movingAverage = circshift(movingAverage,-1);
        movingAverage(end) = abs(Sp_k)^2;
        s.IsqQsqAvg = mean(movingAverage);
    end
    %% (i)
    fc = 1575.42*1e6;
    [s.xk,vTheta] = updatePll(s);
    s.vp = -vTheta/(2*pi*fc);
    [vTotal] = updateDll(s);
    vTotal_hist(k+1) = vTotal;
    %% (j) Update Beat Carrier Phase Estimate
    thetaHat = thetaHat+vTheta*Ta;
    vTheta_history(k+1) = vTheta;
    %% (k)
    tstart = tstart -vTotal*Ta+Ta;
    Sk2_history(k+1) = abs(Sp_k^2);
    Sk_history(k+1) = Sp_k;
    Skl_history(k+1) = abs(Sp_k)^2-abs(Sl_k)^2;
    thetaHat_history (k+1) = thetaHat;
    tstart_history(k+1) = tstart;
    % % Plot IQ plane
    % figure;
    % plot(real(Sk_history(k)), imag(Sk_history(k)), 'o', 'MarkerSize', 5, 'LineWidth', 1.5);
    % hold on;
    % plot([-max(abs(real(Sk_history(k)))) max(abs(real(Sk_history(k))))], [0 0], 'r--'); % Real axis
    % grid on;
    % axis equal;
    % xlabel('Re(S_{p,k})');
    % ylabel('Im(S_{p,k})');
    % title('Phasor Plot of S_{p,k}');

end
% figure,
% plot(Sk2_history)
% title('Sk2')
figure,
Time = [acquisitionStartTime+tsFine(g):Ta:acquisitionStartTime+length(vTheta_history)*Ta]';
plot(Time,-vTheta_history/(2*pi))
title(['fD for PRN: ', num2str(prnFine(g))]);
% ylim([-2235 -2205])
xlim([0 80])
ylabel ('Hz')
xlabel('Time (sec)')

grid on
% figure,
% plot(Skl_history)
% title('Skp-Skl')

% % Plot IQ plane
% figure;
% plot(real(Sk_history), imag(Sk_history), 'o', 'MarkerSize', 5, 'LineWidth', 1.5);
% hold on;
% plot([-max(abs(real(Sk_history))) max(abs(real(Sk_history)))], [0 0], 'r--'); % Real axis
% grid on;
% axis equal;
% xlabel('Re(S_{p,k})');
% ylabel('Im(S_{p,k})');
% title('Phasor Plot of S_{p,k}');
% hold off;