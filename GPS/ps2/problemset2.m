clear all
clc
%%% Problem set 2

%% 1a
clear all;
close all;
clc;

syms f
min = -2;
max = 2;
z=double(int(sinc(f)^2,min,max));

% delt =0.01;
% f = min:delt:max;
% Px = 0;
% for i = 1:length(f)
%    Px = Px+delt*sinc(f(i))^2;
% end

delt =0.001;
f = min:delt:max;
Px = 0;
for i = 2:length(f)
    if isnan(sin(pi*f(i))/(pi*f(i)))
        Px = Px+delt;
        graph(i)=Px;
    else
        Px = Px +delt*(sin(pi*f(i))/(pi*f(i)))^2;
        graph(i)=Px;
    end
end
fprintf('\n');
fprintf('%.3f%%',Px*100)
fprintf('\n');

%% 1 b

Rx = sinc(f);
plot(f,Rx), hold on,
title('Autocorrelation, Rx'),
xlabel('Frequency Hz')

%% Check 2
% b)
W = 2*10^6; % Arbitrary Bandwidth
t=-3/(W):1*10^-9:3/(W);
y=sinc(W*t);
plot(t,y)



%% 6
close all;
clear all;
% Load Big Dish Data
load(['C:\Users\gsh04\Desktop\2024-Fall\GPS\ps2\ps2MatlabFiles\' ...
    'prn31_22apr03_01hrs40min00sec_gmt_fl1_46_08mhz_250msec.mat']);

% Define the number of samples for the window
num_samples = 400;

% Find 400 samples that makes a square
% Loop through different windows of Y to find a well-aligned square
% for start_idx = 1:1000:length(Y)-num_samples
%     % Extract 400 samples starting at start_idx
%     Y_window = Y(start_idx:start_idx + num_samples - 1);
%
%     % Plot the real and imaginary components in the complex plane
%     figure;
%     plot(real(Y_window), imag(Y_window), 'o');
%     xlabel('Real(Y)');
%     ylabel('Imaginary(Y)');
%     title(sprintf('Window starting at sample %d', start_idx));
%     grid on;
%
%     pause(0.5); % Pause to inspect each plot
% end

% Plot 400 samples of Real and Imag part separately respect to time domain
% figure(1)
Y_window = Y(520001:520001+400-1);
% plot(real(Y_window),imag(Y_window),'o');
% title(sprintf('Window starting at sample %d', 52001));
% grid on;
% figure(2)
% subplot(2,1,1)
% plot(real(Y_window));
% title('Real Component');
% grid on;
% subplot(2,1,2)
% plot(imag(Y_window));
% title('img Component');

% Create Time vecotr to plot the 400 samples and find chip rate, Tc
t= 0:(46.08e6)^-1:400*(46.08e6)^-1-(46.08e6)^-1;
figure(2)
plot(t,real(Y_window));
title('Real Component');

% Real (P(Y) for the range of data I chose)
% Find Average Peaks and Std of them
[realpks] = findpeaks(real(Y_window));
realamplitude=mean(abs(realpks));
realstd=std(abs(realpks));
% Find C/No Ratio using amp and std
T = (46.08e6)^-1; % Sampling Rate
realCNRatio=10*log10(realamplitude^2/(4*realstd^2*T));

% Imag (C/A for the range of data I chose)
% Find Average Peaks and Std of them
[imagpks] = findpeaks(imag(Y_window));
imagamplitude=mean(abs(imagpks));
imagstd=std(abs(imagpks));
% Find C/No Ratio using amp and std
T = (46.08e6)^-1; % Sampling Rate
imagCNRatio=10*log10(imagamplitude^2/(4*imagstd^2*T));


% Compare the amplitude of C/A to P(Y)
ampratio= imagamplitude/realamplitude;

% figure(3)
% plot(t,imag(Y_window));
% toitle('Imag Component');

%% 7
close all;
clear all;

% Generate binary sequence with ±1 values
code = sign(randn([1, 2^14]));

% Oversampling setup
M = 10.1; % Oversampling factor
delChip = 1/M; % Oversampling interval
delOffset = 0; % No offset
Np = length(code); % Length of binary sequence
Ns = round(M * Np); % Number of samples after oversampling

% Oversample the binary sequence
codeOS = oversampleSpreadingCode(code, delChip, delOffset, Ns, Np);

% Power spectrum estimation with pwelch
[ospx, osf] = pwelch(codeOS, hann(2^10), [], 2^10, 1/delChip); % one-sided result
% Convert power spectrum to dB scale
ospxdb = 10*log10(ospx);
% Plot onesided result first
figure(1)
plot(osf, ospxdb)
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
title('Power Spectral Density of Oversampled Binary Sequence')
% Zoom in on the main lobe and first few sidelobes
xlim([0, 3]) % Adjust based on the main lobe location
hold on,
%Plot twosided result in the same figure
[ospx, osf] = pwelch(codeOS, hann(2^10), [], 2^10, 1/delChip,"twosided"); % two-sided result
ospxdb = 10*log10(ospx);
plot(osf, ospxdb)
legend('onesided', 'twosided')










clear all
close all
% Change Code Length
% Generate binary sequence with ±1 values
code = sign(randn([1, 2^9]));


% Oversampling setup
M = 10.1; % Oversampling factor
delChip = 1/M; % Oversampling interval
delOffset = 0; % No offset
Np = length(code); % Length of binary sequence
Ns = round(M * Np); % Number of samples after oversampling
% Oversample the binary sequence
codeOS = oversampleSpreadingCode(code, delChip, delOffset, Ns, Np);
% Power spectrum estimation with pwelch
[ospx, osf] = pwelch(codeOS, hann(2^10), [], 2^10, 1/delChip); % one-sided result
% Convert power spectrum to dB scale
ospxdb = 10*log10(ospx);
% Plot onesided result first
figure(1)
plot(osf, ospxdb)
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
title('Power Spectral Density of Oversampled Binary Sequence with Different Code Length'),
hold on,

% Generate binary sequence with ±1 values
code = sign(randn([1, 2^14]));


% Oversampling setup
M = 10.1; % Oversampling factor
delChip = 1/M; % Oversampling interval
delOffset = 0; % No offset
Np = length(code); % Length of binary sequence
Ns = round(M * Np); % Number of samples after oversampling
% Oversample the binary sequence
codeOS = oversampleSpreadingCode(code, delChip, delOffset, Ns, Np);
% Power spectrum estimation with pwelch
[ospx, osf] = pwelch(codeOS, hann(2^10), [], 2^10, 1/delChip); % one-sided result
% Convert power spectrum to dB scale
ospxdb = 10*log10(ospx);
% Plot onesided result first
plot(osf, ospxdb)
hold on,

% Generate binary sequence with ±1 values
code = sign(randn([1, 2^20]));

% Oversampling setup
M = 10.1; % Oversampling factor
delChip = 1/M; % Oversampling interval
delOffset = 0; % No offset
Np = length(code); % Length of binary sequence
Ns = round(M * Np); % Number of samples after oversampling
% Oversample the binary sequence
codeOS = oversampleSpreadingCode(code, delChip, delOffset, Ns, Np);
% Power spectrum estimation with pwelch
[ospx, osf] = pwelch(codeOS, hann(2^10), [], 2^10, 1/delChip); % one-sided result
% Convert power spectrum to dB scale
ospxdb = 10*log10(ospx);
% Plot onesided result first
plot(osf, ospxdb);
legend('2^9','2^{14}','2^{20}');





clear all
close all
% Change the oversampling interval
% Generate binary sequence with ±1 values
code = sign(randn([1, 2^14]));


% Oversampling setup
M = 5.1; % Oversampling factor
delChip = 1/M; % Oversampling interval
delOffset = 0; % No offset
Np = length(code); % Length of binary sequence
Ns = round(M * Np); % Number of samples after oversampling
% Oversample the binary sequence
codeOS = oversampleSpreadingCode(code, delChip, delOffset, Ns, Np);
% Power spectrum estimation with pwelch
[ospx, osf] = pwelch(codeOS, hann(2^10), [], 2^10, 1/delChip); % one-sided result
% Convert power spectrum to dB scale
ospxdb = 10*log10(ospx);
% Plot onesided result first
figure(1)
plot(osf, ospxdb)
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
title('Power Spectral Density of Oversampled Binary Sequence with Different Oversampling Interval'),
hold on,

% Oversampling setup
M = 10.1; % Oversampling factor
delChip = 1/M; % Oversampling interval
delOffset = 0; % No offset
Np = length(code); % Length of binary sequence
Ns = round(M * Np); % Number of samples after oversampling
% Oversample the binary sequence
codeOS = oversampleSpreadingCode(code, delChip, delOffset, Ns, Np);
% Power spectrum estimation with pwelch
[ospx, osf] = pwelch(codeOS, hann(2^10), [], 2^10, 1/delChip); % one-sided result
% Convert power spectrum to dB scale
ospxdb = 10*log10(ospx);
plot(osf, ospxdb);
hold on,

% Oversampling setup
M = 20.1; % Oversampling factor
delChip = 1/M; % Oversampling interval
delOffset = 0; % No offset
Np = length(code); % Length of binary sequence
Ns = round(M * Np); % Number of samples after oversampling
% Oversample the binary sequence
codeOS = oversampleSpreadingCode(code, delChip, delOffset, Ns, Np);
% Power spectrum estimation with pwelch
[ospx, osf] = pwelch(codeOS, hann(2^10), [], 2^10, 1/delChip); % one-sided result
% Convert power spectrum to dB scale
ospxdb = 10*log10(ospx);
plot(osf, ospxdb);

legend('5.1', '10.1', '20.1');







clear all
close all
% Change the nfft interval
% Generate binary sequence with ±1 values
code = sign(randn([1, 2^14]));


% Oversampling setup
M = 10.1; % Oversampling factor
delChip = 1/M; % Oversampling interval
delOffset = 0; % No offset
Np = length(code); % Length of binary sequence
Ns = round(M * Np); % Number of samples after oversampling
% Oversample the binary sequence
codeOS = oversampleSpreadingCode(code, delChip, delOffset, Ns, Np);
% Power spectrum estimation with pwelch
[ospx, osf] = pwelch(codeOS, hann(2^10), [], 2^5, 1/delChip); % one-sided result
% Convert power spectrum to dB scale
ospxdb = 10*log10(ospx);
% Plot onesided result first
figure(1)
plot(osf, ospxdb)
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
title('Power Spectral Density of Oversampled Binary Sequence with Different nfft'),
hold on,

[ospx, osf] = pwelch(codeOS, hann(2^10), [], 2^8, 1/delChip); % one-sided result
% Convert power spectrum to dB scale
ospxdb = 10*log10(ospx);
plot(osf, ospxdb);
hold on,

[ospx, osf] = pwelch(codeOS, hann(2^10), [], 2^20, 1/delChip); % one-sided result
% Convert power spectrum to dB scale
ospxdb = 10*log10(ospx);
plot(osf, ospxdb);
legend('2^5','2^8','2^{20}')

%% 8 Lfrs Sequence
clear all;
clc;
n = 4;
ciVec = [2,3]';
a0Vec = [0,0,1]';
[lfsrSeq] = generateLfsrSequence(n,ciVec,a0Vec);


%% 9 Correlation Experiment
% a)
clear;clc;
%---- Multiple Runs
numberofRuns=10000;
Rseq12k0 = zeros(numberofRuns,1);
for i = 1:numberofRuns

    %----- Setup
    nStages = 10;                  % Number of stages in LFSR
    Tc = 1e-3/1023;               % Chip interval in seconds
    delChip = 3/217;              % Sampling interval in chips
    delOffset  = 0;               % Offset of first sample
    delt = delChip*Tc;            % Sampling interval in seconds
    fs = 1/delt;                  % Sampling frequency in Hz
    Np = 2^nStages - 1;           % Period of the sequence in chips
    Nr = 20;                      % Number of repetitions of the sequence
    Ns = round(Nr*Np/delChip);    % Number of samples of the sequence
    % codeType:
    % rand ---- Sequence derived from Matlab randn function
    % pi ------ Sequence derived from the digits of pi
    % mseq ---- Maximal-length sequence with n = nStages
    codeType = 'rand';

    %----- Generate codes
    X1 = zeros(Np,1);
    X2 = zeros(Np,1);
    if(strcmp(codeType,'rand'))
        X1 = sign(sign(randn(Np,1)) + 0.1);
        X2 = sign(sign(randn(Np,1)) + 0.1);
    elseif(strcmp(codeType,'pi'))
        [sPi,vPi] = pi2str(2*Np);
        X1 = vPi(1:Np) >= 5;
        X1 = 2*X1 - 1;
        X2 = vPi(Np+1:2*Np) >= 5;
        X2 = 2*X2 - 1;
    elseif(strcmp(codeType,'mseq'))
        ciVec1 = [9, 4]';
        ciVec2 = [9, 2]';
        a0Vec1 = [1;zeros(nStages-1,1)];
        a0Vec2 = ones(nStages,1);
        X1 = generateLfsrSequence(nStages,ciVec1,a0Vec1);
        X2 = generateLfsrSequence(nStages,ciVec2,a0Vec2);
        X1 = 2*X1 - 1;
        X2 = 2*X2 - 1;
    else
        error('Unrecognized code type');
    end

    %----- Compute the sequence crosscorrelation
    [Rseq12,iiVecSeq] = ccorr(X1,X2);

    %----- Store R(k=0)
    Rseq12k0(i)=Rseq12(1);
end
format shortG
disp(var(Rseq12))



%---------------------------------------------------------------------------

% b)

clear;clc;
%----- Setup
nStages = 10;                  % Number of stages in LFSR
Tc = 1e-3/1023;               % Chip interval in seconds
delChip = 3/217;              % Sampling interval in chips
delOffset  = 0;               % Offset of first sample
delt = delChip*Tc;            % Sampling interval in seconds
fs = 1/delt;                  % Sampling frequency in Hz
Np = 2^nStages - 1;           % Period of the sequence in chips
Nr = 20;                      % Number of repetitions of the sequence
Ns = round(Nr*Np/delChip);    % Number of samples of the sequence 
% codeType:
% rand ---- Sequence derived from Matlab randn function
% pi ------ Sequence derived from the digits of pi
% mseq ---- Maximal-length sequence with n = nStages
codeType = 'mseq';

%----- Generate codes
X1 = zeros(Np,1);
X2 = zeros(Np,1);
if(strcmp(codeType,'rand'))
  X1 = sign(sign(randn(Np,1)) + 0.1);
  X2 = sign(sign(randn(Np,1)) + 0.1);
elseif(strcmp(codeType,'pi'))
  [sPi,vPi] = pi2str(2*Np);
  X1 = vPi(1:Np) >= 5;
  X1 = 2*X1 - 1;
  X2 = vPi(Np+1:2*Np) >= 5;
  X2 = 2*X2 - 1;
elseif(strcmp(codeType,'mseq'))
ciVec1 = [10, 7]';  
ciVec2 = [10, 9, 7, 6]';
ciVec3 = [10, 9, 8, 7, 6, 5, 4, 1]';
ciVec4 = [10, 9, 8, 7, 6, 5, 4, 1]';
ciVec5 = [10, 9, 8, 6, 5, 1]';
ciVec6 = [10, 9, 8, 6, 4, 3]';
a0Vec1 = ones(nStages,1);
X1 = generateLfsrSequence(nStages,ciVec1,a0Vec1);
X2 = generateLfsrSequence(nStages,ciVec2,a0Vec1);
X3 = generateLfsrSequence(nStages,ciVec3,a0Vec1);
X4 = generateLfsrSequence(nStages,ciVec4,a0Vec1);
X5 = generateLfsrSequence(nStages,ciVec5,a0Vec1);
X6 = generateLfsrSequence(nStages,ciVec6,a0Vec1);
X1 = 2*X1 - 1;
X2 = 2*X2 - 1;
X3 = 2*X3 - 1;
X4 = 2*X4 - 1;
X5 = 2*X5 - 1;
X6 = 2*X6 - 1;
else
  error('Unrecognized code type');
end

%----- Compute the sequence autocorrelation
[Rseq1,iiVecSeq] = ccorr(X1,X1);
[Rseq2,iiVecSeq] = ccorr(X2,X2);
[Rseq3,iiVecSeq] = ccorr(X3,X3);
[Rseq4,iiVecSeq] = ccorr(X4,X4);
[Rseq5,iiVecSeq] = ccorr(X5,X5);
[Rseq6,iiVecSeq] = ccorr(X6,X6);

figure(1);clf;
subplot(611)
plot(iiVecSeq,Rseq1);
grid on;
ylabel('R_{X1}');
title('autocorrelation')
subplot(612)
plot(iiVecSeq,Rseq2);
grid on;
ylabel('R_{X2}');
title('autocorrelation')
subplot(613)
plot(iiVecSeq,Rseq3);
grid on;
ylabel('R_{X3}');
title('autocorrelation')
subplot(614)
plot(iiVecSeq,Rseq4);
grid on;
ylabel('R_{X4}');
title('autocorrelation')
subplot(615)
plot(iiVecSeq,Rseq5);
grid on;
ylabel('R_{X5}');
title('autocorrelation')
subplot(616)
plot(iiVecSeq,Rseq6);
grid on;
ylabel('R_{X6}');
title('autocorrelation')

%---------------------------------------------------------------------------

[Rseq12,iiVecSeq] = ccorr(X1,X2);
[Rseq13,iiVecSeq] = ccorr(X1,X3);
[Rseq14,iiVecSeq] = ccorr(X1,X4);
[Rseq15,iiVecSeq] = ccorr(X1,X5);
[Rseq16,iiVecSeq] = ccorr(X1,X6);
[Rseq23,iiVecSeq] = ccorr(X2,X3);
[Rseq24,iiVecSeq] = ccorr(X2,X4);
[Rseq25,iiVecSeq] = ccorr(X2,X5);
[Rseq26,iiVecSeq] = ccorr(X2,X6);
[Rseq34,iiVecSeq] = ccorr(X3,X4);
[Rseq35,iiVecSeq] = ccorr(X3,X5);
[Rseq36,iiVecSeq] = ccorr(X3,X6);
[Rseq45,iiVecSeq] = ccorr(X4,X5);
[Rseq46,iiVecSeq] = ccorr(X4,X6);
[Rseq56,iiVecSeq] = ccorr(X5,X6);

minBoundCrossCorr= sqrt(Np);
figure(2);clf;
maxValues=[minBoundCrossCorr; max(Rseq12); max(Rseq13); max(Rseq14) ;max(Rseq15) ;max(Rseq16);
    max(Rseq23); max(Rseq24); max(Rseq25); max(Rseq26); max(Rseq34); max(Rseq35);
    max(Rseq36); max(Rseq45) ;max(Rseq46) ;max(Rseq56)];

b= bar(maxValues);
b.FaceColor = 'flat';
b.CData(1,:) = [.5 0 .5];
title('Maximum Crosscorrealtion')
xlabel('Lag (samples)');
grid on;


%-------------------------------------------------


% %c)
% 
% clear;clc;
% %----- Setup
% nStages = 10;                  % Number of stages in LFSR
% Tc = 1e-3/1023;               % Chip interval in seconds
% delChip = 3/217;              % Sampling interval in chips
% delOffset  = 0;               % Offset of first sample
% delt = delChip*Tc;            % Sampling interval in seconds
% fs = 1/delt;                  % Sampling frequency in Hz
% Np = 2^nStages - 1;           % Period of the sequence in chips
% Nr = 20;                      % Number of repetitions of the sequence
% Ns = round(Nr*Np/delChip);    % Number of samples of the sequence 
% % codeType:
% % rand ---- Sequence derived from Matlab randn function
% % pi ------ Sequence derived from the digits of pi
% % mseq ---- Maximal-length sequence with n = nStages
% codeType = 'pi';
% 
% %----- Generate codes
% X1 = zeros(Np,1);
% X2 = zeros(Np,1);
% if(strcmp(codeType,'rand'))
%   X1 = sign(sign(randn(Np,1)) + 0.1);
%   X2 = sign(sign(randn(Np,1)) + 0.1);
% elseif(strcmp(codeType,'pi'))
%   [sPi,vPi] = pi2str(2*Np);
%   X1 = vPi(1:Np) >= 5;
%   X1 = 2*X1 - 1;
%   X2 = vPi(Np+1:2*Np) >= 5;
%   X2 = 2*X2 - 1;
% elseif(strcmp(codeType,'mseq'))
% ciVec1 = [9, 4]';  
% ciVec2 = [9, 2]';
% a0Vec1 = [1;zeros(nStages-1,1)];
% a0Vec2 = ones(nStages,1);
% X1 = generateLfsrSequence(nStages,ciVec1,a0Vec1);
% X2 = generateLfsrSequence(nStages,ciVec2,a0Vec2);
% X1 = 2*X1 - 1;
% X2 = 2*X2 - 1;
% else
%   error('Unrecognized code type');
% end
% 
% %----- Oversample code
% X1os = oversampleSpreadingCode(X1,delChip,delOffset,Ns,Np);
% X2os = oversampleSpreadingCode(X2,delChip,delOffset,Ns,Np);
% 
% %----- Compute autocorrelation 
% [R1,iiVec] = ccorr(X1os,X1os);
% [R2,iiVec] = ccorr(X2os,X2os);
% 
% %----- Compute crosscorrelation 
% [R12,iiVec] = ccorr(X1os,X2os);
% 
% figure(1);clf;
% subplot(211)
% plot(iiVec,R1/Ns);
% grid on;
% ylabel('R_{X1}');
% title('X1 and X2 autocorrelation')
% subplot(212)
% plot(iiVec,R2/Ns);
% grid on;
% ylabel('R_{X2}');
% xlabel('Lag (samples)');
% figure(2);clf;
% plot(iiVec,R12/Ns);
% title('X1 and X2 crosscorrelation')
% xlabel('Lag (samples)');
% grid on;
% ylabel('R_{X1,X2}');
% 
% %----- Compute power spectra
% % The scaling here ensures that sum(Si*delf) = 1 W, as expected, for i = 1, 2.
% S1 = abs(delt*fft(X1os)).^2/(Ns*delt);
% S2 = abs(delt*fft(X2os)).^2/(Ns*delt);
% S12 = abs(delt*fft(R12)).^2/(Ns*delt);
% delf = 1/(delt*Ns);
% fVec = [0:Ns-1]'*(delf);
% %----- Scale for 1-Hz delf for plotting
% % We'll present the power spectra in units of dBW/Hz, so we need to scale
% % the spectra so that they match a delf = 1 Hz frequency spacing.
% S1 = S1*delf;
% S2 = S2*delf;
% S12 = S12*delf;
% 
% figure(3);clf;
% subplot(211)
% plot(fVec/1e3,10*log10(S1));
% grid on;
% xlim([0 30]);
% ylim([-100,0]);
% title('X1 and X2 power spectral densities')
% ylabel('S_{X1}(f) (dBW/Hz)');
% subplot(212)
% plot(fVec/1e3,10*log10(S2));
% grid on;
% xlim([0 30]);
% ylim([-100,0]);
% ylabel('S_{X2}(f) (dBW/Hz)');
% xlabel('Frequency (kHz)');