clear all
clc
%%% Problem set 2

%% 1a

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

% Generate binary sequence with ±1 values
code = randi([-1, 1], [1, 2^14]);
code(~code) = -1; % Ensure values are only ±1

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
code = randi([-1, 1], [1, 2^9]);
code(~code) = -1; % Ensure values are only ±1
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
code = randi([-1, 1], [1, 2^14]);
code(~code) = -1; % Ensure values are only ±1
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
code = randi([-1, 1], [1, 2^20]);
code(~code) = -1; % Ensure values are only ±1
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
code = randi([-1, 1], [1, 2^14]);
code(~code) = -1; % Ensure values are only ±1
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
% Change the oversampling interval
% Generate binary sequence with ±1 values
code = randi([-1, 1], [1, 2^14]);
code(~code) = -1; % Ensure values are only ±1
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

