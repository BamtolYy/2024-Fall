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
load(['C:\Users\gsh04\Desktop\2024-Fall\GPS\ps2\ps2MatlabFiles\' ...
    'prn31_22apr03_01hrs40min00sec_gmt_fl1_46_08mhz_250msec.mat']);

% Define the number of samples for the window
num_samples = 400;

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
figure(1)
Y_window = Y(520001:520001+400-1);
plot(real(Y_window),imag(Y_window),'o');
grid on;
figure(2)
subplot(2,1,1)
title('Real Component');
plot(real(Y_window));
grid on;
subplot(2,1,2)
title('img Component');
plot(imag(Y_window));
figure()
% Take a larger set of samples (e.g., 4000 samples)
Y_large = Y(1:4000);

% Plot the real and imaginary parts of the larger sample
figure(2);
plot(real(Y_large), imag(Y_large), 'o');
xlabel('Real Part');
ylabel('Imaginary Part');
title('4000 Complex Samples of the GPS Signal (Doughnut Shape)');
grid on;

