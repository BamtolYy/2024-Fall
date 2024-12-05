close all; clear all; clc;

%% a) Parameters
Ta = 10 / 1000;      % seconds
Bn = 10;             % Hz

%% b) Define Loop Filter
loopOrder = 3;
[s.Ad, s.Bd, s.Cd, s.Dd, Bn_act] = configureLoopFilter(Bn, Ta, loopOrder);

%% Generate Fictitious Phase Time History
fs = 10000;                 % Hz; Fictitious sampling rate of signal
T = 1 / fs;                 % seconds; Sampling interval
t = (0:T:1)';               % Time vector
fd = 1;                     % Frequency of phase history
PhaseHist = 2 * pi * fd * t; % Raw Phase History (linearly increasing phase)

%% Average over sub-interval Ta
Nk = floor(Ta / T);         % Number of samples in an accumulation
N_sub = floor(length(t) / Nk); % Number of sub-intervals in the signal
xk0 = zeros(loopOrder - 1, 1);
s.xk = xk0;                 % Initial state of the loop filter
thetahat = zeros(length(t), 1); % Estimated phase
deltaThetak = 0;            % Initial phase error
s.Ip = cos(deltaThetak);    % Initial in-phase value
s.Qp = sin(deltaThetak);    % Initial quadrature value

%% Main Loop for Phase Tracking
for ii = 1:N_sub
    % Define start and end indices for the current accumulation interval
    start = (ii - 1) * Nk + 1;
    endpt = ii * Nk;

    % Update PLL state using the previous error signal
    [xkp1, vk] = updatePll(s);
    s.xk = xkp1;

    % Update thetahat over the current interval
    for jj = start:endpt
        if jj == start
            thetahat(jj) = thetahat(jj) + vk * T; % First point based on previous thetahat
        else
            thetahat(jj) = thetahat(jj - 1) + vk * T;
        end
    end

    % Compute the phase difference (error) over the current interval
    deltaTheta = PhaseHist(start:endpt) - thetahat(start:endpt);
    deltaThetak = deltaTheta(end); % Take the last value of deltaTheta for the next update

    % Update `Ip` and `Qp` using the new phase error
    s.Ip = cos(deltaThetak);
    s.Qp = sin(deltaThetak);
end

%% Plot Results
plot(t, PhaseHist, t, thetahat, '--')
legend('True', 'Estimate')
xlabel('Time (s)')
ylabel('Phase (rad)')
title('True vs Estimated Phase')
