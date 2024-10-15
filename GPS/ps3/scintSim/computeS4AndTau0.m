function [S4,tau0] = computeS4AndTau0(zkhist,tkhist)
% computeS4AndTau0 : Compute the scintillation index S4 and the decorrelation
% time tau0 corresponding to the input complex channel
% response function time history zkhist.
%
%
% INPUTS
%
% zkhist ----- Nt-by-1 vector containing the normalized complex scintillation
% time history in the form of averages over Ts with sampling
% interval Ts. zkhist(kp1) is the average over tk to tkp1.
%
% tkhist ----- Nt-by-1 vector of time points corresponding to zkhist.
%
%
% OUTPUTS
%
% S4 --------- Intensity scintillation index of the scintillation time history
% in zkhist, equal to the mean-normalized standard deviation of
% the intensity abs(zkhist).^2.
%
% tau0 ------- The decorrelation time of the scintillation time history in
% zkhist, in seconds.
%
%
%+------------------------------------------------------------------------------+
% References:v
%
%
%+==============================================================================+

% Compute the intensity
I   = abs(zkhist).^2;
I_mean  = mean(I);
S4 = std(I)/I_mean;


% Subtract the mean intensity to focus on fluctuations
I_fluct = I - mean(I);

% Calculate the autocorrelation of the intensity fluctuations
[acf, lags] = xcorr(I_fluct, 'coeff');  % Autocorrelation of intensity fluctuations

% Extract positive lags (to avoid negative lag calculations)
acf = acf(lags >= 0);
lags = lags(lags >= 0);

% Find the time where the autocorrelation drops to 1/e
e_frac = exp(-1);
tau0_index = find(acf <= e_frac, 1);  % First index where ACF falls to 1/e
tau0 = tkhist(tau0_index);  % Time at which ACF is 1/e

