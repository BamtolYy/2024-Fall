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
% References:
%
%
%+==============================================================================+

% Compute the intensity
I   = abs(zkhist).^2;
I_mean  = mean(I);
S4 = std(I)/I_mean;

zkhist_centered = zkhist - mean(zkhist);
[R, lags] = xcorr(zkhist_centered, 'normalized');
figure
plot(lags,R)
R = R(lags >= 0);
lags = lags(lags >= 0);
e_frac = exp(-1);
tau0_index = find(R <= e_frac, 1);  % First index where ACF falls to 1/e
tau0 = tkhist(tau0_index); 
xlim([-100 100]);


