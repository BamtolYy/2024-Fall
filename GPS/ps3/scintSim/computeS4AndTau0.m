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
I2_mean = mean(I.^2);
S4 = sqrt((I2_mean-I_mean^2)/I_mean^2)



tau0 = sqrt(1-S4^2) /(1-sqrt(1-S4^2))

