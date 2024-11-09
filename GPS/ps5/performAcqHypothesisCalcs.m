function [pZ_H0,pZ_H1,lambda0,Pd,ZVec] = performAcqHypothesisCalcs(s)
% performAcqHypothesisCalcs : Calculate the null-hypothesis and alternative
%                             hypothesis probability density functions and the
%                             decision threshold corresponding to GNSS signal
%                             acquisition with the given inputs.
% Z is the acquisition statistic:
%
%
% Z =        N  |      |^2
%           sum |  Sk  |
%           k=1 |      |
%
%               _           _
%           N  |             |
%   =      sum | Ik^2 + Qk^2 |
%          k=1 |_           _|
%
%
%
%
%
%
% where Sk = rhok + nk
%          = Ik + j*Qk
% and nk = nIk + j*nQk
%
%
% with nIk ~ N(0,1), nQk ~ N(0,1), E[nIk nIi] = E[nQk nQi] = 1 for k = i and 0
% for k != i, and E[nIk nQi] = 0 for all k,i. The amplitude rhok is related
% to familiar parameters Nk, Abark, and sigma_IQ by rhok =
% (Nk*Abark)/(2*sigma_IQ), i.e., it is the magnitude of the usual complex
% baseband phasor normalized by sigma_IQ.
%
% Under H0, the statistic Z is distributed as a chi-square distribution with
% 2*N degrees of freedom; under H1, it is distributed as a noncentral
% chi-square distribution with lambda = N*rhok^2 and 2*N degrees of freedom.
%
% The total number of cells in the search grid is assumed to be nCells =
% nCodeOffsets*nFreqOffsets, where nFreqOffsets = 2*fMax*Ta and Ta = Na*T is
% the total coherent accumulation time. Here, Na is the average value of the
% number of samples in each accumulation, Nk.
%
% INPUTS
%
% s-------- A structure containing the following fields:
%
%       C_N0dBHz------- Carrier to noise ratio in dB-Hz.
%
%       Ta------------- Coherent accumulation interval, in seconds.
%
%       N-------------- The number of accumulations summed noncoherently to
%                       get Z.
%
%       fMax----------- Frequency search range delimiter. The total
%                       frequency search range is +/- fMax.
%
%       nCodeOffsets--- Number of statistically independent code offsets in
%                       the search range.
%
%       PfaAcq--------- The total acquisition false alarm probability.
%                       This is the probability that the statistic Z
%                       exceeds the threshold lambda in any one of the
%                       search cells under the hypothesis H0. One can
%                       derive the false alarm probability for *each*
%                       search cell from PfaAcq. This procedure is
%                       straightforward if we assume that the detection
%                       statistics from the search cells are independent
%                       of one another.
%       ZMax----------- The maximum value of Z that will be considered.
%
%       delZ----------- The discretization interval used for the
%                       independent variable Z. The full vector of Z
%                       values considered is thus ZVec = [0:delZ:ZMax].
%
%
% OUTPUTS
%
% pZ_H0---------- The probability density of Z under hypothesis H0.
%
% pZ_H1---------- The probability density of Z under hypothesis H1.
%
% lambda0-------- The detection threshold.
%
% Pd------------- The probability of detection.
%
% Zvec----------- The vector of Z values considered.
%
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+

CN0 = 10^(s.C_N0dBHz/10);    % Convert to linear scale

rho=sqrt(CN0*2*s.Ta*s.N);
lambda = s.N*rho^2;
Zvec = [0:s.delZ:s.ZMax];
pZ_H0=chi2pdf(Zvec,2*s.N);
pZ_H1=ncx2pdf(Zvec,2*s.N,lambda);
lambda0=1;
Pd= 1- ncx2cdf(lambda0, 2* s.N,lambda0);

