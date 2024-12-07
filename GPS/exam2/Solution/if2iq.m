function [IVec,QVec] = if2iq(xVec,T,fIF)
% IF2IQ : Convert intermediate frequency samples to baseband I and Q samples.
%
% Let x(n) = I(n*T)*cos(2*pi*fIF*n*T) − Q(n*T)*sin(2*pi*fIF*n*T) be a
% discrete−time bandpass signal centered at the user−specified intermediate
% frequency fIF, where T is the bandpass sampling interval.Then this
% function converts the bandpass samples to quadrature samples from a complex
% discrete−time baseband representation of the form xl(m*Tl) = I(m*Tl) +
% j*Q(m*Tl), where Tl = 2*T.
%
%
% INPUTS
%
% xVec −−−−−−−− N−by−1 vector of intermediate frequency samples with
%               sampling interval T.
%
% T −−−−−−−−−−− Sampling interval of intermediate frequency samples, in
%               seconds.
%
% fIF −−−−−−−−− Intermediate frequency of the bandpass signal, in Hz.
%
% OUTPUTS
%
% IVec −−−−−−−− N/2−by−1 vector of in−phase baseband samples.
%
% QVec −−−−−−−− N/2−by−1 vector of quadrature baseband samples.
%
%
%+−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−+
% References:
%
%
%+==============================================================================+

%−−−−− Setup
N = 2*floor(length(xVec)/2);
xVec = xVec(1:N);
tIfVec = [0:N-1]'*T;
%−−−−− Downconvert
cVec = sqrt(2)*cos(2*pi*fIF*tIfVec);
sVec = -sqrt(2)*sin(2*pi*fIF*tIfVec);
bbcVec = xVec.*cVec;
bbsVec = xVec.*sVec;
%−−−−− Filter and resample
IVec = decimate(bbcVec,2);
QVec = decimate(bbsVec,2);