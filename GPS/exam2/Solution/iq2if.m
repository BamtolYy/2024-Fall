function [xVec] = iq2if(IVec,QVec,Tl,fIF)
% IQ2IF : Convert baseband I and Q samples to intermediate frequency samples.
%
% Let xl(m*Tl) = I(m*Tl) + j*Q(m*Tl) be a discrete−time baseband
% representation of a bandpass signal. This function converts xl(n) to a
% discrete−time bandpass signal x(n) = I(n*T)*cos(2*pi*fIF*n*T) −
% Q(n*T)*sin(2*pi*fIF*n*T) centered at the user−specified intermediate
% frequency fIF, where T = Tl/2.
%
%
% INPUTS
%
% IVec −−−−−−−− N−by−1 vector of in−phase baseband samples.
%
% QVec −−−−−−−− N−by−1 vector of quadrature baseband samples.
%
% Tl −−−−−−−−−− Sampling interval of baseband samples (complex sampling
%               interval), in seconds.
%
% fIF −−−−−−−−− Intermediate frequency to which the baseband samples will
%               be up−converted, in Hz.
%
%
% OUTPUTS
%
% xVec −−−−−−−− 2*N−by−1 vector of intermediate frequency samples with
%               sampling interval T = Tl/2.
%
%+−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−+
% References:
%
%
%+==============================================================================+
%−−−−− Setup
N = length(IVec);
if(length(QVec) ~= N)
    error('IVec and QVec must be of same length');
end
T = Tl/2;
tIfVec = [0:2*N-1]'*T;
%−−−−− Interpolate and resample the baseband quadrature signals
IVecInterp = interp(IVec,2);
QVecInterp = interp(QVec,2);
%−−−−− Upconvert
cVec = sqrt(2)*cos(2*pi*fIF*tIfVec);
sVec = sqrt(2)*sin(2*pi*fIF*tIfVec);
xVec = IVecInterp.*cVec - QVecInterp.*sVec;