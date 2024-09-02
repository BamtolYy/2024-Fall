function [fDVec,tVec] = ...
    simulateTrainDoppler(fc, vTrain, t0, x0, xObs, dObs, delt, N, vs)
% simulateTrainDoppler : Simulate the train horn Doppler shift scenario.
%
% INPUTS
%
% fc------ train horn frequency, in Hz
%
% vTrain-- constant along-track train speed, in m/s seconds
%
% t0------ time at which train passed the along-track coordinate x0, in
%          seconds
%
% x0------ scalar along-track coordinate of train at time t0, in meters
%
% xObs---- scalar along-track coordinate of observer, in meters
%
% dObs---- scalar cross-track coordinate of observer, in meters (i.e., 
%          shortest distance of observer from tracks)
%
% delt---- measurement interval, in seconds
%
% N------- number of measurements
%
% vs------ speed of sound, in m/s
%
%
% OUTPUTS
%
% fDVec--- N-by-1 vector of apparent Doppler frequency shift measurements as
%          sensed by observer at the time points in tVec
%
% tVec---- N-by-1 vector of time points starting at t0 and spaced by delt
%          corresponding to the measurements in fDVec
%
%+------------------------------------------------------------------------------+
% References: main.pdf (Dr. Todd Humphreys)
%
%
% Author: Bonsuck Koo
%+==============================================================================+



vlos = dt(dObs-xobs)vs;
fr = fc/(1+vlos/vs);
fDVec = fr- fc;
