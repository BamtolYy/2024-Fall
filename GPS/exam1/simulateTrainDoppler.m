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



% [t,rTrain] = ode45(@(t,x) trainODE(t,x,vTrain),tspan,x0 );


tspan = [t0:delt:t0+N*delt]';
rObs=[xObs,dObs];
fDVec=zeros(length(tspan)-1,1);
rTrain=zeros(length(tspan),2);
TOF = zeros(length(tspan),1);
vlos = zeros(length(tspan),1);

j=1;

for j = 1:length(tspan)
    %Find TOF through iteration
    error = 9999;
    i=1;
    TOF_guess(1)=0;
    rTrain(j,1) = x0 + vTrain*tspan(j);
    while error > 0.0000001
        r_tTOF=TOF_guess(i)*vTrain;
        % From iterative method: norm[r_obs(t)-r_train(t-TOF)] gives you
        % the distance from the observer to the Train.
        % Divide this by Speed of sound or the speed that the horn travels at.
        % The result gives you the first guess of the TOF or the time of flight
        % of the horn travel from the train to the observer.

        TOF_guess(i+1)=norm(rObs-(rTrain(j,:)-[r_tTOF 0]))/vs;
        error = abs((TOF_guess(i+1)-TOF_guess(i))/TOF_guess(i+1)*100);
        i=i+1;
    end
    TOF(j) = TOF_guess(end);
end

for jj = 1:length(tspan)
    alpha = atan2(rTrain(jj,2)-dObs, (rTrain(jj,1)-TOF(jj)*vTrain)-xObs);
    vlos(jj) = cos(alpha)*vTrain;
    fr = fc/(1+vlos(jj)/vs);
    fDVec(jj) = fr- fc;
end






tVec=tspan;
