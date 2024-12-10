function [Ske,Skp,Skl] = correlate(Y,fIF,tstart,vTheta,thetaHatStart,teml,prn,T,Ta); 

% T = sampling Time
fs = 1/T;

Nk = round(Ta/T);

jk = floor(tstart/T);

tVec = [jk:jk+Nk-1]'*T;


[codeEarly, codePrompt, codeLate] = generateEarlyPromptLateCodes(prn, tstart, teml, fs, Nk);


ThetaVec = [2*pi*fIF*tVec+thetaHatStart+vTheta*[0:Nk-1]'*T];
carrierVec = exp(-i*ThetaVec);
lVeckEarly = carrierVec.*codeEarly;
lVeckPrompt = carrierVec.*codePrompt;
lVeckLate = carrierVec.*codeLate;

temlt = teml*1e-3/1023;
jkp = floor(tstart/T);
jke = floor((tstart-temlt/2)/T);
jkl = floor((tstart+temlt/2)/T);
xVecPrompt = Y(jkp:jkp+Nk-1);
xVecEarly = Y(jke:jke+Nk-1);
xVecLate = Y(jkl:jkl+Nk-1);
% xVecPrompt = Y(jkp:jkp+Nk-1);
% xVecEarly = Y(jke:jke+Nk-1);
% xVecLate = Y(jkl:jkl+Nk-1);

Skp = sum(xVecPrompt.*lVeckPrompt);
Ske = sum(xVecEarly.*lVeckEarly);
Skl = sum(xVecLate.*lVeckLate);
disp('prompt')
abs(Skp)^2
abs(Ske)^2
abs(Skl)^2