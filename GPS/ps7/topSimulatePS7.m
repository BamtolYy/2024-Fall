clear; clc;
%% Coarse Search
% Coarse Search Parameter
fDRange = [-5000:100:5000];
prn = [1:37];
NC = 10;                    % Noncoherent sum number
Ta = 0.001;                 % Accumulation time in seconds
coarsets = zeros(1, length(prn));
coarsefD = zeros(1, length(prn));
% Estimate
for m = 1:length(prn)
    [ts, fD] = acquisition(prn(m),fDRange,NC,Ta);
    coarsets(m) = ts;
    coarsefD(m) = fD;
end


%% Fine Search
prnFine =  find(~isnan(coarsefD));
Ta =0.01;
NC = 8;
for h = 1:length(prnFine)
    fDmaxFine = coarsefD(prnFine(h))+50;
    fDminFine = coarsefD(prnFine(h))-50;
    fDRangeFine = [fDminFine:1:fDmaxFine];
    [ts, fD] = acquisition(prnFine(h),fDRangeFine,NC,Ta);
    tsFine(h) = ts;
    fDFine(h) = fD;
end

