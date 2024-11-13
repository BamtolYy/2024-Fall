% topPerformAcqHypothesisCalcs
%
% Top-level script for performing acquisition calculations
%% change N
%----- Setup
clear; clc;
CN0 = [27,30,33,36,39,42,45,48,51];
Pd=1;
testN = flip([1:500],2);
for n = 1:length(CN0)
    Pd=1;
    for k = 1: length(testN)
        s.C_N0dBHz = CN0(n);
        s.N = testN(k);
        s.PfaAcq = 0.0001;
        s.Ta = 0.001;
        s.fMax = 7000;
        s.nCodeOffsets = 1023*5;
        s.ZMax = 1000;
        s.delZ = 0.1;

        %----- Execute
        [pZ_H0,pZ_H1,lambda0,Pd,ZVec] = performAcqHypothesisCalcs(s);
        if Pd <= 0.95
            N(n) =testN(k)+1;
            break;
        end
    end
end

% Plot Results
figure,
plot(CN0,[N,1,1])
xlabel('C/N0')
ylabel('Number of Samples(N)')
title('Number of Samples required to Achieve P_D = 0.95 or more')



%% Change Ta
%----- Setup
clear; clc;
CN0 = linspace(5,50,10);
testTa = flip([0.0001:0.0001:0.1]);
for n = 1:length(CN0)
    Pd=1;
    for k = 1: length(testTa)
        s.C_N0dBHz =CN0(n);
        s.N = 1;
        s.PfaAcq = 0.001;
        s.Ta = testTa(k);
        s.fMax = 7000;
        s.nCodeOffsets = 1023*5;
        s.ZMax = 1000;
        s.delZ = 0.1;

        %----- Execute
        [pZ_H0,pZ_H1,lambda0,Pd,ZVec] = performAcqHypothesisCalcs(s);
        if Pd <= 0.95
            T(n) =testTa(k);
            break;
        end
    end
end

% Plot Results
figure,
plot(CN0,T)
xlabel('C/N0')
ylabel('Accumulation Time (s) ')
title('Accumulation Time required to achieve P_D = 0.95 or more')


fprintf(['Coherent integration is more efficient than noncoherent integration. \n' ...
    'For C/N0 = 30, you require more 99 accumulations for the noncoherent \n' ...
    'integration with accumulation time = 0.001 second, which means the total \n' ...
    'accumulation time is 0.099 seconds. In contrast, coherernt integration \n' ...
    'required only one accumulation for 0.0331 second to achieve PD = 0.95.'])

% 
% %----- Visualize the results
% figure(2);
% [pmax,iimax] = max(pZ_H1);
% Zmax = ZVec(iimax);
% clf;
% ash = area(ZVec,pZ_H0);
% set(get(ash,'children'), 'facecolor', 'g', 'linewidth', 2, 'facealpha', 0.5);
% hold on;
% ash = area(ZVec,pZ_H1);
% set(get(ash,'children'), 'facecolor', 'b', 'linewidth', 2, 'facealpha', 0.5);
% linemax = 1/5*max([pZ_H0(:);pZ_H1(:)]);
% line([lambda0,lambda0],[0,linemax], 'linewidth', 2, 'color', 'r');
% xlim([0 max(Zmax*2,lambda0*1.5)]);
% ylabel('Probability density');
% xlabel('Z');
% fs = 12;
% title('GNSS Acquisition Hypothesis Testing Problem');
% disp(['Probability of acquisition false alarm (PfaAcq): ' ...
%     num2str(s.PfaAcq)]);
% disp(['Probability of detection (Pd): ' num2str(Pd)]);
% text(lambda0,linemax*1.1, ['\lambda_0 = ' num2str(lambda0) ], ...
%     'fontsize',fs);
% shg