%% Exam 1 Problem 7
clear all;
close all;
clc

nStages = 9;        % Number of stages in LFSR
Np = 2^nStages - 1;  % Maximum length sequence

ciVec1 = [9,4]';
ciVec2 = [9,6,5,3]';
ciVec3 = [9,8,5,3]';
ciVec4 = [9,6,4,3]';
ciVec5 = [9,8,7,3]';
ciVec6 = [9,2]';
a0Vec1 = [1;zeros(nStages-1,1)];
a0Vec2 = [0;1;1;zeros(nStages-3,1)];
a0Vec4 = [1;0;ones(nStages-2,1)];

X1 = generateLfsrSequence(nStages,ciVec1,a0Vec1);
X2 = generateLfsrSequence(nStages,ciVec2,a0Vec1);
X3 = generateLfsrSequence(nStages,ciVec3,a0Vec1);
X4 = generateLfsrSequence(nStages,ciVec4,a0Vec1);
X5 = generateLfsrSequence(nStages,ciVec5,a0Vec1);
X6 = generateLfsrSequence(nStages,ciVec6,a0Vec1);

X1 = 2*X1 - 1;
X2 = 2*X2 - 1;
X3 = 2*X3 - 1;
X4 = 2*X4 - 1;
X5 = 2*X5 - 1;
X6 = 2*X6 - 1;

% G12 = generateGoldLfsrSequence(nStages,ciVec1,ciVec2,a0Vec1,a0Vec2);
% G14 = generateGoldLfsrSequence(nStages,ciVec1,ciVec4,a0Vec1,a0Vec4);
% G24 = generateGoldLfsrSequence(nStages,ciVec2,ciVec4,a0Vec2,a0Vec4);
% G12 = 2*G12 - 1;
% G14 = 2*G14 - 1;
% G24 = 2*G24 - 1;
% [Gseq1214,iiVecSeq] = ccorr(G12,G14);
% [Gseq1224,iiVecSeq] = ccorr(G12,G24);
% [Gseq2414,iiVecSeq] = ccorr(G24,G14);
% figure;clf;
% plot(iiVecSeq,Gseq1214);
% title(['Potential G12 and G14 crosscorrelation'])
% xlabel('Lag (samples)');
% grid on;
% figure;clf;
% plot(iiVecSeq,Gseq1224);
% title(['Potential G24 and G14 crosscorrelation'])
% xlabel('Lag (samples)');
% grid on;
% figure;clf;
% plot(iiVecSeq,Gseq2414);
% title(['Potential G2 and G3 crosscorrelation'])
% xlabel('Lag (samples)');
% grid on;

% Calculate t(n)
if mod(nStages, 2) == 1  % n is odd
    t_n = 1 + 2^((nStages + 1) / 2);
else  % n is even
    t_n = 1 + 2^((nStages + 2) / 2);
end

% Define preferred correlation values
preferredCorrelations = [-t_n, -1, t_n - 2]



%%

[R12,iiVecSeq] = ccorr(X1,X2);
[R14,iiVecSeq] = ccorr(X1,X4);
[R24,iiVecSeq] = ccorr(X2,X4);
figure;clf;
plot(iiVecSeq,R12);
title(['Potential X1 and X2 crosscorrelation'])
ylabel('R_{X1,X2}');
xlabel('Lag (samples)');
grid on;
figure;clf;
plot(iiVecSeq,R14);
title(['Potential X1 and X4 crosscorrelation'])
ylabel('R_{X1,X4}');
xlabel('Lag (samples)');
grid on;
figure;clf;
plot(iiVecSeq,R24);
title(['Potential X2 and X4 crosscorrelation'])
ylabel('R_{X2,X4}');
xlabel('Lag (samples)');
grid on;


