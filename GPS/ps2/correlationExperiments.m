% correlationExperiments.m
%
% Experiment with properties of pseudorandom sequences.


clear;clc;
%----- Setup
nStages = 10;                  % Number of stages in LFSR
Tc = 1e-3/1023;               % Chip interval in seconds
delChip = 3/217;              % Sampling interval in chips
delOffset  = 0;               % Offset of first sample
delt = delChip*Tc;            % Sampling interval in seconds
fs = 1/delt;                  % Sampling frequency in Hz
Np = 2^nStages - 1;           % Period of the sequence in chips
Nr = 20;                      % Number of repetitions of the sequence
Ns = round(Nr*Np/delChip);    % Number of samples of the sequence
% codeType:
% rand ---- Sequence derived from Matlab randn function
% pi ------ Sequence derived from the digits of pi
% mseq ---- Maximal-length sequence with n = nStages
z = {'rand','pi','mseq'};

for i = 1:3
    codeType = z(i);
%----- Generate codes
X1 = zeros(Np,1);
X2 = zeros(Np,1);
if(strcmp(codeType,'rand'))
    X1 = sign(sign(randn(Np,1)) + 0.1);
    X2 = sign(sign(randn(Np,1)) + 0.1);
elseif(strcmp(codeType,'pi'))
    [sPi,vPi] = pi2str(2*Np);
    X1 = vPi(1:Np) >= 5;
    X1 = 2*X1 - 1;
    X2 = vPi(Np+1:2*Np) >= 5;
    X2 = 2*X2 - 1;
elseif(strcmp(codeType,'mseq'))
    ciVec1 = [9, 4]';
    ciVec2 = [9, 2]';
    a0Vec1 = [1;zeros(nStages-1,1)];
    a0Vec2 = ones(nStages,1);
    X1 = generateLfsrSequence(nStages,ciVec1,a0Vec1);
    X2 = generateLfsrSequence(nStages,ciVec2,a0Vec2);
    X1 = 2*X1 - 1;
    X2 = 2*X2 - 1;
else
    error('Unrecognized code type');
end

%----- Compute the sequence autocorrelation
[Rseq1,iiVecSeq] = ccorr(X1,X1);
[Rseq2,iiVecSeq] = ccorr(X2,X2);

%----- Compute the sequence crosscorrelation
[Rseq12,iiVecSeq] = ccorr(X1,X2);

%----- Oversample code
X1os = oversampleSpreadingCode(X1,delChip,delOffset,Ns,Np);
X2os = oversampleSpreadingCode(X2,delChip,delOffset,Ns,Np);

%----- Compute autocorrelation
[R1,iiVec] = ccorr(X1os,X1os);
[R2,iiVec] = ccorr(X2os,X2os);

%----- Compute crosscorrelation
[R12,iiVec] = ccorr(X1os,X2os);

%----- Compute power spectra
% The scaling here ensures that sum(Si*delf) = 1 W, as expected, for i = 1, 2.
S1 = abs(delt*fft(X1os)).^2/(Ns*delt);
S2 = abs(delt*fft(X2os)).^2/(Ns*delt);
S12 = abs(delt*fft(R12)).^2/(Ns*delt);
delf = 1/(delt*Ns);
fVec = [0:Ns-1]'*(delf);

%----- Scale for 1-Hz delf for plotting
% We'll present the power spectra in units of dBW/Hz, so we need to scale
% the spectra so that they match a delf = 1 Hz frequency spacing.
S1 = S1*delf;
S2 = S2*delf;
S12 = S12*delf;


period = round(Np/delChip);
%----- Plot
% figure;clf;
% subplot(211)
% plot(iiVec,R1/Ns);
% grid on;
% ylabel('R_{X1}');
% title('X1 and X2 autocorrelation')
% xlim([-period/2,period/2])
% subplot(212)
% plot(iiVec,R2/Ns);
% grid on;
% ylabel('R_{X2}');
% xlabel('Lag (samples)');
% xlim([-period/2,period/2])
% hold off

% figure;clf;
% plot(iiVec,R12/Ns);
% title('X1 and X2 crosscorrelation')
% xlabel('Lag (samples)');
% grid on;
% ylabel('R_{X1,X2}');
% xlim([-period/2,period/2])
% hold off

figure;clf;
subplot(211)
plot(fVec/1e3,10*log10(S1));
grid on;
xlim([0 30]);
ylim([-100,0]);
title('X1 and X2 power spectral densities')
ylabel('S_{X1}(f) (dBW/Hz)');
subplot(212)
plot(fVec/1e3,10*log10(S2));
grid on;
xlim([0 30]);
ylim([-100,0]);
ylabel('S_{X2}(f) (dBW/Hz)');
xlabel('Frequency (kHz)');
hold off


% What is the ratio of max(Rx1) to max(Rx1x2)?
% if(strcmp(codeType,'rand'))
%     R1_period_rand = R1(1:period);
%     R12_period_rand = R12(1:period);
%     Ratio_rand= max(R1_period_rand)/max(R12_period_rand);
%     format longG
%     disp(max(R1_period_rand))
%     disp(max(R12_period_rand))
% elseif(strcmp(codeType,'pi'))
%     R1_period_pi = R1(1:period);
%     R12_period_pi = R12(1:period);
%     Ratio_pi= max(R1_period_pi)/max(R12_period_pi);
%     disp(max(R1_period_pi))
%     disp(max(R12_period_pi))
% elseif(strcmp(codeType,'mseq'))
%     R1_period_mseq = R1(1:period);
%     R12_period_mseq = R12(1:period);
%     Ratio_mseq= max(R1_period_mseq)/max(R12_period_mseq);
%     disp(max(R1_period_mseq))
%     disp(max(R12_period_mseq))
% end

if(strcmp(codeType,'rand'))
    Ratio_rand= max(R1)/max(R12);
    format longG
    disp(max(R1))
    disp(max(R12))
elseif(strcmp(codeType,'pi'))
    Ratio_pi= max(R1)/max(R12);
    disp(max(R1))
    disp(max(R12))
elseif(strcmp(codeType,'mseq'))
    Ratio_mseq= max(R1)/max(R12);
    disp(max(R1))
    disp(max(R12))
end

end

varof30power1=10*log10(var(S1(1:30)))
varof30power2=10*log10(var(S2(1:30)))

Nsep=fVec/delf
