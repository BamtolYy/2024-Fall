close all; clear all; clc;
%% a) Generate 5ms worth of signal and its replica

% ----Genererate Code
% LFSR Parameters:
nStages      = 10;
ciVec1       = [10, 3]';
ciVec2       = [10, 9, 8, 6, 3, 2,]';
a0Vec1       = ones(nStages,1);
a0Vec2       = ones(nStages,1);
% Oversampling Parameters:
fs = 40e6/7;                % Sampling Rate
Tc = 1e-3/1023;             % Chip interval in seconds
T  = 1/fs;                  % Sampling time interval in seconds
delChip = T/Tc;             % Sampling interval in chips
Np = 2^nStages - 1;         % Period of the sequence in chips
Ns = ceil(0.005/T);               % Number of Samples
codeOS = zeros(Ns,1);
G2tab = [2,6];
for j = 1
    [GoldSeq] = generateGoldLfsrSequenceCA(nStages,ciVec1,ciVec2,a0Vec1,...
        a0Vec2,G2tab(j,:));
    % Make code +1/-1 not +1/0
    GoldSeq = 2*GoldSeq - 1;
    % Oversample Code: It makes sense to oversample code, since the code
    % embedded within the signal is sampled at a higher rate than its chip
    % rate. Assuming that the code I generate is sampled at the chip rate,
    % oversampling my code I generated at the rate the signal is sampled
    % will allow my code to correlate with the code embedded in the signal
    GoldSeqOS = oversampleSpreadingCode(GoldSeq,delChip,0,Ns,Np);
    codeOS(:,j) = GoldSeqOS;
end

% ----Generate Noise
sigman = sqrt(1/2);
a      = sigman.*randn(Ns,1);
b      = sigman.*randn(Ns,1);
nj     = a +1i*b;

% ----x(j)
CN0   = 10^(45/10);
A     = sqrt(4*sigman^2*T*CN0);
% Delay code by 1e-3 seconds
ts    = 1e-3;                           % True delay in seconds
shift = floor(ts/T);
xj    = A*circshift(codeOS,shift) + nj;

% ----l(j)
delt_s        = 0.25e-6;
shiftReplica = floor((ts + delt_s)/T);            % Estimate Delay in seconds
lj           = circshift(codeOS,shiftReplica);

%% b) Compute Rtilde_xl(m)
m        =  -10:1:10;
Rxltilde = zeros(length(m),1);
jk = 1024;
for ii = 1:length(m)
    Nk       = ceil(Np*4/delChip);
    Rxltilde(ii) = 1/Nk*sum(xj(jk:jk+Nk-1).*lj(jk-m(ii):jk+Nk-1-m(ii)));
end

%% c) Nearest-Sample Estimator

[~,mIndex] = max(abs(Rxltilde));
delt_ns  = -m(mIndex)*T;
error_ts = abs(delt_ns - delt_s);

fprintf('Coarse %ct_s estimate error: %g %cs. \n',948-32,error_ts*1e6, 956 )

%% d) Least Squares Estimator
dt           = linspace(-T,T,1024);
normRxltilde = zeros(length(m),1);
Rc           = zeros(length(m),length(dt));
for kk = 1:length(m)
    for jj = 1:length(dt)
        Rc(kk,jj) = 1-abs(m(kk)*T-m(mIndex)*T-dt(jj))/Tc;
    end
    normRxltilde(kk)=abs(Rxltilde(kk))/abs(Rxltilde(mIndex));
end

MSE = zeros(length(dt),1);
for aa  = 1:length(dt)
    MSE(aa) = 1/21*sum(abs(normRxltilde-Rc(:,aa)).^2);
end

[~,minIndex] = min(MSE);
dt_r  = dt(minIndex);
dt_ls = -m(mIndex)*T +dt_r;
error_tls = abs(dt_ls-delt_s);

fprintf('Finer %ct_s estimate error: %g %cs. \n',948-32,error_tls*1e6, 956 )
