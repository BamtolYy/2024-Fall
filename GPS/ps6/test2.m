close all; clear all; clc;
%% a) Generate 5ms worth of signal and its replica

% ----x(j)
% LFSR Parameters:
nStages      = 10;
ciVec1       = [10, 3]';
ciVec2       = [10, 9, 8, 6, 3, 2,]';
a0Vec1       = ones(nStages,1);
a0Vec2       = ones(nStages,1);
G2tab = [2,6];
% Oversampling Parameters:
fs = 40e6/7;                % Sampling Rate
Tc = 1e-3/1023;             % Chip interval in seconds
T  = 1/fs;                  % Sampling time interval in seconds
delChip = T/Tc;             % Sampling interval in chips
Np = 2^nStages - 1;         % Period of the sequence in chips
Ns = ceil(0.005/T);         % Number of Samples
ts    = 1e-3;               % True delay in seconds
delOffset = ts/Tc;          % Offset number in chips
xj = zeros(Ns,1);
lj = zeros(Ns,1);
% Signal Parameters
sigman = sqrt(1/2);
CN0   = 10^(45/10);
A     = sqrt(4*sigman^2*T*CN0);
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
    GoldSeqOS = oversampleSpreadingCode(GoldSeq,delChip,delOffset,Ns,Np);
    xj(:,j) = GoldSeqOS;
end

% Generate Noise
a      = sigman.*randn(Ns,1);
b      = sigman.*randn(Ns,1);
nj     = a +1i*b;
xj    = A*xj +nj;

% ----l(j)
delt_s        = 0.25e-6;                    % seconds
delOffSetReplica  = (ts + delt_s)/Tc;            % Estimate Delay in chips
GoldSeq       = [];
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
    GoldSeqOS = oversampleSpreadingCode(GoldSeq,delChip,delOffSetReplica,Ns,Np);
    lj(:,j) = GoldSeqOS;
end



%% b) Compute Rtilde_xl(m)
m        =  -10:10;
Rxltilde = zeros(length(m),1);
jk       = 1024;
Nk       = floor(Np*4/delChip);
for ii = 1:length(m)
    Rxltilde(ii) = 1/Nk*sum(xj(jk:jk+Nk-1).*lj(jk-m(ii):jk+Nk-1-m(ii)));
end

%% c) Nearest-Sample Estimator

[~,mIndex] = max(abs(Rxltilde));
delt_ns  = m(mIndex)*T;
error_ts = abs(delt_ns - delt_s);
disp(['3c) and 3d)'])
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
dt_ls = m(mIndex)*T +dt_r;
error_tls = abs(dt_ls-delt_s);
fprintf('Finer %ct_s estimate error: %g %cs. \n',948-32,error_tls*1e6, 956 )
disp(['-----------------------------------------------------------------'])
%% e) Repeat 1000 times
disp(['3e)'])
ensemble = 1000;
for pp = 1:ensemble
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
        GoldSeqOS = oversampleSpreadingCode(GoldSeq,delChip,delOffset,Ns,Np);
        xj(:,j) = GoldSeqOS;
    end
    xj    = A*xj +nj;


    % ----l(j)
    delt_s        = 0.25e-6;
    shiftReplica = floor((ts + delt_s)/T);            % Estimate Delay in samples
    for j = 1
        [GoldSeq] = generateGoldLfsrSequenceCA(nStages,ciVec1,ciVec2,a0Vec1,...
            a0Vec2,G2tab(j,:));
        % Make code +1/-1 not +1/0
        GoldSeq = 2*GoldSeq - 1;
        GoldSeqOS = oversampleSpreadingCode(GoldSeq,delChip,delOffSetReplica,Ns,Np);
        lj(:,j) = GoldSeqOS;
    end



    m        =  -10:1:10;
    Rxltilde = zeros(length(m),1);
    jk = 1024;
    for ii = 1:length(m)
        Nk       = floor(Np*4/delChip);
        Rxltilde(ii) = 1/Nk*sum(xj(jk:jk+Nk-1).*lj(jk-m(ii):jk+Nk-1-m(ii)));
    end
    [~,mIndex] = max(abs(Rxltilde));
    delt_ns  = m(mIndex)*T;
    sampleEstimatorError(pp) = abs(delt_ns - delt_s);
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
    dt_ls = m(mIndex)*T +dt_r;
    leastSquaresEstimatorError(pp) = abs(dt_ls-delt_s);
end

% Calculate mean squared error
sigmahat_ns2 = mean(sampleEstimatorError.^2);
sigmahat_ls2 = mean(leastSquaresEstimatorError.^2);
RMSE_ns = sqrt(sigmahat_ns2) * physconst('LightSpeed');
RMSE_ls = sqrt(sigmahat_ls2) * physconst('LightSpeed');

fprintf('RMSE (Sample-Level Estimator): %g meters\n', RMSE_ns);
fprintf('RMSE (Least-Squares Estimator): %g meters\n', RMSE_ls);

fprintf(['The difference between the two is significant. They differ by more\n' ...
    'than 17 meters.\n'])
disp(['-----------------------------------------------------------------'])
%% f) CRLB
disp(['3f)'])
Ta = Nk*T;
betams2 = 1.275*10^(13);       % (rad/s)^2
CRLB = (2*betams2*Ta*CN0)^(-1);
RMSE_distance = sqrt(CRLB)*physconst('LightSpeed');
fprintf('CRLB: %g meters\n', CRLB);
fprintf('RMSE distance L1: %g meters\n', RMSE_distance);
fprintf(['The RMSE from the least squares estimator is very close to the RMSE distance.\n' ...
    'It is only about 0.2 meters greater than RMSE distance.' ...
    'In contrast, the sample-level \nestimator is much greater than ' ...
    'the RMSE distance.\n'])
disp(['-----------------------------------------------------------------'])

%% g) New CRLB
disp(['3g)'])
betams2New = (2*pi*fs)^2/12;    % (rad/s)^2
CRLBNew = (2*betams2New*Ta*CN0)^(-1);
RMSE_distanceNew = sqrt(CRLBNew)*physconst('LightSpeed');
fprintf('New CRLB: %g meters\n', CRLBNew);
fprintf('New RMSE distance: %g meters\n', RMSE_distanceNew);
fprintf(['This new signal has much lower CRLB than the L1 signal. ' ...
    'It is about 3 to 4 meters \nless than that of L1 signal. This means ' ...
    'it can provide better position estimates, \nsince the theoretical lower' ...
    'limit of the variance is smaller. However, it may be more suseptible\n' ...
    'to noise because the bandwidth is increased.\n'])
disp(['-----------------------------------------------------------------'])
