%% Parameters
N = 100;
Ta = 0.001; 
% True Alpha Parameters
f     = 20;
rho   = 1;
theta = pi/4;
fs    = 2*f;
T     = 1/fs;


%% Generate Sk
nI = randn(1,N);
nQ = randn(1,N);

k = 0:N-1;
nk = nI +1i*nQ;
Sk = rho*exp(1i*(2*pi*f*k*Ta+theta))+nk;

%% Estimate

% frequency
% FFT of the signal Sk
Nf = length(Sk);
S_fft = fft(Sk);

% Frequency vector
fv = (0:Nf-1) * (fs / Nf);  % Create a frequency vector corresponding to FFT bins

% Find the frequency with the maximum magnitude
[~, max_idx] = max(abs(S_fft));
f_ml = fv(max_idx); 