clear all
clc
%%% Problem set 2

%% 1a

syms f
min = -2;
max = 2;
z=double(int(sinc(f)^2,min,max));

% delt =0.01;
% f = min:delt:max;
% Px = 0;
% for i = 1:length(f)
%    Px = Px+delt*sinc(f(i))^2;
% end

delt =0.001;
f = min:delt:max;
Px = 0;
for i = 2:length(f)
    if isnan(sin(pi*f(i))/(pi*f(i)))
        Px = Px+delt;
        graph(i)=Px;
    else
        Px = Px +delt*(sin(pi*f(i))/(pi*f(i)))^2;
        graph(i)=Px;
    end
end
fprintf('\n');
fprintf('%.3f%%',Px*100)
fprintf('\n');

%% 1 b

Rx = sinc(f);
plot(f,Rx), hold on,
title('Autocorrelation, Rx'),
xlabel('Frequency Hz')

