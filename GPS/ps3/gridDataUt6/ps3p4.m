clear all;
close all;
clc;

load("channel.mat")
M= channel';  % transpose to get in column
iidum = find(M(:,14)==7 & M(:,10) == 1);
txid07 = M(iidum,:);

% Extract L1 and L2 data
jjdum=find(txid07(:,13) ==0);
kkdum=find(txid07(:,13) == 2);
L1=txid07(jjdum,:);
L2=txid07(kkdum,:);

% Match L1 and L2 that were taken simultaneously and extract their
% pseudoranges and Carrier to Noise ratio
p  = []; % Initiate Pseudorange container
cn = []; % Initiate Carrier to Noise Ratio container
for i= 1:length(L1)
    idx = find( L2(:,2) == L1(i,2));
    if ~isempty(idx)
        p  = [p;L1(i,8), L2(idx,8)];
        cn = [cn;L1(i,9), L2(idx,9)];
    end
end




%%
f_L1 = 1575.42 * 10^6; % Hz
f_L2 = 1227.6 * 10^6; % Hz

% Code Derived Ionospheric Delay
I_L1_p = f_L2^2/(f_L1^2-f_L2^2)*(p(:,2)-p(:,1));
p_IF = f_L1^2*p(:,1)/(f_L1^2-f_L2^2)-f_L2^2*p(:,2)/(f_L1^2-f_L2^2);

% Carrier Derived Ionospheric Delay
I_L1_cn = f_L2^2/(f_L1^2-f_L2^2)*(cn(:,1)-cn(:,2));
%% Plot

t= L2(:,2)-L2(1,2);% Create Time vecotr based on measurements from channel.mat
plot(t,I_L1_p,t,I_L1_cn);
plot(t,I_L1_cn);
hold on,


title('Code-Derived Ionospheric Delay');
xlabel('Time (seconds)');
ylabel('Ionospheric Delay (meters)');

