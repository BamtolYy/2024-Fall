close all; clear all; clc;
%% Generate Continuous Filters
Bn = 10;    % Hz

% 1st Order
k1real = 4*Bn;
H1real = tf(k1real,[1 k1real]);
% 2nd Order
k2real = 8/3*Bn;
a2real = k2real/2;
H2real = tf([k2real a2real*k2real],[1 k2real k2real*a2real]);
% 3rd Order
a3real = 1.2*Bn;
b3real = a3real^2/2;
k3real = 2*a3real;
H3real = tf([k3real k3real*a3real k3real*b3real],[1 k3real k3real*a3real k3real*b3real]);

%% Set up and Discretize filter loops at different orders
Bn   = 10;    % Hz
T    = [1,10,20,40]./1000;

for ii = 1:length(T)
    t = [0:T(ii):1];
    ustep = ones(length(t),1);

    slope = 1;
    uramp = slope*t';
    u = [ustep,uramp];

    Az   = tf([1 1],[2 0],T(ii));
    NCO  = tf(1,[1 0]);
    NCOd = c2d(NCO,T(ii));
    % 1st Order
    k1 = 4*Bn;
    D1 = tf(k1,1);
    Dd1 = c2d(D1,T(ii));
    H1 = feedback(Az*NCOd*Dd1,1);
    % 2nd Order
    k2 = 8/3*Bn;
    a2 = k2/2;
    D2 = tf(k2*[1 a2],[1 0]);
    Dd2 = c2d(D2,T(ii));
    H2 = feedback(Az*NCOd*Dd2,1);
    % 3rd Order
    a3 = 1.2*Bn;
    b3 = a3^2/2;
    k3 = 2*a3;
    D3 = tf(k3*[1 a3 b3],[1 0 0]);
    Dd3 = c2d(D3,T(ii));
    H3 = feedback(Az*NCOd*Dd3,1);
    figure,
    opts = bodeoptions;
    opts.XLimMode = 'manual';
    opts.XLim = {[1 6*10^1]};
    bodeplot(H1,H2,H3,opts)
    hold on
    bode(H1real,'b--',H2real,'r--',H3real,'y--')
    title(['Frequnecy Response with Discretization Interval of ',num2str(T(ii)*1000), ' ms'])
    legend('1st Order Loop Filter','2nd Order Loop Filter','3rd Order Loop Filter','1st order Cont.','2nd order Cont.','3rd order Cont.','Location','southwest')
    switch ii
        case 1
            ystep1(:,1)=lsim(H1,u(:,1),t);
            ystep1(:,2)=lsim(H2,u(:,1),t);
            ystep1(:,3)=lsim(H3,u(:,1),t);
            yramp1(:,1)=lsim(H1,u(:,2),t);
            yramp1(:,2)=lsim(H2,u(:,2),t);
            yramp1(:,3)=lsim(H3,u(:,2),t);
            %% Calculate Actual Loop Noise
            % Calculate Bn_act
            walias = pi/T(ii);
            wvec = [0:10000]'*(walias/10000);
            [magvec,phsvec] = bode(H1,wvec);
            magvec = magvec(:);
            Bn_act1(1) = sum(magvec.^2)*mean(diff(wvec))/(2*pi*(magvec(1,1)^2));
            [magvec2,phsvec2] = bode(H2,wvec);
            magvec2 = magvec2(:);
            Bn_act1(2) = sum(magvec2.^2)*mean(diff(wvec))/(2*pi*(magvec2(1,1)^2));
            [magvec3,phsvec3] = bode(H3,wvec);
            magvec3 = magvec3(:);
            Bn_act1(3) = sum(magvec3.^2)*mean(diff(wvec))/(2*pi*(magvec3(1,1)^2));
        case 2
            ystep10(:,1)=lsim(H1,u(:,1),t);
            ystep10(:,2)=lsim(H2,u(:,1),t);
            ystep10(:,3)=lsim(H3,u(:,1),t);
            yramp10(:,1)=lsim(H1,u(:,2),t);
            yramp10(:,2)=lsim(H2,u(:,2),t);
            yramp10(:,3)=lsim(H3,u(:,2),t);
            %% Calculate Actual Loop Noise
            % Calculate Bn_act
            walias = pi/T(ii);
            wvec = [0:10000]'*(walias/10000);
            [magvec,phsvec] = bode(H1,wvec);
            magvec = magvec(:);
            Bn_act10(1) = sum(magvec.^2)*mean(diff(wvec))/(2*pi*(magvec(1,1)^2));
            [magvec2,phsvec2] = bode(H2,wvec);
            magvec2 = magvec2(:);
            Bn_act10(2) = sum(magvec2.^2)*mean(diff(wvec))/(2*pi*(magvec2(1,1)^2));
            [magvec3,phsvec3] = bode(H3,wvec);
            magvec3 = magvec3(:);
            Bn_act10(3) = sum(magvec3.^2)*mean(diff(wvec))/(2*pi*(magvec3(1,1)^2));
        case 3
            ystep20(:,1)=lsim(H1,u(:,1),t);
            ystep20(:,2)=lsim(H2,u(:,1),t);
            ystep20(:,3)=lsim(H3,u(:,1),t);
            yramp20(:,1)=lsim(H1,u(:,2),t);
            yramp20(:,2)=lsim(H2,u(:,2),t);
            yramp20(:,3)=lsim(H3,u(:,2),t);
            %% Calculate Actual Loop Noise Bandwidth
            % Calculate Bn_act
            walias = pi/T(ii);
            wvec = [0:10000]'*(walias/10000);
            [magvec,phsvec] = bode(H1,wvec);
            magvec = magvec(:);
            Bn_act20(1) = sum(magvec.^2)*mean(diff(wvec))/(2*pi*(magvec(1,1)^2));
            [magvec2,phsvec2] = bode(H2,wvec);
            magvec2 = magvec2(:);
            Bn_act20(2) = sum(magvec2.^2)*mean(diff(wvec))/(2*pi*(magvec2(1,1)^2));
            [magvec3,phsvec3] = bode(H3,wvec);
            magvec3 = magvec3(:);
            Bn_act20(3) = sum(magvec3.^2)*mean(diff(wvec))/(2*pi*(magvec3(1,1)^2));
        case 4
            ystep40(:,1)=lsim(H1,u(:,1),t);
            ystep40(:,2)=lsim(H2,u(:,1),t);
            ystep40(:,3)=lsim(H3,u(:,1),t);
            yramp40(:,1)=lsim(H1,u(:,2),t);
            yramp40(:,2)=lsim(H2,u(:,2),t);
            yramp40(:,3)=lsim(H3,u(:,2),t);
            %% Calculate Actual Loop Noise
            % Calculate Bn_act
            walias = pi/T(ii);
            wvec = [0:10000]'*(walias/10000);
            [magvec,phsvec] = bode(H1,wvec);
            magvec = magvec(:);
            Bn_act40(1) = sum(magvec.^2)*mean(diff(wvec))/(2*pi*(magvec(1,1)^2));
            [magvec2,phsvec2] = bode(H2,wvec);
            magvec2 = magvec2(:);
            Bn_act40(2) = sum(magvec2.^2)*mean(diff(wvec))/(2*pi*(magvec2(1,1)^2));
            [magvec3,phsvec3] = bode(H3,wvec);
            magvec3 = magvec3(:);
            Bn_act40(3) = sum(magvec3.^2)*mean(diff(wvec))/(2*pi*(magvec3(1,1)^2));
    end
end

%% Plot step response at different dicretization intervals
t1=0:T(1):1;
t10=0:T(2):1;
t20=0:T(3):1;
t40=0:T(4):1;

figure
title(['Step Response of First Order Loop Filter with Varying Discretization'])
hold on,
plot(t40,ustep,'--',t1,ystep1(:,1),t10,ystep10(:,1),t20,ystep20(:,1),t40,ystep40(:,1))
legend('Reference Input','1 ms','10 ms','20 ms', '40 ms')
xlabel('Time (second)')
ylabel('y(t)')

figure
title(['Step Response of 2nd Order Loop Filter with Varying Discretization'])
hold on,
plot(t40,ustep,'--',t1,ystep1(:,2),t10,ystep10(:,2),t20,ystep20(:,2),t40,ystep40(:,2))
legend('Reference Input','1 ms','10 ms','20 ms', '40 ms')
xlabel('Time (second)')
ylabel('y(t)')

figure
title(['Step Response of 3rd Order Loop Filter with Varying Discretization'])
hold on,
plot(t40,ustep,'--',t1,ystep1(:,3),t10,ystep10(:,3),t20,ystep20(:,3),t40,ystep40(:,3))
legend('Reference Input','1 ms','10 ms','20 ms', '40 ms')
xlabel('Time (second)')
ylabel('y(t)')

%% Plot ramp response at different dicretization intervals

figure
title(['Ramp Response of First Order Loop Filter with Varying Discretization'])
hold on,
plot(t40,uramp,'--',t1,yramp1(:,1),t10,yramp10(:,1),t20,yramp20(:,1),t40,yramp40(:,1))
legend('Reference Input','1 ms','10 ms','20 ms', '40 ms')
xlabel('Time (second)')
ylabel('y(t)')

figure
title(['Ramp Response of 2nd Order Loop Filter with Varying Discretization'])
hold on,
plot(t40,uramp,'--',t1,yramp1(:,2),t10,yramp10(:,2),t20,yramp20(:,2),t40,yramp40(:,2))
legend('Reference Input','1 ms','10 ms','20 ms', '40 ms')
xlabel('Time (second)')
ylabel('y(t)')

figure
title(['Ramp Response of 3rd Order Loop Filter with Varying Discretization'])
hold on,
plot(t40,uramp,'--',t1,yramp1(:,3),t10,yramp10(:,3),t20,yramp20(:,3),t40,yramp40(:,3))
legend('Reference Input','1 ms','10 ms','20 ms', '40 ms')
xlabel('Time (second)')
ylabel('y(t)')


%% Good "Rule of Thumb" for BnT
disp('----------------------------------------------------------------------')
fprintf(['A good rule of thumb for the maximum value of the product BnT would be\n' ...
    '0.01, since the phase and gain of the system begin to be distorted significantly\n' ...
    'for discretization intervals greater than 10. Especially, the phase of the system\n' ...
    'with 20 ms and 40 ms system drop at a much higher rate than the continuous system.\n' ...
    'Also, 40 ms system shows a much higher resonant frequency than the origina system.\n' ...
    'In contrast, the 1 ms system exhibits a relatively similar drop in phase.\n'])
disp('----------------------------------------------------------------------')

%% Actual Loop Noise
fprintf(['Actual Loop Noise Bandwidths\n'])
disp(['             |    1 ms     |     10 ms     |     20 ms     |     40 ms     |'])
disp(['First Order:    ', num2str(Bn_act1(1)),'Hz       ', num2str(Bn_act10(1)),'Hz       ',num2str(Bn_act20(1)),'Hz       ',num2str(Bn_act40(1)),'Hz    '])
disp(['Second Order:   ', num2str(Bn_act1(2)),'Hz       ', num2str(Bn_act10(2)),'Hz     ',num2str(Bn_act20(2)),'Hz       ',num2str(Bn_act40(2)),'Hz    '])
disp(['Third Order:    ', num2str(Bn_act1(3)),'Hz       ', num2str(Bn_act10(3)),'Hz     ',num2str(Bn_act20(3)),'Hz       ',num2str(Bn_act40(3)),'Hz    '])
disp('\n----------------------------------------------------------------------')

%% What is the code doing and why does it accurately estimates the actual bandwidth?

fprintf(['We are essentially integrating power of the signal over all the frequencies upto Nyquist frequency.\n' ...
    'It needs to be divided by 2*pi, since we caclulated the value with radians per second and the unit \n' ...
    'we want is Hz. response of the system.It is then normalized by the zero frequency.\n' ...
    'Then, this gives us the bandwidth of the noise. This calculation accurately estimates the actual loop noise bandwidth beccause it \n' ...
    'takes the powers of the actual signal and integrating them over the effective frequencies as the definition \n' ...
    'of the Noise bandwidth states. The normalization ensures that the bandwidth is normalized to the DC gain of the system, and \n' ...
    'eliminate any DC gain effect of the system at 0 frequency (Isolate noise effect only in the frequnecy response).'])
disp('\n----------------------------------------------------------------------')

%% Explain foh and tustin method
fprintf(['The first order hold method uses "triangular pulse" as opposed to "rectangular pulse" in the zero order hold discretization method.\n' ...
    'It will not give a constant value from one sample to another as in the zoh method, foh will give a linear value from one sample to another.\n' ...
    'The tustin method is a bilinear transformation that maps the function in s-plane to z plane using z = exp(sT), where T is the discretization interval.'])
disp('\n----------------------------------------------------------------------')

%% Experiment with foh
% First order
figure
for jj =1 : 4
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'foh');
    Ddfoh = c2d(D1,T(jj),'foh');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    waliasfoh = pi/T(jj);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvecfoh,phsvecfoh] = bode(Hfoh,wvecfoh);
    mag1(:,jj) = magvecfoh;
    pha1(:,jj) = phsvecfoh;
    magvecfoh = magvecfoh(:);
    Bn_actfoh(jj) = sum(magvecfoh.^2)*mean(diff(wvecfoh))/(2*pi*(magvecfoh(1,1)^2));
    if jj <4
        opts = bodeoptions;
        opts.XLimMode = 'manual';
        opts.XLim = {[1 6*10^1]};
        bodeplot(Hfoh,opts)
        hold on
    elseif jj == 4
        bode(Hfoh)
        hold on
        bode(H1real,'r--')
        legend('1 ms','10 ms','20 ms','40 ms', 'Continous','Location','southwest')
        sgtitle('1st Order Filter FOH Method')
    end
end

% 2nd order
figure
for jj =1 : 4
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'foh');
    Dd2foh = c2d(D2,T(jj),'foh');
    Hfoh2 = feedback(Azfoh*NCOdfoh*Dd2foh,1);
    waliasfoh = pi/T(jj);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvecfoh,phsvecfoh] = bode(Hfoh2,wvecfoh);
    mag2(:,jj) = magvecfoh;
    pha2(:,jj) = phsvecfoh;
    magvecfoh = magvecfoh(:);
    Bn_actfoh2(jj) = sum(magvecfoh.^2)*mean(diff(wvecfoh))/(2*pi*(magvecfoh(1,1)^2));
    if jj <4
        opts = bodeoptions;
        opts.XLimMode = 'manual';
        opts.XLim = {[1 6*10^1]};
        bodeplot(Hfoh2,opts)
        hold on
    elseif jj == 4
        bode(Hfoh2)
        hold on
        bode(H2real,'r--')
        legend('1 ms','10 ms','20 ms','40 ms', 'Continous','Location','southwest')
        sgtitle('2nd Order Filter FOH Method')
    end
end

% 3rd order
figure
for jj =1 : 4
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'foh');
    Dd3foh = c2d(D3,T(jj),'foh');
    Hfoh3 = feedback(Azfoh*NCOdfoh*Dd3foh,1);
    waliasfoh = pi/T(jj);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvecfoh,phsvecfoh] = bode(Hfoh3,wvecfoh);
    mag3(:,jj) = magvecfoh;
    pha3(:,jj) = phsvecfoh;
    magvecfoh = magvecfoh(:);
    Bn_actfoh3(jj) = sum(magvecfoh.^2)*mean(diff(wvecfoh))/(2*pi*(magvecfoh(1,1)^2));
    if jj <4
        opts = bodeoptions;
        opts.XLimMode = 'manual';
        opts.XLim = {[1 6*10^1]};
        bodeplot(Hfoh3,opts)
        hold on
    elseif jj == 4
        bode(Hfoh3)
        hold on
        bode(H3real,'r--')
        legend('1 ms','10 ms','20 ms','40 ms', 'Continous','Location','southwest')
        sgtitle('3rd Order Filter FOH Method')
    end
end
disp('----------------------------------------------------------------------')

fprintf(['Actual Loop Noise Bandwidths from First order discretization method\n'])
disp(['             |    1 ms     |     10 ms     |     20 ms     |     40 ms     |'])
disp(['First Order:    ', num2str(Bn_actfoh(1)),'Hz       ', num2str(Bn_actfoh(2)),'Hz        ',num2str(Bn_actfoh(3)),'Hz       ',num2str(Bn_actfoh(4)),'Hz    '])
disp(['Second Order:   ', num2str(Bn_actfoh2(1)),'Hz       ', num2str(Bn_actfoh2(2)),'Hz        ',num2str(Bn_actfoh2(3)),'Hz       ',num2str(Bn_actfoh2(4)),'Hz    '])
disp(['Third Order:    ', num2str(Bn_actfoh3(1)),'Hz       ', num2str(Bn_actfoh3(2)),'Hz        ',num2str(Bn_actfoh3(3)),'Hz       ',num2str(Bn_actfoh3(4)),'Hz    '])
disp('\n----------------------------------------------------------------------')

%% Experiment with Tustin

% 1st order
figure
for jj =1 : 4
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'tustin');
    Ddfoh = c2d(D1,T(jj),'tustin');
    Htustin = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    waliasfoh = pi/T(jj);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvectustin,phsvectustin] = bode(Htustin,wvecfoh);
    magvectustin = magvectustin(:);
    Bn_acttustin(jj) = sum(magvectustin.^2)*mean(diff(wvecfoh))/(2*pi*(magvectustin(1,1)^2));
    if jj <4
        opts = bodeoptions;
        opts.XLimMode = 'manual';
        opts.XLim = {[1 6*10^1]};
        bodeplot(Htustin,opts)

        hold on
    elseif jj == 4
        bode(Htustin)
        hold on
        bode(H1real,'r--')
        legend('1 ms','10 ms','20 ms','40 ms', 'Continous','Location','southwest')
        sgtitle('1st Order Filter Tustin Method')
    end
end

% 2nd order
figure
for jj =1 : 4
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'tustin');
    Ddfoh2 = c2d(D2,T(jj),'tustin');
    Htustin2 = feedback(Azfoh*NCOdfoh*Ddfoh2,1);
    waliasfoh = pi/T(jj);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvectustin,phsvectustin] = bode(Htustin2,wvecfoh);
    magvectustin2 = magvectustin(:);
    Bn_acttustin2(jj) = sum(magvectustin2.^2)*mean(diff(wvecfoh))/(2*pi*(magvectustin2(1,1)^2));
    if jj <4
        opts = bodeoptions;
        opts.XLimMode = 'manual';
        opts.XLim = {[1 6*10^1]};
        bodeplot(Htustin2,opts)
        hold on
    elseif jj == 4
        bode(Htustin2)
        hold on
        bode(H2real,'r--')
        legend('1 ms','10 ms','20 ms','40 ms', 'Continous','Location','southwest')
        sgtitle('2nd Order Filter Tustin Method')
    end
end

% 3rd order
figure
for jj =1 : 4
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'tustin','');
    Ddfoh3 = c2d(D3,T(jj),'tustin');
    Htustin3 = feedback(Azfoh*NCOdfoh*Ddfoh3,1);
    waliasfoh = pi/T(jj);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvectustin,phsvectustin] = bode(Htustin3,wvecfoh);
    magvectustin3 = magvectustin(:);
    Bn_acttustin3(jj) = sum(magvectustin3.^2)*mean(diff(wvecfoh))/(2*pi*(magvectustin3(1,1)^2));
    if jj <4
        opts = bodeoptions;
        opts.XLimMode = 'manual';
        opts.XLim = {[1 6*10^1]};
        bodeplot(Htustin3,opts)
        hold on
    elseif jj == 4
        bode(Htustin3)
        hold on
        bode(H3real,'r--')
        legend('1 ms','10 ms','20 ms','40 ms', 'Continous','Location','southwest')
        sgtitle('3rd Order Filter Tustin Method')
    end
end

fprintf(['Actual Loop Noise Bandwidths from Tustin discretization method\n'])
disp(['             |    1 ms     |     10 ms     |     20 ms     |     40 ms     |'])
disp(['First Order:    ', num2str(Bn_acttustin(1)),'Hz       ', num2str(Bn_acttustin(2)),'Hz        ',num2str(Bn_acttustin(3)),'Hz       ',num2str(Bn_actfoh(4)),'Hz    '])
disp(['Second Order:   ', num2str(Bn_acttustin2(1)),'Hz       ', num2str(Bn_acttustin2(2)),'Hz        ',num2str(Bn_acttustin2(3)),'Hz       ',num2str(Bn_actfoh(4)),'Hz    '])
disp(['Third Order:    ', num2str(Bn_acttustin3(1)),'Hz       ', num2str(Bn_acttustin3(2)),'Hz        ',num2str(Bn_acttustin3(3)),'Hz       ',num2str(Bn_actfoh(4)),'Hz    '])
disp('\n----------------------------------------------------------------------')

%% Step Response of foh and Tustin
% FOH
figure
for jj =1 : 4
    t=0:T(jj):1;
    ustep = ones(length(t),1);
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'foh');
    Ddfoh = c2d(D1,T(jj),'foh');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    if jj <4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,ustep,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Reference','Location','northeast')
        title('1st Order Filter FOH Method')
    end
end

% 2nd order
figure
for jj =1 : 4
    t=0:T(jj):1;
    ustep = ones(length(t),1);
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'foh');
    Ddfoh2 = c2d(D2,T(jj),'foh');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh2,1);
    if jj <4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,ustep,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Reference','Location','northeast')
        title('2nd Order Filter FOH Method')
    end
end

% 3rd order
figure
for jj =1 : 4
    t=0:T(jj):1;
    ustep = ones(length(t),1);
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'foh');
    Ddfoh3 = c2d(D3,T(jj),'foh');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh3,1);
    if jj <4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,ustep,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Reference','Location','northeast')
        title('3rd Order Filter FOH Method')
    end
end

% Tustin
figure
for jj =1 : 4
    t=0:T(jj):1;
    ustep = ones(length(t),1);
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'tustin');
    Ddfoh = c2d(D1,T(jj),'tustin');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    if jj <4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,ustep,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Reference','Location','northeast')
        title('1st Order Filter Tustin Method')
    end
end

% 2nd order
figure
for jj =1 : 4
    t=0:T(jj):1;
    ustep = ones(length(t),1);
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'tustin');
    Ddfoh2 = c2d(D2,T(jj),'tustin');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh2,1);
    if jj <4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,ustep,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Reference','Location','northeast')
        title('2nd Order Filter Tustin Method')
    end
end

% 3rd order
figure
for jj =1 : 4
    t=0:T(jj):1;
    ustep = ones(length(t),1);
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'tustin');
    Ddfoh3 = c2d(D3,T(jj),'tustin');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh3,1);
    if jj <4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,ustep,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,ustep,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Reference','Location','northeast')
        title('3rd Order Filter Tustin Method')
    end
end

%% Ramp Response of foh and tustin
% FOH
figure
for jj =1 : 4
    t=0:T(jj):1;
    uramp = slope*t';
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'foh');
    Ddfoh = c2d(D1,T(jj),'foh');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    if jj <4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,uramp,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Reference','Location','northeast')
        title('1st Order Filter FOH Method')
    end
end

% 2nd order
figure
for jj =1 : 4
    t=0:T(jj):1;
    uramp = slope*t';
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'foh');
    Ddfoh2 = c2d(D2,T(jj),'foh');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh2,1);
    if jj <4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,uramp,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Reference','Location','northeast')
        title('2nd Order Filter FOH Method')
    end
end

% 3rd order
figure
for jj =1 : 4
    t=0:T(jj):1;
    uramp = slope*t';
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'foh');
    Ddfoh3 = c2d(D3,T(jj),'foh');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh3,1);
    if jj <4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,uramp,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Reference','Location','northeast')
        title('3rd Order Filter FOH Method')
    end
end

% Tustin
figure
for jj =1 : 4
    t=0:T(jj):1;
    uramp = slope*t';
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'tustin');
    Ddfoh = c2d(D1,T(jj),'tustin');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    if jj <4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,uramp,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Reference','Location','northeast')
        title('1st Order Filter Tustin Method')
    end
end

% 2nd order
figure
for jj =1 : 4
    t=0:T(jj):1;
    uramp = slope*t';
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'tustin');
    Ddfoh2 = c2d(D2,T(jj),'tustin');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh2,1);
    if jj <4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,uramp,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Reference','Location','northeast')
        title('2nd Order Filter Tustin Method')
    end
end

% 3rd order
figure
for jj =1 : 4
    t=0:T(jj):1;
    uramp = slope*t';
    Azfoh   = tf([1 1],[2 0],T(jj));
    NCOdfoh = c2d(NCO,T(jj),'tustin');
    Ddfoh3 = c2d(D3,T(jj),'tustin');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh3,1);
    if jj <4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        hold on
    elseif jj == 4
        yfohstep = lsim(Hfoh,uramp,t);
        plot(t,yfohstep)
        ylabel('Output y(t)')
        xlabel('Time (s)')
        hold on
        plot(t,uramp,'--')
        legend('1 ms','10 ms','20 ms','40 ms','Location','northeast')
        title('3rd Order Filter Tustin Method')
    end
end
fprintf(['Both discretization methods provide better frequency and step \n' ...
    'responses that more accurately tracks the continous system. Their \n' ...
    'actual bandwidth is much closer to Bn, the target bandwidth than that\n' ...
    'of zero order hold discretization method.\n' ...
    'A loop filter designed in the z-domain could not have an arbitrarily\n' ...
    'high BnT. It will have limits, since poles outside the unit circle\]n' ...
    'will be unstable. However, the maximum of BnT could be higher with \n' ...
    'z-domain design, since z-domain design is free of discretization discrepancy\n' ...
    'from continous to discretization. A designer can fine tune the controller to \n' ...
    'have a better BnT than BnT from the rule of thumb'])