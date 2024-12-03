close all; clear all; clc;
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
    title(['Frequnecy Response with Discretization Interval of ',num2str(T(ii)*1000), ' ms'])
    legend('1st Order Loop Filter','2nd Order Loop Filter','3rd Order Loop Filter','Location','southwest')
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
disp('----------------------------------------------------------------------')

%% What is the code doing and why does it accurately estimates the actual bandwidth?

fprintf(['We are essentially integrating power of the signal over all the frequencies upto Nyquist frequency.\n' ...
    'It needs to be divided by 2*pi, since we caclulated the value with radians per second and the unit \n' ...
    'we want is Hz. response of the system.It is then normalized by the zero frequency.\n' ...
    'Then, this gives us the bandwidth of the noise. This calculation accurately estimates the actual loop noise bandwidth beccause it \n' ...
    'takes the powers of the actual signal and integrating them over the effective frequencies as the definition \n' ...
    'of the Noise bandwidth states. The normalization ensures that the bandwidth is normalized to the DC gain of the system, and \n' ...
    'eliminate any DC gain effect of the system at 0 frequency (Isolate noise effect only in the frequnecy response).'])
disp('----------------------------------------------------------------------')

%% Explain foh and tustin method
fprintf(['The first order hold method uses "triangular pulse" as opposed to "rectangular pulse" in the zero order hold discretization method.\n' ...
    'It will not give a constant value from one sample to another as in the zoh method, foh will give a linear value from one sample to another.\n' ...
    'The tustin method is a bilinear transformation that maps the function in s-plane to z plane using z = exp(sT), where T is the discretization interval.'])
disp('----------------------------------------------------------------------')

%% Experiment with foh
% 20ms
for jj =1 : 3
    Azfoh   = tf([1 1],[2 0],T(3));
    NCOdfoh = c2d(NCO,T(3),'foh');
    Ddfoh = c2d(D1,T(3),'foh');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    waliasfoh = pi/T(3);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvecfoh,phsvecfoh] = bode(Hfoh,wvecfoh);
    magvecfoh = magvecfoh(:);
    Bn_act20foh(jj) = sum(magvecfoh.^2)*mean(diff(wvecfoh))/(2*pi*(magvecfoh(1,1)^2));
end
% 10ms
for jj =1 : 3
    Azfoh   = tf([1 1],[2 0],T(2));
    NCOdfoh = c2d(NCO,T(2),'foh');
    Ddfoh = c2d(D1,T(2),'foh');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    waliasfoh = pi/T(2);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvecfoh,phsvecfoh] = bode(Hfoh,wvecfoh);
    magvecfoh = magvecfoh(:);
    Bn_act10foh(jj) = sum(magvecfoh.^2)*mean(diff(wvecfoh))/(2*pi*(magvecfoh(1,1)^2));
end
% 1 ms
for jj =1 : 3
    Azfoh   = tf([1 1],[2 0],T(1));
    NCOdfoh = c2d(NCO,T(1),'foh');
    Ddfoh = c2d(D1,T(1),'foh');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    waliasfoh = pi/T(1);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvecfoh,phsvecfoh] = bode(Hfoh,wvecfoh);
    magvecfoh = magvecfoh(:);
    Bn_act1foh(jj) = sum(magvecfoh.^2)*mean(diff(wvecfoh))/(2*pi*(magvecfoh(1,1)^2));
end
% 40 ms
for jj =1 : 3
    Azfoh   = tf([1 1],[2 0],T(4));
    NCOdfoh = c2d(NCO,T(4),'foh');
    Ddfoh = c2d(D1,T(4),'foh');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    waliasfoh = pi/T(4);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvecfoh,phsvecfoh] = bode(Hfoh,wvecfoh);
    magvecfoh = magvecfoh(:);
    Bn_act40foh(jj) = sum(magvecfoh.^2)*mean(diff(wvecfoh))/(2*pi*(magvecfoh(1,1)^2));
end
fprintf(['Actual Loop Noise Bandwidths from First order discretization method\n'])
disp(['             |    1 ms     |     10 ms     |     20 ms     |     40 ms     |'])
disp(['First Order:    ', num2str(Bn_act1foh(1)),'Hz       ', num2str(Bn_act10foh(1)),'Hz        ',num2str(Bn_act20foh(1)),'Hz       ',num2str(Bn_act40foh(1)),'Hz    '])
disp(['Second Order:   ', num2str(Bn_act1foh(2)),'Hz       ', num2str(Bn_act10foh(2)),'Hz        ',num2str(Bn_act20foh(2)),'Hz       ',num2str(Bn_act40foh(2)),'Hz    '])
disp(['Third Order:    ', num2str(Bn_act1foh(3)),'Hz       ', num2str(Bn_act10foh(3)),'Hz        ',num2str(Bn_act20foh(3)),'Hz       ',num2str(Bn_act40foh(3)),'Hz    '])
disp('----------------------------------------------------------------------')

%% Experiment with Tustin

% 20ms
for jj =1 : 3
    Azfoh   = tf([1 1],[2 0],T(3));
    NCOdfoh = c2d(NCO,T(3),'tustin');
    Ddfoh = c2d(D1,T(3),'tustin');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    waliasfoh = pi/T(3);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvecfoh,phsvecfoh] = bode(Hfoh,wvecfoh);
    magvecfoh = magvecfoh(:);
    Bn_act20tustin(jj) = sum(magvecfoh.^2)*mean(diff(wvecfoh))/(2*pi*(magvecfoh(1,1)^2));
end
% 10ms
for jj =1 : 3
    Azfoh   = tf([1 1],[2 0],T(2));
    NCOdfoh = c2d(NCO,T(2),'tustin');
    Ddfoh = c2d(D1,T(2),'tustin');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    waliasfoh = pi/T(2);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvecfoh,phsvecfoh] = bode(Hfoh,wvecfoh);
    magvecfoh = magvecfoh(:);
    Bn_act10tustin(jj) = sum(magvecfoh.^2)*mean(diff(wvecfoh))/(2*pi*(magvecfoh(1,1)^2));
end
% 1 ms
for jj =1 : 3
    Azfoh   = tf([1 1],[2 0],T(1));
    NCOdfoh = c2d(NCO,T(1),'tustin');
    Ddfoh = c2d(D1,T(1),'tustin');
    Hfoh = feedback(Azfoh*NCOdfoh*Ddfoh,1);
    waliasfoh = pi/T(1);
    wvecfoh = [0:10000]'*(waliasfoh/10000);
    [magvecfoh,phsvecfoh] = bode(Hfoh,wvecfoh);
    magvecfoh = magvecfoh(:);
    Bn_act1tustin(jj) = sum(magvecfoh.^2)*mean(diff(wvecfoh))/(2*pi*(magvecfoh(1,1)^2));
end
% 40 ms
for jj =1 : 3
    Aztustin   = tf([1 1],[2 0],T(4));
    NCOdtustin = c2d(NCO,T(4),'tustin');
    Ddtustin = c2d(D1,T(4),'tustin');
    Htustin = feedback(Aztustin*NCOdtustin*Ddtustin,1);
    waliastustin = pi/T(4);
    wvectustin = [0:10000]'*(waliastustin/10000);
    [magvectustin,phsvectustin] = bode(Htustin,wvectustin);
    magvectustin = magvectustin(:);
    Bn_act40tustin(jj) = sum(magvectustin.^2)*mean(diff(wvectustin))/(2*pi*(magvectustin(1,1)^2));
end
fprintf(['Actual Loop Noise Bandwidths from First order discretization method\n'])
disp(['             |    1 ms     |     10 ms     |     20 ms     |     40 ms     |'])
disp(['First Order:    ', num2str(Bn_act1tustin(1)),'Hz       ', num2str(Bn_act10tustin(1)),'Hz        ',num2str(Bn_act20tustin(1)),'Hz       ',num2str(Bn_act40foh(1)),'Hz    '])
disp(['Second Order:   ', num2str(Bn_act1tustin(2)),'Hz       ', num2str(Bn_act10tustin(2)),'Hz        ',num2str(Bn_act20tustin(2)),'Hz       ',num2str(Bn_act40foh(2)),'Hz    '])
disp(['Third Order:    ', num2str(Bn_act1tustin(3)),'Hz       ', num2str(Bn_act10tustin(3)),'Hz        ',num2str(Bn_act20tustin(3)),'Hz       ',num2str(Bn_act40foh(3)),'Hz    '])
disp('----------------------------------------------------------------------')
