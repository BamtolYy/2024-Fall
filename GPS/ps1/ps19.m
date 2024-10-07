h= 400;
v= sqrt(398600.4/(6378+h));
z= [0:0.01:pi/2];

rdot  = sin(-z)*v;
rddot = -v^2/h*(cos(z)).^3;
figure(1)
sgtitle('Velocity and Acceleration of Line of Sight Vector of  Satellite to Observer')
subplot(121)
plot(z,rdot);
xlabel('Zenith Angle [rad^{-1}]' )
ylabel('Speed [km/s]')
subplot(122)
plot(z,rddot);
xlabel('Zenith Angle [rad^{-1}] ' )
ylabel('Acceleration [km^2/s]')



%%

fc = 10*10^9;
vs = physconst('lightspeed')/1000; %km/s

for i = 1:length(z)

    fr(i) = fc/(1+rdot(i)/vs);
    fd(i) = fr(i)-fc;
end

figure(3)
plot (z,fd/1000);
maxfd=max(abs(fd/1000))
