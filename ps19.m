h= 400;
v= -sqrt(398600.4/(6378+h));
z= [0:0.01:pi/2];

rdot  = sin(z)*v;
rddot = v^2/h*(cos(z)).^3;
figure(1)
plot(z,rdot);
figure (2)
plot(z,rddot);


%%

fc = 10;
vs = 3*10^6; %km/s

for i = 1:length(z)

    fr(i) = fc/(1+rdot(i)/vs);
    fd(i) = fr(i)-fc;

end

figure(3)
plot (z,fd*1000000);
