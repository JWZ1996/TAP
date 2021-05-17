pidr
close all
Ts = 0.001;
fs = 1/Ts;
load('time.mat');
load('Yzad.mat');
[y1,t1]=lsim(G12,Yzad(1,:),time);
[y2,t2]=lsim(G21,Yzad(2,:),time);

figure;
plot(t2, y2, t2, Yzad(2,:));
title('Wyjście 2 - temperatura')
xlabel('Czas symulacji [min]')
ylabel('T [K]')

figure;
plot(t1, y1, t1, Yzad(1,:));
title('Wyjście 1 - stężenie substancji A')
xlabel('Czas symulacji [min]')
ylabel('Ca [kmol/m3]') 

