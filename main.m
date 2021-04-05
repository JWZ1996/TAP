close all;
clear all;

% ==== Plik gÅ‚owny =====
global g m K mol kmol cal min ro cp k E_R h a b ro cp k E_R h a b V Fin CAin Fc Tin Tcin Ca T;
consts
% === Symulacja obiektu - RK4 ===
step = 0.01;
Ca = 0.042953 * kmol/m^3;
T =  431.176 * K;

CAin_vect = [2.4, 2.2, 2, 1.8, 1.6, 1.4].*kmol;
Fc_vect = [18,17,16,15,14,13,12];

for iter = 1:1:length(CAin_vect)
    CAin = CAin_vect(iter);
    [y, t] = rk4(@dCa, @dT, Ca, T, step);
    txt = ['CAin = ',num2str(CAin_vect(iter))];
    figure(1)
    plot(t,y(1,:),'DisplayName',txt);
    title('Stê¿enie Ca w funkcji czasu - skok CAin')
    hold on
end
legend show
hold off

for iter = 1:1:length(CAin_vect)
    CAin = CAin_vect(iter);
    [y, t] = rk4(@dCa, @dT, Ca, T, step);
    txt = ['CAin = ',num2str(CAin_vect(iter))];
    figure(2)
    plot(t,y(2,:),'DisplayName',txt);
    title('Temperatura T w funkcji czasu - skok CAin')
    hold on
end
legend show
hold off
CAin = 2;
for iter = 1:1:length(Fc_vect)
    Fc = Fc_vect(iter);
    [y, t] = rk4(@dCa, @dT, Ca, T, step);
    txt = ['Fc = ',num2str(Fc_vect(iter))];
    figure(3)
    plot(t,y(1,:),'DisplayName',txt);
    title('Stê¿enie Ca w funkcji czasu - skok Fc')
    hold on
end
legend show
hold off



for iter = 1:1:length(Fc_vect)
    Fc = Fc_vect(iter);
    [y, t] = rk4(@dCa, @dT, Ca, T, step);
    txt = ['Fc = ',num2str(Fc_vect(iter))];
    figure(4)
    plot(t,y(2,:),'DisplayName',txt);
    title('Temperatura T w funkcji czasu - skok Fc')
    hold on
end
legend show
hold off
Fc = 15;
% =============================
for iter = 1:1:length(CAin_vect)
    CAin = CAin_vect(iter);
    [y, t] = rk4(@dCaLin, @dTLin, Ca, T, step);
    txt = ['CAin = ',num2str(CAin_vect(iter))];
    figure(5)
    plot(t,y(1,:),'DisplayName',txt);
    title('Stê¿enie Ca w funkcji czasu - skok CAin')
    hold on
end
legend show
hold off

for iter = 1:1:length(CAin_vect)
    CAin = CAin_vect(iter);
    [y, t] = rk4(@dCaLin, @dTLin, Ca, T, step);
    txt = ['CAin = ',num2str(CAin_vect(iter))];
    figure(6)
    plot(t,y(2,:),'DisplayName',txt);
    title('Temperatura T w funkcji czasu - skok CAin')
    hold on
end
legend show
hold off
CAin = 2;
for iter = 1:1:length(Fc_vect)
    Fc = Fc_vect(iter);
    [y, t] = rk4(@dCaLin, @dTLin, Ca, T, step);
    txt = ['Fc = ',num2str(Fc_vect(iter))];
    figure(7)
    plot(t,y(1,:),'DisplayName',txt);
    title('Stê¿enie Ca w funkcji czasu - skok Fc')
    hold on
end
legend show
hold off

for iter = 1:1:length(Fc_vect)
    Fc = Fc_vect(iter);
    [y, t] = rk4(@dCaLin, @dTLin, Ca, T, step);
    txt = ['Fc = ',num2str(Fc_vect(iter))];
    figure(8)
    plot(t,y(2,:),'DisplayName',txt);
    title('Temperatura T w funkcji czasu - skok Fc')
    hold on
end
legend show
hold off
Fc = 15;
% =============================
% 1 step - 0.01 min
  Ts = 10; %  -> sampling time = 0.1min = 6s

for iter = 1:1:length(CAin_vect)
    CAin = CAin_vect(iter);
    [y, t] = rk4Discrete(@dCaLin, @dTLin, Ca, T, step,Ts);
    txt = ['CAin = ',num2str(CAin_vect(iter)/1000)];
    figure(9)
    plot(t,y(1,:),'DisplayName',txt);
    title('Stê¿enie Ca w funkcji czasu - skok CAin')
    hold on
end
legend show
hold off

for iter = 1:1:length(CAin_vect)
    CAin = CAin_vect(iter);
    [y, t] = rk4Discrete(@dCaLin, @dTLin, Ca, T, step,Ts);
    txt = ['CAin = ',num2str(CAin_vect(iter)/1000)];
    figure(10)
    plot(t,y(2,:),'DisplayName',txt);
    title('Temperatura T w funkcji czasu - skok CAin')
    hold on
end
legend show
hold off
CAin = 2;
for iter = 1:1:length(Fc_vect)
    Fc = Fc_vect(iter);
    [y, t] = rk4Discrete(@dCaLin, @dTLin, Ca, T, step,Ts);
    txt = ['Fc = ',num2str(Fc_vect(iter))];
    figure(11)
    plot(t,y(1,:),'DisplayName',txt);
    title('Stê¿enie Ca w funkcji czasu - skok Fc')
    hold on
end
legend show
hold off

for iter = 1:1:length(Fc_vect)
    Fc = Fc_vect(iter);
    [y, t] = rk4Discrete(@dCaLin, @dTLin, Ca, T, step,Ts);
    txt = ['Fc = ',num2str(Fc_vect(iter))];
    figure(12)
    plot(t,y(2,:),'DisplayName',txt);
    title('Temperatura T w funkcji czasu - skok Fc')
    hold on
end
legend show
hold off
Fc = 15;
function [y,t] = rk4(dCa,dT,Ca,T,step)
t=0:step:5;
y(:,1) = [Ca T];

    for i=1:(length(t)-1)
        k11=dCa(Ca,T);
        k12=dT(Ca,T);

        k21=dCa(Ca+0.5*step,T+0.5*step*k11);
        k22=dT(Ca+0.5*step,T+0.5*step*k12);

        k31=dCa(Ca+0.5*step,T+0.5*step*k21);
        k32=dT(Ca+0.5*step,T+0.5*step*k22);

        k41=dCa(Ca+step,T+step*k31);
        k42=dT(Ca+step,T+step*k32);

        Ca=Ca+(step/6)*(k11+k41+2*(k21+k31));
        T=T+(step/6)*(k12+k42+2*(k22+k32));
        y(:,i+1)=[Ca T];
    end
end

function [y,t] = rk4Discrete(dCa,dT,Ca,T,step,Ts)
t=0:step:5;
y(:,1) = [Ca T];
acc = [0 0];
%wyliczanie wspo³czynników
for i=1:(length(t)-1)
    k11=dCa(Ca,T);
    k12=dT(Ca,T);

    k21=dCa(Ca+0.5*step,T+0.5*step*k11);
    k22=dT(Ca+0.5*step,T+0.5*step*k12);

    k31=dCa(Ca+0.5*step,T+0.5*step*k21);
    k32=dT(Ca+0.5*step,T+0.5*step*k22);

    k41=dCa(Ca+step,T+step*k31);
    k42=dT(Ca+step,T+step*k32);

    Ca=Ca+(step/6)*(k11+k41+2*(k21+k31));
    T=T+(step/6)*(k12+k42+2*(k22+k32));
    
    if(mod(i-1,Ts) == 0)
        acc = [Ca T];
    end
    y(:,i+1)=acc;
end
end

function [dCa] = dCa(Ca, T)
global k E_R V Fin CAin;
dCa = (Fin*CAin - Fin*Ca - V*k*exp(-E_R/T)*Ca)/V;
end

function [dT] = dT(Ca,T)
global  ro cp k E_R h a b V Fin Fc Tin ;
dT = (Fin*ro*cp*Tin - Fin*ro*cp*T + V*h*k*exp(-E_R/T)*Ca - (a*Fc^(b+1)/(Fc+(a*Fc^b/(2*ro*cp))))*(T-Tin))/(V*ro*cp);
end

function [dCa] = dCaLin(Ca, T)
global k E_R V Fin CAin Ca0 T0 F CAin0;
dCa = -((Ca - Ca0)*(F + V*k*exp(-E_R/T0)) + Ca0*F - CAin0*Fin - Fin*(CAin - CAin0) + Ca0*V*k*exp(-E_R/T0) + (Ca0*E_R*V*k*exp(-E_R/T0)*(T - T0))/T0^2)/V;
end

function [dT] = dTLin(Ca,T)
global  ro cp k E_R h a b V Fin Fc Fc0 Tin Ca0 Tin T0 F ;
dT = (a*(Fc - Fc0)*(T0 - Tin)*((exp(log(Fc0)*(b + 1))...
    *((a*b*exp(b*log(Fc0)))/(2*Fc0*cp*ro) + 1))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro))^2 ...
    - (exp(log(Fc0)*(b + 1))*(b + 1))/(Fc0*(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro))))...
    - (T - T0)*((a*exp(log(Fc0)*(b + 1)))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro))...
    + F*cp*ro - (Ca0*E_R*V*h*k*exp(-E_R/T0))/T0^2)...
    - (a*exp(log(Fc0)*(b + 1))*(T0 - Tin))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro))...
    - F*T0*cp*ro + Fin*Tin*cp*ro + Ca0*V*h*k*exp(-E_R/T0) + V*h*k*exp(-E_R/T0)*(Ca - Ca0))/(V*cp*ro);
 end
