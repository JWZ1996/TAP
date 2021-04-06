close all;
clear all;

% ==== Plik gÅ‚owny =====
global g m K mol kmol cal min ro cp k E_R h a b ro cp k E_R h a b V Fin CAin Fc Tin Tcin Ca T;
consts
% === Symulacja obiektu - RK4 ===
time_step = 0.01;
Ca = 0.042953 * kmol/m^3;
T =  431.176 * K;

CAin_vect = [2.4, 2.2, 1.8, 1.6, 1.4 2].*kmol;
Fc_vect = [18,17,16,14,13,12, 15];

%% ============================
%==========    b)   =============
figure(1)
PlotModel('CAin', 'Ca', CAin_vect, @dCa, @dT, time_step, 'Stê¿enie Ca w funkcji czasu - skok CAin', ':');
PlotModel('CAin', 'Ca', CAin_vect, @dCaLin, @dTLin, time_step, 'Stê¿enie Ca w funkcji czasu - skok CAin', '--');
hold off

CAin = 2;

figure(2)
PlotModel('Fc', 'Ca', Fc_vect, @dCa, @dT, time_step, 'Stê¿enie Ca w funkcji czasu - skok Fc', ':');
PlotModel('Fc', 'Ca', Fc_vect, @dCaLin, @dTLin, time_step, 'Stê¿enie Ca w funkcji czasu - skok Fc', '--');
hold off

figure(3)
PlotModel('CAin', 'T', CAin_vect, @dCa, @dT, time_step, 'Temperatura T w funkcji czasu - skok CAin', ':');
PlotModel('CAin', 'T', CAin_vect, @dCaLin, @dTLin, time_step, 'Temperatura T w funkcji czasu - skok CAin', '--');
hold off

figure(4)
PlotModel('Fc', 'T', Fc_vect, @dCa, @dT, time_step, 'Temperatura T w funkcji czasu - skok Fc', ':');
PlotModel('Fc', 'T', Fc_vect, @dCaLin, @dTLin, time_step, 'Temperatura T w funkcji czasu - skok Fc', '--');
hold off

Fc = 15;
%% =============================
%========Model nieliniowy=========
figure(1)
PlotModel('CAin', 'Ca', CAin_vect, @dCa, @dT, time_step, 'Stê¿enie Ca w funkcji czasu - skok CAin', '-');
hold off

figure(2)
PlotModel('CAin', 'T', CAin_vect, @dCa, @dT, time_step, 'Temperatura T w funkcji czasu - skok CAin', '--');
hold off

CAin = 2;

figure(3)
PlotModel('Fc', 'Ca', Fc_vect, @dCa, @dT, time_step, 'Stê¿enie Ca w funkcji czasu - skok Fc', '-');
hold off

figure(4)
PlotModel('Fc', 'T', Fc_vect, @dCa, @dT, time_step, 'Temperatura T w funkcji czasu - skok Fc', '--');
hold off

Fc = 15;
%% =============================
%======Model zlinearyzowany=======
figure(5)
PlotModel('CAin', 'Ca', CAin_vect, @dCaLin, @dTLin, time_step, 'Stê¿enie Ca w funkcji czasu - skok CAin', '-');
hold off

figure(6)
PlotModel('CAin', 'T', CAin_vect, @dCaLin, @dTLin, time_step, 'Temperatura T w funkcji czasu - skok CAin', '-');
hold off

CAin = 2;

figure(7)
PlotModel('Fc', 'Ca', Fc_vect, @dCaLin, @dTLin, time_step, 'Stê¿enie Ca w funkcji czasu - skok Fc', '-');
hold off

figure(8)
PlotModel('Fc', 'T', Fc_vect, @dCaLin, @dTLin, time_step, 'Temperatura T w funkcji czasu - skok Fc', '-');
hold off

Fc = 15;

%% =============================
%=========Model dyskretny=========
% wektor roznych okresow probkowania[w sekundach - wszedzie podstawowa jednostka jest minuta
% zatem przemnozony wektor zostal przez 60 s]

Ts = [0.0001 0.001 0.01 0.1]*60*min;
for Ts_step = Ts
    [cont_ss, cont_tf, discr_ss, discr_tf] = generate_models (Ts_step);
    
%   Generowanie wykresow nalozonych odpowiedzi skokowych dla transmitancji
%   ciaglej oraz zdyskretyzowanej. 
    figure;
    step(cont_tf, discr_tf);
    legend('Ciagla','Dyskretna','Location','SouthEast')
    str = sprintf('Odpowiedzi skokowe modelow Ts: %0.4f', Ts_step);
    title(str);
end


%% =================================================
%=====================FUNKCJE======================
%==================================================

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
 
function PlotModel(InVar, OutVar, vect, dCaIn, dTIn, step, titleText, plotStyle)
    global Ca T CAin Fc;
    for iter = 1:1:length(vect)
        if (strcmp(InVar,'CAin') == 1)
            CAin = vect(iter);
        elseif (strcmp(InVar,'Fc') == 1)
            Fc = vect(iter);   
        end
        if (strcmp(OutVar,'Ca') == 1)
            y_plot = 1;
        elseif (strcmp(OutVar,'T') == 1)
            y_plot = 2;   
        end
           
        [y, t] = rk4(dCaIn, dTIn, Ca, T, step); %cze�� Zan Jembowicz tutaj mam pytanie we� no sp�jrz kurcz� jak tu wstawi� do @dCa wszystko na co mam ochot�
        txt = [InVar ' = ',num2str(vect(iter))];
        plot(t,y(y_plot,:),plotStyle,'DisplayName',txt);
        title(titleText)
        hold on
    end
    legend show
end %o tu pytanie b�dzie Jan

function PlotModelDiscrete(InVar, OutVar, Ts, vect, dCa, dT, step, titleText, plotStyle)
    global Ca T CAin Fc;
    for iter = 1:1:length(vect)
        if (strcmp(InVar,'CAin') == 1)
            CAin = vect(iter);
        elseif (strcmp(InVar,'Fc') == 1)
            Fc = vect(iter);   
        end
        if (strcmp(OutVar,'Ca') == 1)
            y_plot = 1;
        elseif (strcmp(OutVar,'T') == 1)
            y_plot = 2;   
        end
           
        [y, t] = rk4Discrete(@dCa, @dT, Ca, T, step, Ts);
        txt = [InVar ' = ',num2str(vect(iter))];
        plot(t,y(y_plot,:), plotStyle,'DisplayName',txt);
        title(titleText)
        hold on
    end
    legend show
end

