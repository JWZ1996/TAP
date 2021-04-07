close all;
clear all;

% ==== Plik gÅ‚owny =====
global g m K mol kmol cal min ro cp k E_R h a b ro cp k E_R h a b V Fin CAin Fc Tin Tcin Ca T;
consts
% === Symulacja obiektu - RK4 ===
step = 0.01;
Ca = 0.042953 * kmol/m^3;
T =  431.176 * K;

CAin_vect = [2.6, 2.4, 2.2, 1.8, 1.6, 1.4, 2].*kmol;
Fc_vect = [12,13,14,16,17,18, 15];

%% ============================
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
plot(y(1,:),y(2,:), '*')
%% ============================
%==========    b)   =============

figure(1)
PlotModel('CAin', 'Ca', CAin_vect, @dCa, @dT, step, 'St�enie Ca w funkcji czasu - skok CAin', ':');
PlotModel('CAin', 'Ca', CAin_vect, @dCaLin, @dTLin, step, 'St�enie Ca w funkcji czasu - skok CAin', '--');
plotLegend('Ca', CAin_vect)
hold off

CAin = 2; Fc = 15;

figure(2)
PlotModel('Fc', 'Ca', Fc_vect, @dCa, @dT, step, 'St�enie Ca w funkcji czasu - skok Fc', ':');
PlotModel('Fc', 'Ca', Fc_vect, @dCaLin, @dTLin, step, 'St�enie Ca w funkcji czasu - skok Fc', '--');
plotLegend('Fc', Fc_vect)
hold off

CAin = 2; Fc = 15;

figure(3)
PlotModel('CAin', 'T', CAin_vect, @dCa, @dT, step, 'Temperatura T w funkcji czasu - skok CAin', ':');
PlotModel('CAin', 'T', CAin_vect, @dCaLin, @dTLin, step, 'Temperatura T w funkcji czasu - skok CAin', '--');
plotLegend('Ca', CAin_vect)
hold off

CAin = 2; Fc = 15;

figure(4)
PlotModel('Fc', 'T', Fc_vect, @dCa, @dT, step, 'Temperatura T w funkcji czasu - skok Fc', ':');
PlotModel('Fc', 'T', Fc_vect, @dCaLin, @dTLin, step, 'Temperatura T w funkcji czasu - skok Fc', '--');
plotLegend('Fc', Fc_vect)
hold off

CAin = 2; Fc = 15;
%% =============================
%========Model nieliniowy=========
figure(1)
PlotModel('CAin', 'Ca', CAin_vect, @dCa, @dT, step, 'St�enie Ca w funkcji czasu - skok CAin', '-');
hold off

figure(2)
PlotModel('CAin', 'T', CAin_vect, @dCa, @dT, step, 'Temperatura T w funkcji czasu - skok CAin', '--');
hold off

CAin = 2;

figure(3)
PlotModel('Fc', 'Ca', Fc_vect, @dCa, @dT, step, 'St�enie Ca w funkcji czasu - skok Fc', '-');
hold off

figure(4)
PlotModel('Fc', 'T', Fc_vect, @dCa, @dT, step, 'Temperatura T w funkcji czasu - skok Fc', '--');
hold off

Fc = 15;
%% =============================
%======Model zlinearyzowany=======
figure(5)
PlotModel('CAin', 'Ca', CAin_vect, @dCaLin, @dTLin, step, 'St�enie Ca w funkcji czasu - skok CAin', '-');
hold off

figure(6)
PlotModel('CAin', 'T', CAin_vect, @dCaLin, @dTLin, step, 'Temperatura T w funkcji czasu - skok CAin', '-');
hold off

CAin = 2;

figure(7)
PlotModel('Fc', 'Ca', Fc_vect, @dCaLin, @dTLin, step, 'St�enie Ca w funkcji czasu - skok Fc', '-');
hold off

figure(8)
PlotModel('Fc', 'T', Fc_vect, @dCaLin, @dTLin, step, 'Temperatura T w funkcji czasu - skok Fc', '-');
hold off

Fc = 15;

%% =============================
%=========Model dyskretny=========
% 1 step - 0.01 min
Ts = 10; %  -> sampling time = 0.1min = 6s

figure(9)
PlotModelDiscrete('CAin', 'Ca', Ts, CAin_vect, @dCaLin, @dTLin, step, 'St�enie Ca w funkcji czasu - skok CAin', '-');
hold off

figure(10)
PlotModelDiscrete('CAin', 'T', Ts, CAin_vect, @dCaLin, @dTLin, step, 'Temperatura T w funkcji czasu - skok CAin', '-'); 
hold off

CAin = 2;
figure(11)
PlotModelDiscrete('Fc', 'Ca', Ts, Fc_vect, @dCaLin, @dTLin, step, 'St�enie Ca w funkcji czasu - skok Fc', '-');
hold off

figure(12)
PlotModelDiscrete('Fc', 'T', Ts, Fc_vect, @dCaLin, @dTLin, step, 'Temperatura T w funkcji czasu - skok Fc', '-'); 
hold off

Fc = 15;
%%=================================================
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
    colour_vect = [ 0.9294,    0.6941,    0.1255;
                    0,    	   0.4471,    0.7412;
                    0.4941,    0.1843,    0.5569;
                    0.1490,    0.1490,    0.1490;
                    0.6353,    0.0784,    0.1843;
                    1.0000,    0.4118,    0.1608;
                    0.0549,    0.6196,    0.9020];
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
        [y, t] = rk4(dCaIn, dTIn, Ca, T, step);
        %txt = [InVar ' = ',num2str(vect(iter))];
        plot(t,y(y_plot,:),plotStyle,'Color', colour_vect(iter,:),'LineWidth',1.5); %'DisplayName',txt,
        title(titleText)
        xlabel('t [min]')
        if (strcmp(OutVar,'Ca') == 1)
            ylabel('Ca [kmol/m^3]');
        elseif (strcmp(OutVar,'T') == 1)
            ylabel('T [K]');  
        end
        hold on
    end
    %legend show
end

function PlotModelDiscrete(InVar, OutVar, Ts, vect, dCaIn, dTIn, step, titleText, plotStyle)
    global Ca T CAin Fc;
    colour_vect = [ 0.9294,    0.6941,    0.1255;
                    0,    	   0.4471,    0.7412;
                    0.4941,    0.1843,    0.5569;
                    0.1490,    0.1490,    0.1490;
                    0.6353,    0.0784,    0.1843;
                    1.0000,    0.4118,    0.1608;
                    0.0549,    0.6196,    0.9020];
    for iter = 1:1:length(vect)
        if (strcmp(InVar,'CAin') == 1)
            CAin = vect(iter);
        elseif (strcmp(InVar,'Fc') == 1)
            Fc = vect(iter);   
        end
        if (strcmp(OutVar,'Ca') == 1)
            y_plot = 1;
            ylabel('Ca [kmol/m^3]');
        elseif (strcmp(OutVar,'T') == 1)
            y_plot = 2;
            ylabel('T [K]'); 
        end
   
        [y, t] = rk4Discrete(dCaIn, dTIn, Ca, T, step, Ts);
        %txt = [InVar ' = ',num2str(vect(iter))];
        plot(t,y(y_plot,:),plotStyle,'Color',colour_vect(iter,:),'LineWidth',1.5); %'DisplayName',txt,
        title(titleText)
        if (strcmp(OutVar,'CA') == 1)
            ylabel('Ca [kmol/m^3]');
        elseif (strcmp(OutVar,'T') == 1)
            ylabel('T [K]');  
        end
        hold on
    end
    %legend show
end

function plotLegend(InVar, vect)
colour_vect = [ 0.9294,    0.6941,    0.1255;
                0,    	   0.4471,    0.7412;
                0.4941,    0.1843,    0.5569;
                0.1490,    0.1490,    0.1490;
                0.6353,    0.0784,    0.1843;
                1.0000,    0.4118,    0.1608;
                0.0549,    0.6196,    0.9020];
h = zeros(size(colour_vect, 1));
h(1) = plot(NaN,NaN,'-','Color',colour_vect(1,:));
h(2) = plot(NaN,NaN,'-','Color',colour_vect(2,:));
h(3) = plot(NaN,NaN,'-','Color',colour_vect(3,:));
h(4) = plot(NaN,NaN,'-','Color',colour_vect(4,:));
h(5) = plot(NaN,NaN,'-','Color',colour_vect(5,:));
h(6) = plot(NaN,NaN,'-','Color',colour_vect(6,:));
h(7) = plot(NaN,NaN,'-','Color',colour_vect(7,:));
equals = ' = ';
legend(   char(strcat(InVar,equals,{' '},num2str(vect(1)))), char(strcat(InVar,equals,{' '},num2str(vect(2)))), ...
          char(strcat(InVar,equals,{' '},num2str(vect(3)))), char(strcat(InVar,equals,{' '},num2str(vect(4)))), ...
          char(strcat(InVar,equals,{' '},num2str(vect(5)))), char(strcat(InVar,equals,{' '},num2str(vect(6)))), ...
          char(strcat(InVar,equals,{' '},num2str(vect(7)))) ...
          );
end
