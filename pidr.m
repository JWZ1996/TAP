models
clearvars -except model G m1 Gz1 CAin0 Ca0 Fc0 T0 
modelDiscrete = m1;
modelCont = model;
Gz = Gz1;
Gs = G;
clearvars model G m1 Gz1 
CainIN = 1;
FcIN = 2;
CaOUT = 1;
TOUT = 2;

Gs.InputName = {'Cain','Fc'};
Gs.OutputName = {'Ca','T'};

%%
% =========================================================================
% Macierz RGA

[~,k11] = zero(modelCont(1,1));
[~,k12] = zero(modelCont(1,2));
[~,k21] = zero(modelCont(2,1));
[~,k22] = zero(modelCont(2,2));

K = [k11 k21; k12 k22];   

RGA = K.*transpose(inv(K));

% wkradl sie wczesnej blad:
% RGA =
%     0.0071    0.9929
%     0.9929    0.0071
%

% 1.(Cain) Wejscie reguluje 1. (Ca) wyjscie
% 2.(Fc)   Wejscie reguluje 2. (T) wyjscie 

%%
% =========================================================================
% dostrajanie regulatora PI w klasycznej odmianie
C0 = pidstd(1,1);
% Controller1 i Controller2 - dostrojone regulatory PI przy otwartej
% drugiej petli ( wowczas jest tylko jedna petla dla kazdego regulatora i w
% niej transmitancja jest - odpowiednio - z wejscia 1 do wyjsca 2
% (CaIN-T)iz wejscia 2 do wyjscia 1 (Fc-Ca)
C1 = pidtune(Gs(1,2), C0);
C2 = pidtune(Gs(2,1), C0);

R = [C1 0; 0 C2];   % transmitancje regulatorow

%%
% =========================================================================
% Regulatory rozłączone
RGs = R*Gs;

lineCainT = feedback(RGs(CainIN,TOUT),1);
lineFcCa    = feedback(RGs(FcIN,CaOUT),1);

%%
% =========================================================================
% Regulatory złączone, wyłączone odsprzęganie

% transmitancje z wejsc do wyjsc - perfidnie zapier... z wykladu 
% https://studia2.elka.pw.edu.pl/file/21L/103C-ARxxx-MSP-TAP/priv/01RegWielopetlowa%5FSlajdy.pdf
% slajd 12
G12 = R(1,1)*(Gs(1,2)-(Gs(1,2)*Gs(2,1)*R(2,2))/(1+R(2,2)*Gs(2,1)));
G21 = R(2,2)*(Gs(2,1)-(Gs(2,1)*Gs(1,2)*R(1,1)/(1+R(1,1)*Gs(1,2))));


timeVector = 0:0.05:10;

% o ile zmienia sie wartosc zadana w stosunku do punktu pracy
T_step_size = -5;
Ca_step_size=2.5;

figure;
pidTest(lineCainT,Ca0,Ca_step_size,timeVector,true);
hold on;
pidTest(G12,Ca0,Ca_step_size,timeVector,true);
ylabel("Temperatura [K]");
xlabel("Czas [min]");
legend("1-2/Cain-T - Otwarta petla","Zadana wartość", "1-2/Cain-T - Zamknieta petla")
hold off;

figure;
pidTest(lineFcCa,T0,T_step_size,timeVector,true);
hold on;
pidTest(G21,T0,T_step_size,timeVector,true);
ylabel("Stężenie [kmol/m^3]");
xlabel("Czas [min]");
legend("2-1/Fc-Ca - Otwarta petla ","Zadana wartość", "2-1/Fc-Ca - Zamknieta petla")
hold off;


%%
% % =========================================================================
% % Regulatory złączone, włączone odsprzęganie 
 D11=-Gs(1,1)/Gs(2,1);
 D22=-Gs(2,2)/Gs(1,2);
% 
% % nie dziala jeszcze!
% D = [0 D21; D12 0];
G12c = R(1,1)*(Gs(1,2)-(Gs(1,2)*Gs(2,1)*R(2,2) + D22*Gs(2,2) )/(1+R(2,2)*Gs(2,1) + D22*Gs(2,2)) );
G21c = R(2,2)*(Gs(2,1)-(Gs(2,1)*Gs(1,2)*R(1,1) + D11*Gs(1,1) )/(1+R(1,1)*Gs(1,2) + D11*Gs(1,1)) );


figure;
hold on;
pidTest(G21,T0,T_step_size,timeVector,true);
pidTest(G21c,T0,T_step_size,timeVector,true);
legend("2-1 bez odsprzęgania"," Odsprzeganie 2-1");
ylabel("Temperatura [K]");
xlabel("Czas [min]");
hold off;

figure;
hold on;
pidTest(G21,T0,T_step_size,timeVector,true);
pidTest(G12c,Ca0,Ca_step_size,timeVector,true);
legend("1-2 bez odsprzęgania"," Odsprzeganie 1-2");
ylabel("Stężenie [kmol/m^3]");
xlabel("Czas [min]");
hold off;

x = BLTwsklog(1.48e-06,6.85e-05,0.00377,6.85e-05,Gs)

function Lcmax = BLTwsklog(kp1,Ti1,kp2,Ti2,G)
    %podanie modelu obiektu 2x2:
    s=tf('s');
    G11=G(1,1);
    G12=G(1,2);
    G21=G(2,1);
    G22=G(2,2);

    %regulatory PI:
    R1=kp1*(1+(1/(Ti1*s)));
    R2=kp2*(1+(1/(Ti2*s)));

    G11R1=series(G11,R1); %pierwszy składnik sumy
    G22R2=series(G22,R2); %drugi składnik sumy
    G11R1ss=ss(G11R1);
    G22R2ss=ss(G22R2);
    W1=parallel(G11R1ss,G22R2ss); %pierwsze sumowanie
    W1ss=ss(W1);
    GGRR1=series(G11R1ss,G22R2ss); %trzeci składnik sumy
    W2=parallel(GGRR1,W1ss); %drugie sumowanie
    G12G21=series(G12,G21);
    R1R2=series(R1,R2);
    GGRR2=series(G12G21,R1R2);
    GGRR2m=series(GGRR2,-1); %czwarty składnik sumy (suma z plusem)
    GGRR2mss=ss(GGRR2m);
    W=parallel(W2,GGRR2mss); %trzecie sumowanie (ostatnie)
    % ch-ka częstotliwościowa:
    wektw=logspace(-4,0,100); % 100 punktów od 10^-4 do 10^0
    Wfr=freqresp(W,wektw); %wektor ch-ki częstotliwościowej W
    Wabs=abs(Wfr);
    W1=parallel(W,1); %1+W
    Wfr1=freqresp(W1,wektw); %wektor ch-ki częstotliwościowej 1+W
    W1abs=abs(Wfr1);
    Lclogabs=20*(log10(Wabs)-log10(W1abs)); %wektor wartości funkcji wskaźnika
    
    figure;
    semilogx(wektw*100/60 ,Lclogabs(:,:));
    ylabel("Lc [dB]");
    xlabel("Częstotliwość [Hz]");
    Lcmax=max(Lclogabs); %wartość kryterium BLT
end

function [y,t] = pidTest(Gs, outputEquilVal, requestedValue, timeVector,isPlotted)
    opt = stepDataOptions('InputOffset',0,'StepAmplitude',requestedValue);
    y = step(Gs,timeVector,opt);

    y=y+outputEquilVal;

    if isPlotted
        plot(timeVector,y);
        yline(outputEquilVal+requestedValue);
    end
end

