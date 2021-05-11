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

% =========================================================================
% Macierz RGA

[~,k11] = zero(modelCont(1,1));
[~,k12] = zero(modelCont(1,2));
[~,k21] = zero(modelCont(2,1));
[~,k22] = zero(modelCont(2,2));

K = [k11 k21; k12 k22];

RGA = K.*transpose(inv(K));

%RGA =
%    0.9929    0.0071
%    0.0071    0.9929

%
% 1.(Cain) Wejscie reguluje 1. (Ca) wyjscie
% 2.(Fc)   Wejscie reguluje 2. (T) wyjscie 

% =========================================================================
% Nastawy regulatora PI
G11 = Gs(CainIN,CaOUT) + Gs(FcIN,CaOUT);
G22 = Gs(FcIN,TOUT)    + Gs(CainIN,TOUT);

% Nastawy dobrane za pomocą aplikacji pidtuner
feedbackLineStub = tf(1);
feedbackLineStubVec = [tf(1) tf(1)];
timeVector = 0:0.05:5;

kp1 = 0.01542;
Ti1 = 0.3186;

kp2 = -2.802e-5;
Ti2 = 6.852e-5;

R = tf([0 0; 0 pidstd(kp2,Ti2)]);


% =========================================================================
% Regulatory rozłączone
RGs = R*Gs;

lineCainCa = feedback(RGs(CainIN,CaOUT),feedbackLineStub);
lineFcT    = feedback(RGs(FcIN,TOUT),feedbackLineStub);


figure;
pidTest(lineFcT,Ca0,1,timeVector,true);

figure;
pidTest(lineCainCa,T0,1,timeVector,true);
% =========================================================================
% Regulatory złączone, wyłączone odsprzęganie
nonCoupled = feedback(feedback(RGs,feedbackLineStub,FcIN,TOUT,-1),feedbackLineStub,CainIN,CaOUT,-1);


timeVector = 0:0.05:10;

figure;
pidTest(nonCoupled(CainIN,CaOUT) + nonCoupled(FcIN,CaOUT),Ca0,1,timeVector,true);



figure;
pidTest(nonCoupled(FcIN,TOUT) + nonCoupled(CainIN,TOUT),T0,1,timeVector,true);




% =========================================================================
% Regulatory złączone, włączone odsprzęganie
% D21=-Gs(2,1)/Gs(2,2);
% D12=-Gs(1,2)/Gs(1,1);
% 
% D = [0 D21; D12 0];
% 
% coupled = feedback(feedback(R*D*Gs,feedbackLineStub,CainIN,CaOUT,-1),feedbackLineStub,FcIN,TOUT,-1);
% 
% figure;
% pidTest(coupled(CainIN,CaOUT) + coupled(FcIN,CaOUT),Ca0,1,timeVector,true);
% 
% figure;
% pidTest(coupled(FcIN,TOUT)+ coupled(CainIN,CaOUT),T0,1,timeVector,true);



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
    Lcmax=max(Lclogabs); %wartość kryterium BLT
end

function [y,t] = pidTest(Gs, outputEquilVal, requestedValue, timeVector,isPlotted)
    opt = stepDataOptions('InputOffset',0,'StepAmplitude',requestedValue-outputEquilVal);
    y = step(Gs,timeVector,opt);

    y=y+outputEquilVal;

    if isPlotted
        plot(timeVector,y);
        yline(requestedValue);
    end
end


