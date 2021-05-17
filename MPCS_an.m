clear all;
close all;

global Ca0 T0 CAin0 Fc0;

[A, B, C, D] = import_ABCD();
Ts = 0.01; Tsim = 30;
contin=ss(A,B,C,D);
disc=c2d(contin,Ts);
A=disc.A;
B=disc.B;
C=disc.C;
D=disc.D;

% horyzont predykcji; horyzont sterowanie; rozmiar wyjść; rozmiar sterowań
N=10;Nu=3;
ny=2;nu=2;  

%wyznaczanie macierzy algorytmu MPCS
[M,CtAt,CtV]=MPCSmatrices(A,B,C,N,Nu);

%macierze kar
lambda1 = 0.005;
lambda2 = 0.005;
psi1 = 30;
psi2 = 30;

psi=[psi1 psi2];
psi=diag(repmat(psi,1,N));

lambda=[lambda1 lambda2];
lambda=diag(repmat(lambda,1,Nu));

Kmat=(M'*psi*M+lambda)^-1*M'*psi; 
K1=Kmat(1:nu,:); %nu pierwszych wierszy

time=0:Ts:Tsim-Ts;

%inicjalizacja pierwszej iteracji
iter = 2:Tsim/Ts;

U=zeros(nu,length(iter)+1); 
v=zeros(nu,length(iter)+1); 
y=zeros(ny,length(iter)+1);
x=zeros(ny,length(iter)+1);
dU=zeros(ny,length(iter)+1);
Yzad=zeros(N*ny,length(iter)+1);

for iter = 2:Tsim/Ts-1
    %trajektoria zadana
    if(iter<0.1*length(time))
        yzad=[0; 0]; 
        Yzad(1:2:end,iter)=yzad(1);
        Yzad(2:2:end,iter)=yzad(2);
    elseif(iter<0.2*length(time)) 
        yzad=[0; -2]; 
        Yzad(1:2:end,iter)=yzad(1);
        Yzad(2:2:end,iter)=yzad(2);    
    elseif(iter<0.4*length(time)) 
        yzad=[0.01; -2]; 
        Yzad(1:2:end,iter)=yzad(1);
        Yzad(2:2:end,iter)=yzad(2);
    elseif(iter<0.6*length(time)) 
        yzad=[0.01; 2]; 
        Yzad(1:2:end,iter)=yzad(1);
        Yzad(2:2:end,iter)=yzad(2);
	elseif(iter<0.8*length(time)) 
        yzad=[-0.01; -5];
        Yzad(1:2:end,iter)=yzad(1);
        Yzad(2:2:end,iter)=yzad(2);
    else
        yzad=[0; 0]; 
        Yzad(1:2:end,iter)=yzad(1);
        Yzad(2:2:end,iter)=yzad(2);
    end
    
    v(:,iter)=x(:,iter)-A*x(:,iter-1)-B*U(:,iter-1); %W5 10/47, sam dół  
    dU(:,iter)=K1*(Yzad(:,iter)-CtAt*x(:,iter)-CtV*B*U(:,iter-1)-CtV*v(:,iter)); %W5 10/47
    U(:,iter)=U(:,iter-1)+dU(:,iter); %wyrażenie przyrostowego charakteru regulacji
    
    U(1,iter)=min(U(1,iter), 3 - CAin0); %ograniczenie na CAin
    U(1,iter)=max(U(1,iter), 1 - CAin0); %odjęcie wartości dla punktu pracy: przyrost -1 oznacza, że wartość spadnie z 2 do 1
    
    U(2,iter)=min(U(2,iter),30-Fc0); %ograniczenia na Fc
    U(2,iter)=max(U(2,iter),0-Fc0); %ograniczenie do wyschnięcia rur, być może nierozsądne
    
    y(:,iter)=C*x(:,iter)+D*U(:,iter); %równania stanu zlinearyzowane
    x(:,iter+1)=A*x(:,iter)+B*U(:,iter); 
end

subplot(2,2,1);
grid on;
hold on;
plot(time,y(1,:)+Ca0);
hold on
plot(time,Yzad(1,:)+Ca0, 'm--')

title('Wyjście 1 - stężenie substancji A')
xlabel('Czas symulacji [min]')
ylabel('Ca [kmol/m3]')
legend('MPCS', 'PID', 'Wartość zadana', 'bool')

hold off
subplot(2,2,2);
grid on;
hold on;
plot(time,y(2,:)+T0);
hold on;
plot(time,Yzad(2,:)+T0, 'm--')
title('Wyjście 2 - temperatura')
xlabel('Czas symulacji [min]')
ylabel('T [K]')
legend('MPCS', 'PID', 'Wartość zadana')

subplot(2,2,3);
grid on;
hold on;
plot(time,U(1,:)+CAin0);
title('Sterowanie 1 - wejściowe stężenie substancji A')
xlabel('Czas symulacji [min]')
ylabel('CAin [kmol/m3]') 

subplot(2,2,4);
grid on;
hold on;
plot(time,U(2,:)+Fc0);
title('Sterowanie 2 - strumień objętości substancji chłodzącej')
xlabel('Czas symulacji [min]')
ylabel('Fc [m3/min]')

