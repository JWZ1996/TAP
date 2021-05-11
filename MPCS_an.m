close all;
clear all;

global CAin Fc Ca T Ca0 T0 Cain Fc0;

[A, B, C, D] = import_ABCD();
Ts = 0.001; Tsim = 300;
contin=ss(A,B,C,D);
disc=c2d(contin,Ts);
A=disc.A;
B=disc.B;
C=disc.C;
D=disc.D;

%[y, ~] = rk4(@dCa, @dT, Ca, T, step);

% wektory stanu oraz sterowan
u0 = [CAin Fc]; 
%dy = zeros(size(y)); y_free = zeros(size(y)); du = zeros(size(y));
N=30;Nu=25;ny=2;nu=2;   % horyzont predykcji; horyzont sterowanie;rozmiar wyjść; sterowań
Y_zad=[Ca T];

[M,CtAt,CtV]=MPCSmatrices(A,B,C,N,Nu);

psi=eye(N*ny); 
lambda1= 0.001;
lambda2= 0.001;
lambda=[lambda1 lambda2];
lambda=diag(repmat(lambda,1,Nu));

Kmat=(M'*psi*M+lambda)^-1*M'*psi; 
K1=Kmat(1:nu,:);

time=0:Ts:Tsim;
U=zeros(nu,1); 
v=[0;0]; 
y=zeros(ny,1);
x=[0 0; 0 0];
iter = 2:Tsim/Ts;
Yzad=zeros(N*ny,length(iter));

Y_zad=[Ca T];
%Yzad(1:2:end)=Y_zad(1);
%Yzad(2:2:end)=Y_zad(2);

for iter = 2:Tsim/Ts-1
        %% zmiana wartości zadanych
    if(iter<0.1*length(time))
        Yzad(:,iter)=[ones(N,1)*0; ones(N,1)*0];
    elseif(iter<(1/3)*length(time)) 
        yzad=[0; 5]; 
        Yzad(1:2:end,iter)=yzad(1);
        Yzad(2:2:end,iter)=yzad(2);
    elseif(iter<0.6*length(time)) 
        yzad=[0.01; 5]; 
        Yzad(1:2:end,iter)=yzad(1);
        Yzad(2:2:end,iter)=yzad(2);
    elseif(iter<0.8*length(time)) 
        yzad=[0.01; -5]; 
        Yzad(1:2:end,iter)=yzad(1);
        Yzad(2:2:end,iter)=yzad(2);
	elseif(iter<0.9*length(time)) 
        yzad=[-0.01; -5]; 
        Yzad(1:2:end,iter)=yzad(1);
        Yzad(2:2:end,iter)=yzad(2);
    else
        yzad=[-0.01; -5];
        Yzad(1:2:end,iter)=yzad(1);
        Yzad(2:2:end,iter)=yzad(2);
    end
    
    v(:,iter)=x(:,iter)-(A*x(:,iter-1)+B*U(:,iter-1));  
    deltaU=K1*(Yzad(:,iter)-CtAt*x(:,iter)-CtV*B*U(:,iter-1)-CtV*v(:,iter));
    U(:,iter)=U(:,iter-1)+deltaU;
    
    U(1,iter)=min(U(1,iter),4 - CAin); %ograniczenie na Cain
    U(1,iter)=max(U(1,iter),0 - CAin);
    
    U(2,iter)=min(U(2,iter),30-Fc); %ograniczenia na Fc;
    U(2,iter)=max(U(2,iter),5-Fc);
    
    y(:,iter)=C*x(:,iter)+D*U(:,iter);
    x(:,iter+1)=A*x(:,iter)+B*U(:,iter);
    % wymagane ograniczenie na Ca0
    %x(1,iter+1)=max(x(1,iter+1),-Ca0); 
end

iter = 2:Tsim/Ts;

figure;
subplot(2,1,1);
grid on;
hold on;
plot(iter,y(1,:)+Ca0);

subplot(2,1,2);
grid on;
hold on;
plot(iter,y(2,:)+T0);

figure;
subplot(2,1,1);
grid on;
hold on;
stairs(iter,U(1,:));

subplot(2,1,2);
grid on;
hold on;
stairs(iter,U(2,:));
