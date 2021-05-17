function [] = MPCS_numeric(Nu, N, Ts, Tsim, lambda1, lambda2, psi1, psi2)

global Ca0 T0 Tin CAin0 Fc0;
% N=20; %horyzont predykcji
% Nu=15 ; %horyzont sterowania
ny=2; %ilosc wyjsc
nu=2; %ilosc sterowań
% 
% Ts=0.01; %czas probkowania
% Tsim=30; %czas symulacji
time=0:Ts:Tsim-Ts;

[A, B, C, D] = import_ABCD();
% dyskretyzacja parametrow 
contin=ss(A,B,C,D);
disc=c2d(contin,Ts);
A=disc.A;
B=disc.B;
C=disc.C;
D=disc.D;

[M,CtAt,CtV]=MPCSmatrices(A,B,C,N,Nu); % wyznaczenie macierzy MPCS
psi = [psi1 psi2];
psi=diag(repmat(psi,1,N));

% psi=eye(N*ny); %macierz psi, czyli wartosciowanie wyjsc
% 
% lambda1=0.01;
% lambda2=0.01;
lambda=[lambda1 lambda2];
lambda=diag(repmat(lambda,1,Nu));

K=(M'*psi*M+lambda)^(-1)*M'*psi; %macierz K
K1=K(1:nu,:); %nu pierwszych wierszy K

options_quad = optimoptions('quadprog','Display','off');

H = 2*(M'*psi*M+lambda);
J= tril(repmat(eye(nu),Nu));
Aopt = [-J;J;-M;M];

v=[0;0]; %wektor zaklocen stanu
% macierze sterowan (nu x Nu) i przyrostow sterowan (nu x Nu)
U=zeros(nu,Nu);  
% du = zeros(nu,Nu);

% wektor trajektorii swobodnej wyjsc
Y0=zeros(ny*N, 1);

% wektor przewidywanych wartosci wyjsc regulowanych
y=zeros(ny,1);

Y_zad = zeros(N*ny);
x = zeros(ny, N);

% ograniczenia na wejscia sterowane
% U = [Cain; Fc]
Umin= repmat([1-CAin0; 0-Fc0],Nu,1);
Umax= repmat([3-CAin0;30-Fc0], Nu,1);

% ograniczenia na wyjscia regulowane
% Y = [Ca; T]
Ymin=repmat([-0.01; -5], N, 1);
Ymax=repmat([0.01; 5], N, 1);

for iter=2:Tsim/Ts
    % zmiana wartości zadanych
    if(iter<0.1*length(time))
        yzad=[0; 0]; 
        Y_zad(1:2:end,iter)=yzad(1);
        Y_zad(2:2:end,iter)=yzad(2);
    elseif(iter<0.2*length(time)) 
        yzad=[0; -2]; 
        Y_zad(1:2:end,iter)=yzad(1);
        Y_zad(2:2:end,iter)=yzad(2);    
    elseif(iter<0.4*length(time)) 
        yzad=[0.01; -2]; 
        Y_zad(1:2:end,iter)=yzad(1);
        Y_zad(2:2:end,iter)=yzad(2);
    elseif(iter<0.6*length(time)) 
        yzad=[0.01; 2]; 
        Y_zad(1:2:end,iter)=yzad(1);
        Y_zad(2:2:end,iter)=yzad(2);
	elseif(iter<0.8*length(time)) 
        yzad=[-0.01; -5]; 
        Y_zad(1:2:end,iter)=yzad(1);
        Y_zad(2:2:end,iter)=yzad(2);
    else
        yzad=[0; 0]; 
        Y_zad(1:2:end,iter)=yzad(1);
        Y_zad(2:2:end,iter)=yzad(2);
    end    
    v(:,iter) = x(:,iter) - A*x(:,iter-1) - B*U(:,iter-1); %wyznaczanie zaklocen stanu
    Y0(:,iter) = CtAt*x(:,iter)+CtV*B*U(:,iter-1) + CtV*v(:,iter); % trajktoria swobodna na podstawie wektora stanu
    
    f = -2*M'*psi*(Y_zad(:,iter)-Y0(:,iter));
    
    % wyznaczanie wektorow ograniczen online
    U_prev=repmat(U(:,iter-1),Nu,1);
    b = [-Umin+U_prev;Umax-U_prev;-Ymin+Y0(:,iter);Ymax-Y0(:,iter)];
    % ograniczenia na wartosci skokow sterowan
    lb = repmat([-2;-15], Nu,1);
    ub = repmat([3;10], Nu,1);
    
    [du, fval, exitflag, output] = quadprog(H,f,Aopt,b,[],[],lb,ub,[],options_quad);

    if isempty(du)
        break;
    else
        U(:,iter)=U(:,iter-1) + du(1:nu);
        y(:,iter)=C*x(:,iter)+D*U(:,iter);
        x(:,iter+1)=A*x(:,iter)+B*U(:,iter) + v(:,iter);
    end

end

subplot(2,2,1);
grid on;
hold on;
plot(time,y(1,:)+Ca0);
title('Wyjście 1 - stężenie substancji A')
xlabel('Czas symulacji [min]')
ylabel('Ca [kmol/m3]') 

subplot(2,2,2);
grid on;
hold on;
plot(time,y(2,:)+T0);
title('Wyjście 2 - temperatura')
xlabel('Czas symulacji [min]')
ylabel('T [K]') 


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

