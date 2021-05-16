clear all;
% close all;

global Ca0 T0 Tin CAin0 Fc0;
%%
N=10; %horyzont predykcji
Nu=4 ; %horyzont sterowania
ny=2; %ilosc wyjsc
nu=2; %ilosc sterowań

Ts=0.01; %czas probkowania
Tsim=10; %czas symulacji
time=0:Ts:Tsim;

[A, B, C, D] = import_ABCD();
% dyskretyzacja parametrow 
contin=ss(A,B,C,D);
disc=c2d(contin,Ts);
A=disc.A;
B=disc.B;
C=disc.C;
D=disc.D;

[M,CtAt,CtV]=MPCSmatrices(A,B,C,N,Nu); % wyznaczenie macierzy MPCS

psi=eye(N*ny); %macierz psi, czyli wartosciowanie wyjsc

lambda1=0.01;
lambda2=0.01;
lambda=[lambda1 lambda2];
lambda=diag(repmat(lambda,1,Nu));

K=(M'*psi*M+lambda)\M'*psi; %macierz K
K1=K(1:nu,:); %nu pierwszych wierszy K

options_quad = optimoptions('quadprog','Display','off');

H = 2*(M'*psi*M+lambda);
J= tril(repmat(eye(nu),Nu));
Aopt = [-J;J;-M;M];
%%
T_scale_factor = 100; % wspolczynnik skalujacy temperature dla lepszych wlasnosci numerycznych 
v=[0;0]; %wektor zaklocen stanu

% macierze sterowan (nu x Nu) i przyrostow sterowan (nu x Nu)
U=zeros(nu,Nu);  
du = zeros(nu,Nu);
x = du;

% wektor trajektorii swobodnej wyjsc
Y0=zeros(ny*N, 1);

% wektor przewidywanych wartosci wyjsc regulowanych
Y_pred=zeros(ny,1);

Y_zad = repmat([Ca0; T0/T_scale_factor], N,1);
states = zeros(ny, N);

% ograniczenia na wejscia sterowane
% U = [Cain; Fc]
Umin= repmat([0; 0],Nu,1);
Umax= repmat([40;100], Nu,1);

% ograniczenia na wyjscia regulowane
% Y = [Ca; T]
Ymin=repmat([0; 0], N, 1);
Ymax=repmat([50; 450/T_scale_factor], N, 1);

%%
for k=2:Tsim/Ts
        
    v(:,k) = states(:,k) - A*states(:,k-1) - B*U(:,k-1); %wyznaczanie zaklocen stanu
    Y0(:,k) = CtAt*states(:,k)+CtV*B*U(:,k-1) + CtV*v(:,k); % trajktoria swobodna na podstawie wektora stanu
    
    f = -2*M'*psi*(Y_zad-Y0(:,k));
    
    % wyznaczanie wektorow ograniczen online
    U_prev=repmat(U(:,k-1), Nu,1);
    b = [-Umin+U_prev;Umax-U_prev;-Ymin+Y0(:,k);Ymax-Y0(:,k)];
    % ograniczenia na wartosci skokow sterowan
    lb = repmat([-10;-10], Nu,1);
    ub = repmat([10;10], Nu,1);
    
    [du, fval, exitflag, output] = quadprog(H,f,Aopt,b,[],[],lb,ub,[],options_quad);
    
    if isempty(du)
        break;
    else
        U(:,k)=U(:,k-1) + du(1:nu);
        Y_pred(:,k)=C*states(:,k)+D*U(:,k);
        states(:,k+1)=A*states(:,k)+B*U(:,k) + v(:,k);
    end

end
%%

min_len = min([size(Y_pred,2) size(Y0,2)]);
Y_zad = repmat(Y_zad(1:2),1, min_len);


figure;
subplot(1,2,1);
plot(time(1:min_len), Y0(1,:), time(1:min_len), Y_zad(1,:), time(1:min_len), Y_pred(1,:));
legend('Ca swobodne', 'Wartosc zadana', 'Ca przewidziane');
xlabel('Czas');
ylabel('Stezenie Ca [kmol/m3]')
subplot(1,2,2);
plot(time(1:min_len), Y0(2,:)*T_scale_factor, time(1:min_len), Y_zad(2,:)*T_scale_factor, time(1:min_len), Y_pred(2,:)*T_scale_factor);
legend("T swobodne", "Wartosc zadana", 'T przewidziane');
xlabel('Czas');
ylabel('Temperatura T [K]')


