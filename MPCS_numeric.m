clear all;
% close all;

global Ca0 T0 Tin CAin0 Fc0;
%%
N=3 %horyzont predykcji
Nu=1 ; %horyzont sterowania
ny=2; %ilosc wyjsc
nu=2; %ilosc sterowań
sim=1;

Ts=0.01; %czas probkowania

Tsim=5; %czas symulacji

[A, B, C, D] = import_ABCD();

contin=ss(A,B,C,D);
disc=c2d(contin,Ts);
A=disc.A;
B=disc.B;
C=disc.C;
D=disc.D;

[M,CtAt,CtV]=MPCSmatrices(A,B,C,N,Nu); % wyznaczenie macierzy MPCS

psi=eye(N*ny); %macierz psi, czyli wartosciowanie wyjsc

lambda1=0.001;
lambda2=0.001;
lambda=[lambda1 lambda2];
lambda=diag(repmat(lambda,1,Nu));

K=(M'*psi*M+lambda)\M'*psi; %macierz K
K1=K(1:nu,:); %nu pierwszych wierszy K

options_quad = optimoptions('quadprog','Display','off');

H = 2*(M'*psi*M+lambda);
J= tril(repmat(eye(nu),Nu));
%%
Aopt = [-J;J;-M;M];

time=0:Ts:Tsim;

v=[0;0]; %bład

U=zeros(nu,N*nu);
du = zeros(nu,N*nu);
Y0=zeros(ny*N, 1);
Y_pred=zeros(ny,ny*N);
% % 
Y_zad = repmat([Ca0; T0/100], N,1);
% % 
x = zeros(ny, ny*N);
Tin_scaled = Tin/100;
% x(:,1)=[CAin0;Tin_scaled];
% x(:,2) = [CAin0; Tin_scaled];
x = repmat([CAin0;Tin_scaled],1,N*ny);
% 
Umin=[0+CAin0; 0+Fc0];
Umax=[4-CAin0;10-Fc0];
Umin=repmat(Umin, Nu,1);
Umax=repmat(Umax, Nu,1);

Ymin=[0+Ca0; -50+Tin_scaled];
Ymax=[5-Ca0; 80-Tin_scaled];
Ymin=repmat(Ymin, N,1);
Ymax=repmat(Ymax, N,1);

%%
for k=2:Tsim/Ts
    
%      % zmiana wartości zadanych
%     if(k<0.1*length(time))
%         Yzad(:,k)=[ones(N,1)*0; ones(N,1)*0];
%     elseif(k<(1/3)*length(time)) 
%         yzad=[0; 1]; 
%         Yzad(1:2:end,k)=yzad(1);
%         Yzad(2:2:end,k)=yzad(2);
%     elseif(k<0.6*length(time)) 
%         yzad=[0.01; 1]; 
%         Yzad(1:2:end,k)=yzad(1);
%         Yzad(2:2:end,k)=yzad(2);
%     elseif(k<0.8*length(time)) 
%         yzad=[0.01; -1]; 
%         Yzad(1:2:end,k)=yzad(1);
%         Yzad(2:2:end,k)=yzad(2);
% 	elseif(k<0.9*length(time)) 
%         yzad=[-0.01; -1]; 
%         Yzad(1:2:end,k)=yzad(1);
%         Yzad(2:2:end,k)=yzad(2);
%     else
%         yzad=[-0.01; -1];
%         Yzad(1:2:end,k)=yzad(1);
%         Yzad(2:2:end,k)=yzad(2);
%     end
    
    v(:,k) = x(:,k) - A*x(:,k-1) - B*U(:,k-1); %wyznaczanie błędu
    Y0(:,k) = CtAt*x(:,k)+CtV*B*U(:,k-1) + CtV*v(:,k);
    f = -2*M'*psi*(Y_zad-Y0(:,k));
    U_prev=repmat(U(:,k-1), Nu,1);
    b = [-Umin+U_prev;Umax-U_prev;-Ymin+Y0(:,k);Ymax-Y0(:,k)];
    lb = repmat([-1.5;-100], Nu,1);
    up = repmat([1.5;100], Nu,1);
    
    deltaU = quadprog(H,f,Aopt,b,[],[],lb,up,[],options_quad);

    deltaU = deltaU(1:nu);
    
    U(:,k)=U(:,k-1) + deltaU;
  
%     if(sim==1)
    % SYMULACJA LINIOWA 
    Y_pred(:,k)=C*x(:,k)+D*U(:,k);
    x(:,k+1)=A*x(:,k)+B*U(:,k);
%     x(:,k+1)=max(x(:,k+1),-Ca0); %ograniczenie na Ca0 - nie moze byc ujemne
%     elseif(sim==2)
%     % SYMULACJA NIELINIOWA     
%         Cain=U(1,k)+Cain0;
%         Fc=U(2,k)+Fc0;
%         
%         Ca(k) = Ca(k-1) - Ts*(Ca(k-1) - Cain + 10000000000*Ca(k-1)*exp(-83301/(10*T(k-1))));
%         T(k) = T(k-1) - Ts*(T(k-1) - 1300000000000*Ca(k-1)*exp(-83301/(10*T(k-1))) + (129*Fc^(3/2)*(T(k-1) - 310))/(250*(Fc + (129*Fc^(1/2))/500)) - 343);
%         y(:,k)=[Ca(k)-Ca0;T(k)-T0];
% 
%         x(:,k+1)=A*y(:,k)+B*U(:,k); %obliczamy stan z wyjscia, a poniewaz C=1, a D=0 to y=x;
%         x(1,k+1)=max(x(1,k+1),-Ca0); %ograniczenie na Ca0 - nie moze byc ujemne
%     end
%     
%   %% analityczny 
%   va(:,k)=xa(:,k)-(A*xa(:,k-1)+B*Ua(:,k-1));  
%     deltaUa=K1*(Yzad(:,k)-CtAt*xa(:,k)-CtV*B*Ua(:,k-1)-CtV*va(:,k));
%   
%     Ua(:,k)=Ua(:,k-1)+deltaUa;
%     
%     %ograniczenia 
%     Ua(1,k)=min(Ua(1,k),3-Cain0); %ograniczenie na Cain
%     Ua(1,k)=max(Ua(1,k),0-Cain0);
%     
%     Ua(2,k)=min(Ua(2,k),30-Fc0); %ograniczenia na Fc;
%     Ua(2,k)=max(Ua(2,k),5-Fc0);
%     
    
%     if(sim==1)
%     %% SYMULACJA LINIOWA 
%         ya(:,k)=C*xa(:,k)+D*Ua(:,k);
%         xa(:,k+1)=A*xa(:,k)+B*Ua(:,k);
%         xa(1,k+1)=max(xa(1,k+1),-Ca0); %ograniczenie na Ca0 - nie moze byc ujemne
%     elseif(sim==2)
%     %% SYMULACJA NIELINIOWA     
%         Cain=Ua(1,k)+Cain0;
%         Fc=Ua(2,k)+Fc0;
%         
%     %% PYTANIE DO ZESPO�U CZY ZAJMUJEMY SI� NIELINIOW� I CZY KOSRZYSTAMY Z NASZEJ IMPLEMENTACJI CZY Z ICH
%         Ca(k) = Ca(k-1) - Ts*(Ca(k-1) - Cain + 10000000000*Ca(k-1)*exp(-83301/(10*T(k-1))));
%         T(k) = T(k-1) - Ts*(T(k-1) - 1300000000000*Ca(k-1)*exp(-83301/(10*T(k-1))) + (129*Fc^(3/2)*(T(k-1) - 310))/(250*(Fc + (129*Fc^(1/2))/500)) - 343);
%         ya(:,k)=[Ca(k)-Ca0;T(k)-T0];
% 
%         xa(:,k+1)=A*ya(:,k)+B*Ua(:,k); %obliczamy stan z wyjscia, a poniewaz C=1, a D=0 to y=x;
%         xa(1,k+1)=max(xa(1,k+1),-Ca0); %ograniczenie na Ca0 - nie moze byc ujemne
%     end

end
%%
% Ca_free = reshape(Y0(1,1,:),[],1);
% T_free = reshape(Y0(2,1,:),[],1);
% Ca_pred = reshape(Y_pred(1,1,:),[],1);
% T_pred = reshape(Y_pred(2,1,:),[],1);
Y_zad = repmat(Y_zad,size(time));
figure;
subplot(1,2,1);
plot(time(2:end), Y0(1,:), time, Y_zad(1,:), time(2:end), Y_pred(1,:));
legend('Ca swobodne', 'Wartosc zadana', 'Ca przewidziane');
subplot(1,2,2);
plot(time(2:end), Y0(2,:), time, Y_zad(2,:), time(2:end), Y_pred(2,:));
legend("T swobodne", "Wartosc zadana", 'T przewidziane');



% Err1 = immse(y(1,:),Yzad(1,:));
% Err2 = immse(y(2,:),Yzad(2,:)); 
% 
% Err3 = immse(ya(1,:),Yzad(1,:));
% Err4 = immse(ya(2,:),Yzad(2,:)); 
% 
% set(groot, 'DefaultTextInterpreter', 'latex')
% set(groot, 'DefaultLegendInterpreter', 'latex')
% 
% figure(wyjscia_fig);
% subplot(2,1,1);
% grid on;
% hold on;
% plot(time,ya(1,:)+Ca0,'DisplayName','$C_\mathrm{a}$');
% plot(time,y(1,:)+Ca0,'DisplayName','$C_\mathrm{a}$');
% stairs(time,Yzad(1,:)+Ca0,'--','DisplayName','$C_\mathrm{a}^\mathrm{zad}$');
% xlabel('t [s]');
% ylabel('$C_{A}$','Interpreter','Latex');
% legend('Wyjscie obiektu (analityczny)','Wyjscie obiektu (numeryczny)', 'Wartosc zadana');
% legend;
% tit=sprintf('$ C_\\mathrm{a}$ $N=%0.1f$ $Nu=%0.1f$ $\\Lambda_1=%0.3f$ $\\Lambda_2=%0.3f$ $MSE_{num}=%f$ $MSE_{anal}=%f$',N,Nu,lambda1,lambda2,Err1,Err3);
% title(tit);
% 
% 
% subplot(2,1,2);
% grid on;
% hold on;
% plot(time,ya(2,:)+T0,'DisplayName','$T$');
% plot(time,y(2,:)+T0,'DisplayName','$T$');
% stairs(time,Yzad(2,:)+T0,'--','DisplayName','$T^\mathrm{zad}$');
% xlabel('t [s]');
% ylabel('$T$','Interpreter','Latex');
% legend('Wyjscie obiektu (analityczny)','Wyjscie obiektu (numeryczny)', 'Wartosc zadana');
% tit=sprintf('Przebiegi $ T $ $MSE_{num}=%f$ $MSE_{anal}=%f$',Err2,Err4);   
% title(tit);
% filename=sprintf('img/%d_%d_%f_%f_wyj_vs',N,Nu,lambda1,lambda2);
% filename=strrep(filename,'.','_');
% filename=sprintf('%s.png',filename);
% %saveas(gcf,filename,'png');
% set(gcf,'PaperPositionMode','auto')
% print(gcf,filename,'-dpng','-r300'); 
% 
% figure(sterowania_fig);
% subplot(2,1,1);
% grid on;
% hold on;
% stairs(time,Ua(1,:)+Cain0,'DisplayName','$C_\mathrm{ain}$');
% stairs(time,U(1,:)+Cain0,'DisplayName','$C_\mathrm{ain}$');
% xlabel('t [s]');
% ylabel('$C_{Ain}$','Interpreter','Latex');
% legend('Analityczny','Numeryczny');
% tit=sprintf('Przebiegi sterowania $ C_\\mathrm{ain} $ $N=%0.1f$ $Nu=%0.1f$ $\\Lambda=%0.3f$ $\\Lambda_2=%0.3f$',N,Nu,lambda1,lambda2);
% title(tit);
% 
% subplot(2,1,2);
% grid on;
% hold on;
% stairs(time,Ua(2,:)+Fc0,'DisplayName','$F_\mathrm{c}$');
% stairs(time,U(2,:)+Fc0,'DisplayName','$F_\mathrm{c}$');
% xlabel('t [s]');
% ylabel('$F_{c}$','Interpreter','Latex');
% legend('Analityczny','Numeryczny');
% tit=sprintf('Przebiegi sterowania $ F_\\mathrm{c}$');
% title(tit);
% filename=sprintf('img/%d_%d_%f_%f_ster_vs',N,Nu,lambda1,lambda2);
% filename=strrep(filename,'.','_');
% filename=sprintf('%s.png',filename);
% %saveas(gcf,filename,'png');
% set(gcf,'PaperPositionMode','auto')
% print(gcf,filename,'-dpng','-r300'); 
