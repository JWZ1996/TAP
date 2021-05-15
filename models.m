%close all
%clear all
consts
global g m K kmol cal min ro cp k E_R h a b ro cp k E_R h a b V Fin CAin Fc Tin Tcin Ca T F V0 Fin0 CAin0 Fc0 Tin0 Tcin0 Ca0 F0 T0;

Tcin = 310;
Tin = 353;

A = [ -(F + V*k*exp(-E_R/T0))/V,                                                                                             -(Ca0*E_R*k*exp(-E_R/T0))/T0^2;
(h*k*exp(-E_R/T0))/(cp*ro), -((a*exp(log(Fc0)*(b + 1)))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro)) + F*cp*ro - (Ca0*E_R*V*h*k*exp(-E_R/T0))/T0^2)/(V*cp*ro)];
 
B = [Fin/V,                                                                                                                                                                                                                   0;
    0, (a*(T0 - Tcin)*((exp(log(Fc0)*(b + 1))*((a*b*exp(b*log(Fc0)))/(2*Fc0*cp*ro) + 1))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro))^2 - (exp(log(Fc0)*(b + 1))*(b + 1))/(Fc0*(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro)))))/(V*cp*ro)];
  
C = [1 0; 0 1];
D = zeros(2,2); % Układ z samymi wspolczynnikami rzeczywistymi - macierz zerowa

syms S Z;

% Podstawienie wartosci

% Model ciagly
model=ss(A,B,C,D);
s = tf('s');
Gs = C*(s*eye(size(A))-A)\B + D;

G  = tf(model);

% [y,:] =  rk4(@dCa, @dT, Ca, T, step);
% figure;
% % figure; step(model,'r--', G, 'b.');
% legend('Przestrzeń stanów', 'Transmitancja');
% title('Ciagly model w przestrzeni stanów i transmitancji');
%%
% Model dyskretny
[m1,mtf1,Gz1] = buildDiscreteModel(model, 0.1 * min);
[m2,mtf2,Gz2] = buildDiscreteModel(model,  0.5 * min);
[m3,mtf3,Gz3] = buildDiscreteModel(model,    1 * min);

function [m,mtf,Gz] = buildDiscreteModel(in,Ts)
syms S Z;
    m = c2d(in,Ts);
    mtf = tf(m); 
    Gz  = vpa(m.C*(Z*eye(size(m.A))-m.A)\m.B+m.D,4);
    %figure; step(m, 'r--', mtf, 'b.');
    %txt= sprintf('Dyskretny model w przestrzeni stanów i transmitancji, Ts=%.1f min', Ts);
    %legend('Przestrzeń stanów', 'Transmitancja');
    %title(txt);
end


