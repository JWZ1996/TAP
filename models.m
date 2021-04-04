close all
clear all
consts
global g m K kmol cal min ro cp k E_R h a b ro cp k E_R h a b V Fin CAin Fc Tin Tcin Ca T F V0 Fin0 CAin0 Fc0 Tin0 Tcin0 Ca0 F0 T0;
A = [ -(F + V*k*exp(-E_R/T0))/V ,                                                                           -(Ca0*E_R*k*exp(-E_R/T0))/T0^2;
    (h*k*exp(-E_R/T0))/(cp*ro), -((a*exp(log(Fc0)*(b + 1)))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro)) + F*cp*ro - (Ca0*E_R*V*h*k*exp(-E_R/T0))/T0^2)/(V*cp*ro)];
 
B = [ Fin/V,                                                                                                                                                               0;
    0, (a*(T0 - Tin)*((exp(log(Fc0)*(b + 1))*((a*b*exp(b*log(Fc0)))/(2*Fc0*cp*ro) + 1))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro))^2 - (exp(log(Fc0)*(b + 1))*(b + 1))/(Fc0*(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro)))))/(V*cp*ro)];
C = [1 0; 0 1];
D = zeros(2,2); % Układ z samymi wspolczynnikami rzeczywistymi - macierz zerowa

syms S;

% Model ciagly
model=ss(A,B,C,D);
G1 = tf(model);

% Model dyskretny
Ts = 0.1 * min; 
discreteModel = c2d(model,Ts);
Gd = tf(discreteModel);
