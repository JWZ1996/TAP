close all
clear all
consts
syms Ca_sym T_sym
eq1 = CAin*Fin - Ca_sym*F - Ca_sym*V*k*exp(-E_R/T_sym);
eq2 = Fin*Tin*cp*ro - F*T_sym*cp*ro - (Fc^(b + 1)*a*(T_sym - Tin))/(Fc + (Fc^b*a)/(2*cp*ro)) + Ca*V*h*k*exp(-E_R/T_sym);
[Ca0, T0] = solve(eq1,eq2,Ca_sym,T_sym);

