
syms s Ca g m K mol kmol cal ro cp k0 E_R h a b ro cp k0 E_R h a b V Fin CAin Fc Tin Tcin T;
E_R = 8330;
% 1. Równania nieliniowe
eq1 =        V * diff(Ca,t) == (Fin*CAin - Fin*Ca - V*k0*exp(-E_R/T)*Ca);
eq2 = (V*ro*cp) * diff(T,t) == (Fin*ro*cp*Tin - Fin*ro*cp*T + V*h*k0*exp(-E_R/T)*Ca - (a*Fc^(b+1)/(Fc+(a*Fc^b/(2*ro*cp))))*(T-Tin));

% 2. Równania liniowe
eq1lin =        V * diff(Ca,t) - Fin*CAin - Fin*Ca - -(k0*exp(-E_R/405)*(4100625*Ca - 1620*E_R + 656100*V + 4*E_R*T - 656100))/4100625;
eq2lin = (V*ro*cp) * diff(T,t) - (Fin*ro*cp*Tin - Fin*ro*cp*T + ((4*E_R*h*k0*exp(-E_R/405))/4100625 - (2*15^(b + 1)*a*cp*ro)/(30*cp*ro + 15^b*a))*(T - 405) + (4*h*k0*exp(-E_R/405))/25 + h*k0*exp(-E_R/405)*(Ca - 4/25) + (4*h*k0*exp(-E_R/405)*(V - 1))/25 + (2*15^(b + 1)*a*cp*ro*(Tin - 405))/(30*cp*ro + 15^b*a) + (2*15^b*a*cp*ro*(15^b*a + 30*b*cp*ro)*(Fc - 15)*(Tin - 405))/(30*cp*ro + 15^b*a)^2);

linearLaplace = laplace([eq1lin eq2lin],[Ca T CAin Fc],s);

linearMatrix = equationsToMatrix(linearLaplace, [CAin Fc]);