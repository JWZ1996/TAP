clear all;

syms s Vt(t) Tint(t) Tcint(t) CAint(t) Fct(t) Ft(t) Fint(t) Cat(t) Tt(t) Ca g m K mol kmol cal ro cp k0 E_R h a b ro cp k0 E_R h a b V F Fin CAin Fc Tin Tcin T V0 Fin0 CAin0 Fc0 Tin0 Tcin0 Ca0 T0 Fin0 F0;

%Dziedzina s - do obliczen symbolicznych
syms Vs Tins Tcins CAins Fcs Fs Cas Ts Fins
eq1 =        V * diff(Ca,t) == (Fin*CAin - F*Ca - V*k0*exp(-E_R/T)*Ca);
eq2 = (V*ro*cp) * diff(T,t) == (Fin*ro*cp*Tin - F*ro*cp*T + V*h*k0*exp(-E_R/T)*Ca - (a*Fc^(b+1)/(Fc+(a*Fc^b/(2*ro*cp))))*(T-Tin));

% Stad w punkcie pracy (diff(Ca,t) == 0 oraz diff(T,t) == 0):
eq1 = (Fin*CAin - F*Ca - V*k0*exp(-E_R/T)*Ca); % == 0
eq2 = (Fin*ro*cp*Tin - F*ro*cp*T + V*h*k0*exp(-E_R/T)*Ca - (a*Fc^(b+1)/(Fc+(a*Fc^b/(2*ro*cp))))*(T-Tin)); % == 0

% [Ca T]    -> wyjscie
% [CAin Fc] -> wejscie

lin1 = taylor(eq1, [CAin Fc Ca T Tin Tcin], [CAin0 Fc0 Ca0 T0 Tin0 Tcin0],'Order',2);
lin2 = taylor(eq2, [CAin Fc Ca T Tin Tcin], [CAin0 Fc0 Ca0 T0 Tin0 Tcin0],'Order',2);

% Po zlinearyzowaniu:
eq1 =        diff(Cat(t),t) == lin1/V;
eq2 =        diff(Tt(t),t)  == lin2/(V*ro*cp);



