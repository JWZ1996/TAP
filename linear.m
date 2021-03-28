syms s CAint(t) Fct(t) Cat(t) Tt(t) Ca g m K mol kmol cal ro cp k0 E_R h a b ro cp k0 E_R h a b V F Fin CAin Fc Tin Tcin T V0 Fin0 CAin0 Fc0 Tin0 Tcin0 Ca0 T0;

% TODO - dopytac sie piotera czy napewno ok

eq1 =        V * diff(Ca,t) == (Fin*CAin - F*Ca - V*k0*exp(-E_R/T)*Ca);
eq2 = (V*ro*cp) * diff(T,t) == (Fin*ro*cp*Tin - Fin*ro*cp*T + V*h*k0*exp(-E_R/T)*Ca - (a*Fc^(b+1)/(Fc+(a*Fc^b/(2*ro*cp))))*(T-Tin));

% Stad w punkcie pracy (diff(Ca,t) == 0 oraz diff(T,t) == 0):
eq1 = (Fin*CAin - F*Ca - V*k0*exp(-E_R/T)*Ca)/V; % == 0
eq2 = (Fin*ro*cp*Tin - F*ro*cp*T + V*h*k0*exp(-E_R/T)*Ca - (a*Fc^(b+1)/(Fc+(a*Fc^b/(2*ro*cp))))*(T-Tin))/(V*ro*cp); % == 0

% y(t) = [Ca T]    -> wejscie
% u(t) = [CAin Fc] -> wyjscie

lin1 = taylor(eq1, [CAin Fc Ca T], [CAin0 Fc0 Ca0 T0],'Order',2);
lin2 = taylor(eq2, [CAin Fc Ca T], [CAin0 Fc0 Ca0 T0],'Order',2);

% Po zlinearyzowaniu:
eq1 =        diff(Ca,t) == lin1; % (Fin*(CAin - CAin0))/V - (Ca0*F - CAin0*Fin + Ca0*V*k0*exp(-E_R/T0))/V - ((Ca - Ca0)*(F + V*k0*exp(-E_R/T0)))/V - (Ca0*E_R*k0*exp(-E_R/T0)*(T - T0))/T0^2
eq2 =        diff(T,t) ==  lin2; % (h*k0*exp(-E_R/T0)*(Ca - Ca0))/(cp*ro) - ((T - T0)*((a*exp(log(Fc0)*(b + 1)))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro)) + F*cp*ro - (Ca0*E_R*V*h*k0*exp(-E_R/T0))/T0^2))/(V*cp*ro) - ((a*exp(log(Fc0)*(b + 1))*(T0 - Tin))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro)) + F*T0*cp*ro - Fin*Tin*cp*ro - Ca0*V*h*k0*exp(-E_R/T0))/(V*cp*ro) + (a*(Fc - Fc0)*(T0 - Tin)*((exp(log(Fc0)*(b + 1))*((a*b*exp(b*log(Fc0)))/(2*Fc0*cp*ro) + 1))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro))^2 - (exp(log(Fc0)*(b + 1))*(b + 1))/(Fc0*(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro)))))/(V*cp*ro)

