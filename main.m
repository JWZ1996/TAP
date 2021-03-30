% ==== Plik gÅ‚owny =====
global g m K mol kmol cal min ro cp k0 E_R h a b ro cp k0 E_R h a b V Fin CAin Fc Tin Tcin Ca T;
consts
% === Symulacja obiektu - RK4 ===
step = 0.01;
Ca = 0.16 * kmol/m^3;
T =  405 * K;
result = rk4(Ca,T,step);

function [y] = rk4(Ca,T,step)
global g m K mol kmol cal min ro cp k0 E_R h a b ro cp k0 E_R h a b V Fin CAin Fc Tin Tcin;
t=0:step:20;
y(:,1) = [Ca T];
%wyliczanie wspo³czynników
for i=1:(length(t)-1)
    k11=dCa(Ca,T);
    k12=dT(Ca,T);

    k21=dCa(Ca+0.5*step,T+0.5*step*k11);
    k22=dT(Ca+0.5*step,T+0.5*step*k12);

    k31=dCa(Ca+0.5*step,T+0.5*step*k21);
    k32=dT(Ca+0.5*step,T+0.5*step*k22);

    k41=dCa(Ca+step,T+step*k31);
    k42=dT(Ca+step,T+step*k32);

    Ca=Ca+(step/6)*(k11+k41+2*(k21+k31));
    T=T+(step/6)*(k12+k42+2*(k22+k32));
    y(:,i+1)=[Ca T];
end
plot(y(1,:),y(2,:));
title('Plot 1')
end

function [dCa] = dCa(Ca, T)
global g m K mol kmol cal min ro cp k0 E_R h a b ro cp k0 E_R h a b V Fin CAin Fc Tin Tcin;
dCa = (Fin*CAin - Fin*Ca - V*k0*exp(-E_R/T)*Ca)/V;
end

function [dT] = dT(Ca,T)
global g m K mol kmol cal min ro cp k0 E_R h a b ro cp k0 E_R h a b V Fin CAin Fc Tin Tcin;
dT = (Fin*ro*cp*Tin - Fin*ro*cp*T + V*h*k0*exp(-E_R/T)*Ca - (a*Fc^(b+1)/(Fc+(a*Fc^b/(2*ro*cp))))*(T-Tin))/(V*ro*cp);
end
function [dCa] = dCaLin(Ca, T)
global g m K mol kmol cal min ro cp k0 E_R h a b ro cp k0 E_R h a b V Fin CAin Fc Tin Tcin;
dCa = -((Ca - Ca0)*(F + V*k0*exp(-E_R/T0)) + Ca0*F - CAin0*Fin - Fin*(CAin - CAin0) + Ca0*V*k0*exp(-E_R/T0) + (Ca0*E_R*V*k0*exp(-E_R/T0)*(T - T0))/T0^2)/V
end

function [dT] = dTLin(Ca,T)
global g m K mol kmol cal min ro cp k0 E_R h a b ro cp k0 E_R h a b V Fin CAin Fc Tin Tcin;
dT = ((a*exp(log(Fc0)*(b + 1)))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro)) + Fin*cp*ro)*(Tin - Tin0) - (T - T0)*((a*exp(log(Fc0)*(b + 1)))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro)) + F*cp*ro - (Ca0*E_R*V*h*k0*exp(-E_R/T0))/T0^2) + a*(Fc - Fc0)*(T0 - Tin0)*((exp(log(Fc0)*(b + 1))*((a*b*exp(b*log(Fc0)))/(2*Fc0*cp*ro) + 1))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro))^2 - (exp(log(Fc0)*(b + 1))*(b + 1))/(Fc0*(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro)))) - (a*exp(log(Fc0)*(b + 1))*(T0 - Tin0))/(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro)) - F*T0*cp*ro + Fin*Tin0*cp*ro + Ca0*V*h*k0*exp(-E_R/T0) + V*h*k0*exp(-E_R/T0)*(Ca - Ca0)/(V*ro*cp)
end
