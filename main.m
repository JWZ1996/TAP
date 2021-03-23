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
