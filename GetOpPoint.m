close all
clear all
consts

digits(64);
res = fsolve(@system,[Ca0 T0]);

function res = system(CaT)
global g m K kmol cal min ro cp k E_R h a b ro cp k E_R h a b V Fin CAin Fc Tin Tcin Ca T F V0 Fin0 CAin0 Fc0 Tin0 Tcin0 Ca0 F0 T0;

    Ca_sym = CaT(1);
    T_sym  = CaT(2);
    
    res(1) = (Fin*CAin)/V - F*Ca_sym/V-k*exp(-E_R/T_sym)*Ca_sym;
    res(2) = Fin*Tin/V - F*T_sym/V + h*k*exp(-E_R/T_sym)*Ca_sym/(ro*cp) - (a*(Fc)^(b+1))/(Fc+(a*(Fc)^(b))/(2*ro*cp)) * (T_sym-Tcin) / (V*ro*cp);
end