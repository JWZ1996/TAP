syms Ca g m K mol kmol cal ro cp k0 E_R h a b ro cp k0 E_R h a b V Fin CAin Fc Tin Tcin T;

% TODO - dopytac sie piotera czy napewno ok

f1 = V*k0*exp(-E_R/T)*Ca;
f2 = V*h*k0*exp(-E_R/T)*Ca - (a*Fc^(b+1)/(Fc+(a*Fc^b/(2*ro*cp))))*(T-Tin);

ptDefault = [0.16 405];
lin1 = getLinearized(f1,ptDefault);
lin2 = getLinearized(f2,ptDefault);

ptCalculated = []
lin1_2 = getLinearized(f1,ptCalculated);
lin2_2 = getLinearized(f2,ptCalculated);


function lin = getLinearized(f,pkt)
syms Ca g m K mol kmol cal ro cp k0 E_R h a b ro cp k0 E_R h a b V Fin CAin Fc Tin Tcin T;
lin = taylor(f,[Ca T],pkt,'Order',2);
end