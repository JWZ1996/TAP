syms Ca g m K mol kmol cal ro cp k0 E_R h a b ro cp k0 E_R h a b V Fin CAin Fc Tin Tcin T;

% TODO - dopytac sie piotera czy napewno ok

f1 = -V*k0*exp(-E_R/T)*Ca;
f2 = V*h*k0*exp(-E_R/T)*Ca - (a*Fc^(b+1)/(Fc+(a*Fc^b/(2*ro*cp))))*(T-Tin);

ptDefaultF1 = [1 0.16 405];
ptDefaultF2 = [1 0.16 405 310 15];

lin1 = taylor(f1,[V Ca T],        ptDefaultF1,'Order',2);
lin2 = taylor(f2,[V Ca T Tcin Fc],ptDefaultF2,'Order',2);

lin1Simple = simplify(lin1);
lin2Simple = simplify(lin2);