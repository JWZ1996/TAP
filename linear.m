clear all;

syms A B C D s Vt(t) Tint(t) Tcint(t) CAint(t) Fct(t) Ft(t) Fint(t) Cat(t) Tt(t) Ca g m K mol kmol cal ro cp k E_R h a b ro cp k E_R h a b V F Fin CAin Fc Tin Tcin T V0 Fin0 CAin0 Fc0 Tin0 Tcin0 Ca0 T0 Fin0 F0;

%Dziedzina s - do obliczen symbolicznych
syms Vs Tins Tcins CAins Fcs Fs Cas Ts Fins
eq1 =        V * diff(Ca,t) == (Fin*CAin - F*Ca - V*k*exp(-E_R/T)*Ca);
eq2 = (V*ro*cp) * diff(T,t) == (Fin*ro*cp*Tin - F*ro*cp*T + V*h*k*exp(-E_R/T)*Ca - (a*Fc^(b+1)/(Fc+(a*Fc^b/(2*ro*cp))))*(T-Tin));

% Stad w punkcie pracy (diff(Ca,t) == 0 oraz diff(T,t) == 0):
eq1 = (Fin*CAin - F*Ca - V*k*exp(-E_R/T)*Ca); % == 0
eq2 = (Fin*ro*cp*Tin - F*ro*cp*T + V*h*k*exp(-E_R/T)*Ca - (a*Fc^(b+1)/(Fc+(a*Fc^b/(2*ro*cp))))*(T-Tin)); % == 0

% [Ca T]    -> wyjscie
% [CAin Fc] -> wejscie

lin1 = taylor(eq1, [CAin Fc Ca T], [CAin0 Fc0 Ca0 T0],'Order',2);
lin2 = taylor(eq2, [CAin Fc Ca T], [CAin0 Fc0 Ca0 T0],'Order',2);

lin1 = lin1/V;
lin2 = lin2/(V*ro*cp);

% Po zlinearyzowaniu:
eq1 =        diff(Ca,t) == lin1/V;
eq2 =        diff(Tt,t) == lin2/(V*ro*cp);

% Wymiary macierzy: 2 wejscia, 2 wyjscia 
% n = 2 zmienne zalezne 
% p = 2 wyjscia
% q = 4 wejscia 
% A: nxn = 2x2
% B: nxp = 2x2
% C: qxn = 2x2
% D: qxp = 2x2 

A = sym('A%d%d', [2 2]);
B = sym('B%d%d', [2 2]);
C = [1 0; 0 1];
D = zeros(2,2); % Układ z samymi wspolczynnikami rzeczywistymi - macierz zerowa

A = insertVal( A,1,1,lin1,Ca );
A = insertVal( A,1,2,lin1,T  );
A = insertVal( A,2,1,lin2,Ca );
A = insertVal( A,2,2,lin2,T  );

B = insertVal( B,1,1,lin1,CAin );
B = insertVal( B,1,2,lin1,Fc  );
B = insertVal( B,2,1,lin2,CAin );
B = insertVal( B,2,2,lin2,Fc  );

% Macierz transmitancji
%
% 1. Bierzecie macierze A i B i podstawiacie przy uzyciu subs albo
% korzystacie z nich w trybie niesymbolicznym - wtedy wystarczy wkleic ich
% zawartosc w kod
% 2. ss2tf zwroci maciesz transmitancji - koniec
% output = ss2tf(A,B,C,D);
%
%
% Model dyskretny - dyskretna macierz transmitancji
% 1. Bierzecie macierze A i B i podstawiacie przy uzyciu subs albo
% korzystacie z nich w trybie niesymbolicznym - wtedy wystarczy wkleic ich
% zawartosc w kod
% 2. Kod powinien wygladac tak:
% 
% Ts to czas probkowania
% d = c2d(ss(A,B,C,D),Ts);
% output = ss2tf(d.A,d.B,d.C,d.D);
% jesli chodzi o odpowiedzi impulsowe:
% ss2tf(A,B,C,D,N) zwroci model z odpowiedzia na N-tym wejsciu (Kuba :))




