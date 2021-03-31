
% ==== Plik zawierajacy globalne stale =====
global g m K mol kmol cal min ro cp k0 E_R h a b ro cp k0 E_R h a b V Fin CAin Fc Tin Tcin Ca T F V0 Fin0 CAin0 Fc0 Tin0 Tcin0 Ca0 T0;
% ==== Jednostki ====
% TODO ja lubie cukier syntaktyczny, ale jbc sie wyrzuci

g = 1;
m = 1;
K = 1;
mol = 1;
kmol = 1000 * mol;
cal = 1;

% TODO Wszystkie jednostki czasu mamy w minutach 
% Nie do konca mi sie to podoba, do ustalenia
min = 1;

% ==== Stale ====
ro = 1e6 * g/m^3;
cp = 1 * cal/(g*K);
k0 = 1e10 * 1/min;
E_R = 8330.1 * 1/K;
h = 130 * 1e6 * cal/kmol;
a = 0.516*1e6 * cal/(K*m^3);
b = 0.5;

% ==== Punkt pracy ====
V = 1 * m^3;
Fin = 1 * m^3/min;
CAin = 2 * kmol/m^3;
Fc = 15 * m^3/min;
F = Fc;
Tin = 343 * K;
Tcin = 310 * K;
Ca = 0.16 * kmol/m^3;
T = 405 * K;

% ==== Punkt pracy tak jak to powinno byc ====
V0 = 1 * m^3;
Fin0 = 1 * m^3/min;
CAin0 = 2 * kmol/m^3;
Fc0 = 15 * m^3/min;
Tin0 = 343 * K;
Tcin0 = 310 * K;
Ca0 = 0.16 * kmol/m^3;
T0 = 405 * K;





