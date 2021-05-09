function [du, y_free, dy]=predict_output(A, B, C, N, Nu, Y_zad, x, x_prev, u)
[M, CtAt, CtV, At, V] = MPCSmatrices(A, B, C, N, Nu);
[nn nu]=size(M);
ny = nn/N;          % liczba wyjsc regulowanych (Ca i T)

psi = eye(nn);
lambda = eye(nu);
K = [M'*psi*M + lambda]\M'*psi;
K1=K(1:nu, :);      % macierz nu pierwszych podmacierzy
Ke = zeros(nu, ny);
% zmienne pomocnicze - states_part do wyliczania wartosci powiazanych x(k);
% controls_part do wyliczania wartosci powiazanych z u(k-1)
states_part = zeros(nu, 1);
controls_part = zeros(nu, 1);

% https://studia2.elka.pw.edu.pl/file/21L/103C-ARxxx-MSP-TAP/priv/05RegMPCS%5FSlajdy.pdf
% strona 10/47
Vp=eye(nu);
v = x - A*x_prev - B*u;
    for i=1:N
        Ke = Ke + K1(:,(i-1)*ny+1:i*ny);
        K1p = K1(1:nu, (i-1)*ny+1:i*ny);
        A_temp = At((i-1)*ny+1:i*ny, 1:nu);
        states_part = states_part + K1p*C*A_temp*x;
        controls_part = controls_part + K1p*C*Vp*(B*u + v);
        Vp = Vp+A_temp;
    end
    
du = Ke*Y_zad - states_part - controls_part;  

% przy zalozeniu ze Y_zad jest const i zawiera tylko dwie wartosci -
% roszerzenie wektora Y_zad dla celow mnozenia; dwie linijki ponizej sluza
% glownie potwierdzeniu poprawnej implementacji algorytmu
Y_zad = repmat(Y_zad, N, 1);
duu = K1*(Y_zad - CtAt*x - CtV*B*u - CtV*v);    % na to samo wychodzi xD

y_free = CtAt*x + CtV*B*u + CtV*v;
dy = M*du;
dy = dy(1:2);
y_free=y_free(1:2);

end
