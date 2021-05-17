% clear all;
% close all;

N = 50; Nu = 50; Ts = 0.01; Tsim = 30; l1 = 1; l2 = 1; p1 = 1; p2 = 1;

figure('units','normalized','outerposition',[0 0 1 1]) %1
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
MPCS_numeric(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure('units','normalized','outerposition',[0 0 1 1]) %2
N = 30; Nu = 30;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
MPCS_numeric(Nu, N, Ts, Tsim, l1, l2, p1, p2)


figure('units','normalized','outerposition',[0 0 1 1]) %3
l1 = 0.1; l2 = 0.1;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
MPCS_numeric(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure('units','normalized','outerposition',[0 0 1 1]) %4
l1 = 0.01; l2 = 0.01;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
MPCS_numeric(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure('units','normalized','outerposition',[0 0 1 1]) %5
l1 = 0.005; l2 = 0.005;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
MPCS_numeric(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure('units','normalized','outerposition',[0 0 1 1]) %6
N = 20; Nu = 20;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
MPCS_numeric(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure('units','normalized','outerposition',[0 0 1 1]) %7
N = 15; Nu = 10;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
MPCS_numeric(Nu, N, Ts, Tsim, l1, l2, p1, p2)


figure('units','normalized','outerposition',[0 0 1 1]) %8
N = 10; Nu = 3; p1 = 30; p2 = 30; l1 = 0.005; l2 = 0.005;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
MPCS_numeric(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure('units','normalized','outerposition',[0 0 1 1]) %9
N = 10; Nu = 1;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
MPCS_numeric(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure('units','normalized','outerposition',[0 0 1 1]) %10
N = 15; Nu = 3; l1 = 0.001; l2 = 0.001; p1 = 30; p2 = 30;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
MPCS_numeric(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure('units','normalized','outerposition',[0 0 1 1]) %11
N = 30; Nu = 30; l1 = 0.005; l2 = 0.005; p1 = 0.1; p2 = 0,1;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
MPCS_numeric(Nu, N, Ts, Tsim, l1, l2, p1, p2)

h=findobj('type','figure'); % find the handles of the opened figures
folder='D:\Studia\21L\TAP\TAP\obrazki_auto';  % Desination folder
for k=1:numel(h)
  filename=sprintf('fig%d.jpg',numel(h)-k+1);
  file=fullfile(folder,filename);
  saveas(h(k),file);
end
