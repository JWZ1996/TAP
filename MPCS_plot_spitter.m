clear all;
close all;

figure;
Nu = 30; N = 30; Ts = 0.01; Tsim = 30; l1 = 1; l2 = 1; p1 = 1; p2 = 1;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure;
l1 = 0.01; l2 = 0.01;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure;
Nu = 5;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure;
N = 8;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)