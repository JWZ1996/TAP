clear all;
close all;

N = 50; Nu = 50; Ts = 0.001; Tsim = 3; l1 = 1; l2 = 1; p1 = 1; p2 = 1;

figure;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure;
N = 30; Nu = 30;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure;
l1 = 0.1; l2 = 0.1;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure;
l1 = 0.01; l2 = 0.01;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure;
l1 = 0.001; l2 = 0.001;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure;
N = 30; Nu = 30;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure;
N = 15; Nu = 10;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure;
N = 10; Nu = 5;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)

figure;
N = 10; Nu = 4;
MPCS_an_func(Nu, N, Ts, Tsim, l1, l2, p1, p2)
