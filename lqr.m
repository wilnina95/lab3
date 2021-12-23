clear all; clc;close all;
pkg load control
pkg load signal
Ao1 = 3.16*10^(-5);
Ao2 = 1.78*10^(-5);
g = 9.8;
k2 = 1.39*10^(-4);
k2 = 7.88*10^(-5);
a1 = 0.5; %Vp
k1 = 0.033;
a2 =0.45;
hi1 =(k1*a1/k2)^2; 
hi2 =((k2^2)*hi1)/((k3*a2)^2);
A = [-k2/(2*Ao1*(hi1)^(1/2)) 0;
     k2/(2*Ao2*(hi1)^(1/2)) (-k3*a2)/(2*Ao2*(hi2)^(1/2))];
B = [k1/Ao1 ;
     0];
C = [0 1];
Q = 1*C'*C; %se recomienda
R = 0.8;
k = lqr(A,B,Q,R);

eig(A-B*k) %polos en lazo cerrado
B1 = [k1/Ao1 0;0 0];
C1 = [0 1;0 0];


sys = ss(A-B*k,B1*k(1),C1,[0 0;0 0]);%crea un nuevo sistema 
step(sys)


[b,a]=ss2tf(A-B*k,B*k(1),C,0)
g= tf(b,a);
step(g)
bode(g)
