clear all;
% Parameter:
M = 0.5; % kg           mass of the cart
m = 0.5; % kg           mass of the pendulum
b = 0.1; % N/m/sec      coefficient of friction for cart 
l = 0.3; % m            length to pendulum center of mass
I = 0.006; % kg.m^2     mass moment of inertia of the pendulum
g = 9.8; %m.s^2         gravity

% F                     force applied to the cart
% x                     cart position coordinate
% theta                 pendulum angle from vertical 

% Arbeitspunkt: x1=0 x2=0

q = (M+m)*(I+m*l^2)-(m*l)^2;   %simplifies input

num = [m*l/q  0];
den = [1  b*(I+m*l^2)/q  -(M+m)*m*g*l/q  -b*m*g*l/q];
pend=tf(num,den)
% t=0:0.01:5;
% impulse(pend,t)
% axis([0 1 0 60])

d=I*(M+m)+M*m*l^2; %denominator for the A and B matrices
A=[0 1 0 0;m*g*l(M+m)/d 0 0 m*l*b/d;0 0 0 1;-g*m^2*l^2/d 0 0 -b*(I+m*l^2)/d];
 
B=[0;-m*l/d;0;(I+m*l^2)/d];
 
C=[0 0 1 0;1 0 0 0];
 
D=[0;0];
 
system=ss(A,B,C,D)
 
% T=0:0.05:10;
% U=0.2*ones(size(T));
% [Y,T,X]=lsim(G1,U,T);
% plot(T,Y)
% axis([0 2 0 100])

rank(ctrb(A,B))
 rank(obsv(A,C))
K=acker(A,B,5*[-1 -0.5 -0.4 -0.3]);
Cn = [0 0 1 0];
Dn = 0;
V=1/dcgain(ss(A-B*K,B,Cn,Dn));
L = place(A',C',[-30 -31 -32 -33])'
AI=[A zeros(4,1);Cn 0];
BI=[B;0];
CI=[Cn  0];
  
KI=acker(AI,BI,5*[-1 -0.5 -0.4 -0.3 -0.2])
K1=KI(:,1:end-size(Cn,1));
K2=KI(:,end-size(Cn,1)+1:end);

Q = C'*C; 
Q(1,1) = 80;
Q(3,3) = 400;
R = 1; 
Kopt = lqr(A,B,Q,R)
Ac = A-B*Kopt;     %control matrix
Vopt= 1/dcgain(ss(Ac,B,Cn,Dn))


