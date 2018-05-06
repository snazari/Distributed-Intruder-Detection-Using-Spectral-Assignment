%% PentaNet: a graph on five vertices conjured up by me to test the 
%  directional residual generator of Speyer and White.
%  Sam Nazari
%  30-August-2015
clear,
clc

%% Graph Structure
Ag = [0 1 1 1 1;
     1 0 1 1 0;
     1 1 0 0 0;
     1 1 0 0 1;
     1 0 0 1 0];
 
Deg = eye(5);
for i = 1:5
    Deg(i,i) = sum(Ag(i,:));
end

C = eye(5)

L = Deg-Ag

x0 = [0 1 1 0 0]

A = -L

% Condition: Distinct Eigenvalues
eig(A)

%% Construct Fault Vectors
f1 = [1 0 0 0 0]'  % Vertex one is the intruder
f2 = [0 1 0 0 0]'  % Vertex two is the intruder
f3 = [0 0 1 0 0]'  % Vertex three is the intruder
f4 = [0 0 0 1 0]'  % Vertex four is the intruder
f5 = [0 0 0 0 1]'  % Vertex five is the intruder
E = [f1 f2 f3 f4 f5]
E = [f2 f3 f4 f5]

% Choose the agent to be attacked
flt1  = 0
flt2  = 1
flt3  = 0
flt4  = 1
flt5  = 0

% Choose a magnitude for the attack
f1Val = 10
f2Val = 10
f3Val = 10
f4Val = 10
f5Val = 10

% Chose the attack time
tf1   = 12
tf2   = 2
tf3   = 1
tf4   = 3
tf5   = 4
%% Step 1: Construct H, T and A1

H = E*pinv(C*E)
T = eye(5)-H*C
A1= A-H*C*A

%% Step 2: Obtain M1..M5

% Check to see if the faults are output seperable rank(C*E)=dim(A)
rank(C*E)

c1 = (eye(5)-C*f1*pinv(C*f1)*C)
c2 = (eye(5)-C*f2*pinv(C*f2)*C)
c3 = (eye(5)-C*f3*pinv(C*f3)*C)
c4 = (eye(5)-C*f4*pinv(C*f4)*C)
c5 = (eye(5)-C*f5*pinv(C*f5)*C)

k1 = A1*(eye(5)-f1*pinv(C*f1)*C)
k2 = A1*(eye(5)-f2*pinv(C*f2)*C)
k3 = A1*(eye(5)-f3*pinv(C*f3)*C)
k4 = A1*(eye(5)-f4*pinv(C*f4)*C)
k5 = A1*(eye(5)-f5*pinv(C*f5)*C)

m1 = [c1;c1*k1;c1*k1^2;c1*k1^3;c1*k1^4]
m2 = [c2;c2*k2;c2*k2^2;c2*k2^3;c2*k2^4]
m3 = [c3;c3*k3;c3*k3^2;c3*k3^3;c3*k3^4]
m4 = [c4;c4*k4;c4*k4^2;c4*k4^3;c4*k4^4]
m5 = [c5;c5*k5;c5*k5^2;c5*k5^3;c5*k5^4]

v1 = 5-rank(m1)
v2 = 5-rank(m2)
v3 = 5-rank(m3)
v4 = 5-rank(m4)
v5 = 5-rank(m5)
v  = v1+v2+v3+v4+v5

% D will be the same as K1.  D will be computed by hand..
% note that the subpace for each fault actually corresponds to the fault
% vector in this example, ie: v1 = f1 = [1 0 0 0 0]' etc..
% D = 2*eye(5)
D = [-2 1 1 1 0;
    1 -1 1 1 0;
    1 1 0 0 0;
    1 1 0 -1 1;
    1 0 0 1 0]

K1 = D
F  = A1-K1*C
K2 = F*H

K = K1+K2

%% Compute projections

CV = C*E
pCV = CV*inv(CV'*CV)*CV'

%% Run model
sim('pentaNetMDL')

%% plot
figure,

subplot(511)
plot(tout,resf1),grid on,ylim([0 50]),title('Residual Signals')
ylabel('f_1')

subplot(512)
plot(tout,resf2),grid on,ylim([0 50])
ylabel('f_2')

subplot(513)
plot(tout,resf3),grid on,ylim([0 50])
ylabel('f_3')

subplot(514)
plot(tout,resf4),grid on,ylim([0 50])
ylabel('f_4')

subplot(515)
plot(tout,resf5),grid on,ylim([0 50])
ylabel('f_5')
xlabel('Time (sec)')