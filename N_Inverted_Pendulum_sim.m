%% 
close all
clear all
clc

% Select Number of Bars

N = 2; % Number of Bars (1 cart + N bars)

% Create Nonlinear and Linearized System Matrices
N_Inverted_Pendulum_V3(N);

%% Choose Configuration
% Control variable chooses control scheme
% Control = 0 -> LQR Control
% Control = 1 -> Pole Placement Control
Control = 0;

% Mode variable chooses simulation mode
% Mode = 0 -> Linear System Full-State Feedback 
% Mode = 1 -> Nonlinear System Full-State Feedback
% Mode = 2 -> Linear System Observer-State Feedback 
% Mode = 3 -> Nonlinear System Observer-State Feedback
Mode = 3;

%% System Parameters
l(1) = 0;
l(2) = 1;
l(3) = 1;
for i = 4:(N+1)
    l(i) = 1;
end
m1 = 1e1;
m2 = 1e0;
m3 = 1e0;
m4 = 1e0;
m5 = 1e0;
m6 = 1e0;
m7 = 1e0;
m8 = 1e0;
m9 = 1e0;
m10 = 1e0;
m11 = 1e0;
b(1) = 0.00;
b(2) = 0.000;
for i = 2:(N+1)
    b(i) = b(i-1);
end
g = 9.8;

m = [m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11]';


m = m(1:(N+1));

for i = 1:(N+1)
    I(i) = (1/3)*m(i)*(l(i)^2);
end

p = [m; I'; l'; b'; g];

%% Create State Space Model and Check Controllability

C = [eye(N+1) zeros(N+1)];

sys = ss(A_matrix(p),B_matrix(p),C,zeros(N+1,1))

Ctrb = ctrb(A_matrix(p),B_matrix(p));

rank(Ctrb)

%% Create Controller

if Control == 0
    % LQR Controller
    Q = diag([ones(1,N+1) zeros(1,N+1)]);
    
    R = 1e0;
    
    [K,~,CLP] = lqr(sys,Q,R)
elseif Control == 1
    % Pole Placement Controller
    % Requirements
    Mp = 2/100;
    ep = 2/100;
    ts = 15;
    
    xi = -log(Mp)/sqrt(pi^2 + log(Mp)^2);
    wn = -log(ep)/(xi*ts);
    wd = wn*sqrt(1 - xi^2);
    sigma = xi*wn
    j = sqrt(-1)    
    
    CLP = [(-sigma + j*wd) (-sigma - j*wd) (-10*sigma + j*wd) (-10*sigma - j*wd) (-15*sigma + j*wd) (-15*sigma - j*wd)]'
    
    K = place(A_matrix(p),B_matrix(p),CLP)
end

%% Create Observer
Obsv = obsv(A_matrix(p),C);

rank(Obsv)

L = place(A_matrix(p)',C',(5*real(CLP) + j*imag(CLP)))'

%% Simulation Parameters and Initial Conditions

tspan = [0 15];
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

x0 = [0 (16/180)*pi*ones(1,N) zeros(1,(N+1)) zeros(1,2*(N+1))];

%% Simulate

if Mode == 0
    % Linear Model full-state feedback (observer states not used)
    [t,q] = ode45(@(t,qv) [(A_matrix(p)-B_matrix(p)*K)*qv(1:2*(N+1)); (L*C-B_matrix(p)*K)*qv(1:2*(N+1))+(A_matrix(p)-L*C)*qv(2*(N+1)+1:4*(N+1))], tspan, x0, opts);
elseif Mode == 1
    % Nonlinear Model full-state feedback (observer states not used)
    [t,q] = ode45(@(t,qv) [qv(N+2:2*(N+1)); (Mn(t,qv(1:2*(N+1)),p)^(-1))*Fn(t,qv(1:2*(N+1)),p,-K*qv(1:2*(N+1))); (L*C-B_matrix(p)*K)*qv(1:2*(N+1))+(A_matrix(p)-L*C)*qv(2*(N+1)+1:4*(N+1))], tspan, x0, opts);
elseif Mode == 2
    % Linear Model observer feedback (observer states used)
    [t,q] = ode45(@(t,qv) [A_matrix(p)*qv(1:2*(N+1)) - B_matrix(p)*K*qv(2*(N+1)+1:4*(N+1)); L*C*qv(1:2*(N+1))+(A_matrix(p)-B_matrix(p)*K-L*C)*qv(2*(N+1)+1:4*(N+1))], tspan, x0, opts);
elseif Mode == 3
    % Nonlinear Model observer feedback (observer states used)
    [t,q] = ode45(@(t,qv) [qv(N+2:2*(N+1)); (Mn(t,qv(1:2*(N+1)),p)^(-1))*Fn(t,qv(1:2*(N+1)),p,-K*qv(2*(N+1)+1:4*(N+1))); L*C*qv(1:2*(N+1))+(A_matrix(p)-B_matrix(p)*K-L*C)*qv(2*(N+1)+1:4*(N+1))], tspan, x0, opts);
end

%% Plot Simulation Results

qv = q(:,1:2*(N+1));
qe = q(:,2*(N+1)+1:4*(N+1));

%sim('simulink_model')

% qv = out.yout(:,1:6);
% qe = out.yout(:,7:12);
% t = out.tout;

figure(1)
obj = VideoWriter('System.avi');
Plot_Inverted_Pendulum(t,qv,N,l,obj)
figure(2)
obj = VideoWriter('Observer.avi');
Plot_Inverted_Pendulum(t,qe,N,l,obj)
