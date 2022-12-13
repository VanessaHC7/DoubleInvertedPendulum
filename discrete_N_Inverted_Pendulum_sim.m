close all
clear
clc
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], ...
    [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], ...
    [0.6350, 0.0780, 0.1840]};

% Set parameters

N = 2; % Number of Bars (1 cart + N bars)

N_Inverted_Pendulum(N);

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

A = A_matrix(p); B = B_matrix(p); C = [eye(N+1) zeros(N+1)]; D = zeros(N+1,1);

% Discretize with zero-order hold
dt = 0.01;
Ad = eye(size(A, 1)) + dt*A;
Bd = dt*B;
i = 2;
while true
    am1 = (dt^i/factorial(i))*A^(i-1);
    addend_A = am1*A;
    addend_B = am1*B;
    if max(abs([addend_A addend_B]), [], 'all') < 1e-10
        break
    end
    Ad = Ad + addend_A;
    Bd = Bd + addend_B;
    i = i + 1;
end

sys = ss(Ad, Bd, C, D, dt)

Ctrb = ctrb(Ad, Bd);

rank(Ctrb)

cond(Ctrb)

Obsv = obsv(Ad, C);

rank(Obsv)

cond(Obsv)

%%
% LQR controller

Q = diag([ones(1,N+1) zeros(1,N+1)]);

R = 1e0;

[K,~,CLP] = dlqr(Ad, Bd, Q, R)
abs(eig(Ad - Bd*K))

%%
% Noise characteristics
Gamma = Bd; % wi enters via the control input
Wi = 0.001; % variance of wi (scalar)
Wo = diag([0.001 0.001*ones(1, size(C, 1)-1)]); % covariance of w0 (vector)
% L = place(Ad.', C.', CLP/5).'

%%
% Discrete-time simulation
tspan = [0 20];
n_t = round((tspan(2) - tspan(1))/dt) + 1;
time = linspace(tspan(1), tspan(2), n_t);
trials = 1000;

x0 = [0 (38/180)*pi*ones(1,N) zeros(1,(N+1))].';
xhat0 = zeros(2*(N+1), 1);
n = length(x0);
m = size(C, 1);

x = zeros(n, n_t, trials);
y = zeros(m, n_t, trials);
u = zeros(n_t, trials);
% xhat = zeros(n, n_t, trials);
xhat_m = zeros(n, n_t, trials);
xhat_p = zeros(n, n_t, trials);
Sigma_m = zeros(n, n, n_t, trials);
Sigma_p = zeros(n, n, n_t, trials);
disp("Starting simulation...")
for i = 1:trials
    x(:, 1, i) = x0;
    y(:, 1, i) = C*x0 + mvnrnd(zeros(m, 1), Wo, 1).';
%     xhat_p(:, 1, i) = xhat0;
    xhat_m(:, 1, i) = xhat0;
    [~,Sigma_m(:, :, 1, i),~] = dlqr(Ad.', C.', Gamma*Wi*Gamma.', Wo);
    for k = 1:n_t-1
        L_k = Sigma_m(:, :, k, i)*C.'*inv(C*Sigma_m(:, :, k, i)*C.' + Wo);
        xhat_p(:, k, i) = xhat_m(:, k, i) + L_k*(y(:, k, i) - C*xhat_m(:, k, i));
        Sigma_p(:, :, k, i) = (eye(n) - L_k*C)*Sigma_m(:, :, k, i);
        u(k, i) = -K*xhat_p(:, k, i);
        x(:, k+1, i) = Ad*x(:, k, i) + Bd*u(k, i) + Gamma*mvnrnd(0, Wi, 1).';
        y(:, k+1, i) = C*x(:, k+1, i) + mvnrnd(zeros(m, 1), Wo, 1).';
%         xhat_p(:, k+1, i) = Ad*xhat_p(:, k, i) + Bd*u(k, i) + L*(y(:, k, i) - C*xhat_p(:, k, i));
        xhat_m(:, k+1, i) = Ad*xhat_p(:, k, i) + Bd*u(k, i);
        Sigma_m(:, :, k+1, i) = Ad*Sigma_p(:, :, k, i)*Ad.' + Gamma*Wi*Gamma.';
    end
    L_end = Sigma_m(:, :, end, i)*C.'*inv(C*Sigma_m(:, :, end, i)*C.' + Wo);
    xhat_p(:, end, i) = xhat_m(:, end, i) + L_end*(y(:, end, i) - C*xhat_m(:, end, i));
    Sigma_p(:, :, end, i) = (eye(n) - L_end*C)*Sigma_m(:, :, end, i);
end
disp("Finished!")

for j = 2:N+1
    x(j, :, :) = x(j, :, :)*(180/pi);
    xhat_p(j, :, :) = xhat_p(j, :, :)*(180/pi);
end

%%
% Plot results

disp("Plotting...")
ylabels = {'x [m]', '\theta_1 [degrees]', '\theta_2 [degrees]'};
xlabels = {'', '', 'Time [s]'};
titles = {'Kalman Filter Initial Condition Response (\sigma_v^2 = 0.001, \sigma_w^2 = 0.001)', '', ''};
figure(1)
for j = 1:N+1
    subplot(N+1, 1, j)
    x_avg = median(x(j, :, :), 3);
%     x_avg = mean(x(j, :, :), 3);
    plot(time, x_avg, 'Color', colors{1})
    hold on
    x_5prct = prctile(x(j, :, :), 5, 3);
    x_25prct = prctile(x(j, :, :), 25, 3);
    x_75prct = prctile(x(j, :, :), 75, 3);
    x_95prct = prctile(x(j, :, :), 95, 3);
    fill([time(1:end-1); time(1:end-1); time(2:end); time(2:end)], ...
        [x_95prct(1:end-1); x_5prct(1:end-1); x_5prct(2:end); x_95prct(2:end)], ...
        colors{1}, 'EdgeColor', 'none', 'FaceAlpha', 0.2)
    fill([time(1:end-1); time(1:end-1); time(2:end); time(2:end)], ...
        [x_75prct(1:end-1); x_25prct(1:end-1); x_25prct(2:end); x_75prct(2:end)], ...
        colors{1}, 'EdgeColor', 'none', 'FaceAlpha', 0.2)
    xhat_avg = median(xhat_p(j, :, :), 3);
    plot(time, xhat_avg, 'Color', colors{2})
    xhat_5prct = prctile(xhat_p(j, :, :), 5, 3);
    xhat_25prct = prctile(xhat_p(j, :, :), 25, 3);
    xhat_75prct = prctile(xhat_p(j, :, :), 75, 3);
    xhat_95prct = prctile(xhat_p(j, :, :), 95, 3);
    fill([time(1:end-1); time(1:end-1); time(2:end); time(2:end)], ...
        [xhat_95prct(1:end-1); xhat_5prct(1:end-1); xhat_5prct(2:end); xhat_95prct(2:end)], ...
        colors{2}, 'EdgeColor', 'none', 'FaceAlpha', 0.2)
    fill([time(1:end-1); time(1:end-1); time(2:end); time(2:end)], ...
        [xhat_75prct(1:end-1); xhat_25prct(1:end-1); xhat_25prct(2:end); xhat_75prct(2:end)], ...
        colors{2}, 'EdgeColor', 'none', 'FaceAlpha', 0.2)
    hold off
    xlabel(xlabels{j})%, 'Interpreter', 'latex')
    ylabel(ylabels{j})%, 'Interpreter', 'latex')
    title(titles{j})%, 'Interpreter', 'latex')leg = {'MRAC with Scalar Filter'};
    leg_fill = cell(1, 2*(length(time)-1));
    leg_fill(:) = {''};
    legend(['System State', leg_fill, 'Observer State', leg_fill])
end
disp("Finished!")
