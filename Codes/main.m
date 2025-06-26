clc;
clear;
close all;

tic;
%% Noise generation
nStates = 2;
nRealization = 100;
nSample = 100000;

mu = [0 0];
varNoise = 0.01;
coVarNoise = 0;
coVarMatrix = [  varNoise   coVarNoise
                coVarNoise   varNoise  ];
sigma = coVarMatrix;

noise = noise_generator(mu,sigma,nStates,nSample,nRealization);

%% Control parameters
lambda = 0.8;

delta = 0.1;

eta = 0.99;

%% Data parameters
N = nStates + 10;

%% U0 data
nInputs = 1;
U0 = zeros(nInputs,N);
% U0(:,1:N) = 5*randn(nInputs,N);

for i = 1:N

    rand_signal = rand(nInputs,1);

    if(round(rand_signal(1))==1)

        U0(1,i) = 1;

    else

        U0(1,i) = -1;

    end

end

%% Initial condition for identification
x0_id = randn(nStates,1);

%% Initial condition for simulation
% x0_sim = [-2,4]';
% x0_sim = [-3.5,2]';
% x0_sim = [1.5,2]';
% x0_sim = [0,-4]';
% x0_sim = [-3,4]';
% x0_sim = [0,3.7]';    % For LQR
x0_sim = [-3.5,2.0]';

%% Simulation
[x_o, u_o, P_1_inv, P_2_inv, P_3_inv, P_LQR_DD_2_opt_max, cost_average_o, t] = DD_optimal_controller(N,U0,x0_id,x0_sim,lambda,delta,sigma,noise);
[x_s, u_s, ~, ~, ~, ~, cost_average_s, ~] = DD_safe_controller(N,U0,x0_id,x0_sim,lambda,delta,sigma,noise);
% [x_so, u_so, ~, ~, ~, ~, cost_average_so, ~] = DD_safe_optimal_controller(N,U0,x0_id,x0_sim,lambda,delta,sigma,noise);
% [x_so, u_so, ~, ~, ~, ~, cost_average_so, ~] = DD_safe_optimal_controller_2(N,U0,x0_id,x0_sim,lambda,delta,eta,sigma,noise);
[x_so, u_so, ~, ~, ~, ~, cost_average_so, ~] = DD_safe_optimal_controller_3(N,U0,x0_id,x0_sim,lambda,delta,eta,sigma,noise);

x_o_end = [ sum(x_o(1:2:end,end))/(nRealization/2)
            sum(x_o(2:2:end,end))/(nRealization/2) ];

x_s_end = [ sum(x_s(1:2:end,end))/(nRealization/2)
            sum(x_s(2:2:end,end))/(nRealization/2) ];

x_so_end = [ sum(x_so(1:2:end,end))/(nRealization/2)
             sum(x_so(2:2:end,end))/(nRealization/2) ];

toc;
%% Plot results
xx = sdpvar(nStates,1);

%% Phase plane
figure(1);
subplot(3,1,1);
plot(t,cost_average_o,'b-','LineWidth',1);
grid on;
xlabel('Time [s]','FontSize',20);
ylabel('Cost','FontSize',20);
title('Average cost value for different noise realizations');

subplot(3,1,2);
plot(t,cost_average_s,'b-','LineWidth',1);
grid on;
xlabel('Time [s]','FontSize',20);
ylabel('Cost','FontSize',20);

subplot(3,1,3);
plot(t,cost_average_so,'b-','LineWidth',1);
grid on;
xlabel('Time [s]','FontSize',20);
ylabel('Cost','FontSize',20);

figure(2);
subplot(3,1,1);
plot(t,u_o,'b-','LineWidth',1);
grid on;
xlabel('Time [s]','FontSize',20);
ylabel('u(k)','FontSize',20);
title('Control input for one of the noise realizations');

subplot(3,1,2);
plot(t,u_s,'b-','LineWidth',1);
grid on;
xlabel('Time [s]','FontSize',20);
ylabel('u(k)','FontSize',20);
title('Control input for one of the noise realizations');

subplot(3,1,3);
plot(t,u_so,'b-','LineWidth',1);
grid on;
xlabel('Time [s]','FontSize',20);
ylabel('u(k)','FontSize',20);
title('Control input for one of the noise realizations');

figure(3);
subplot(1,3,1);
F_1 = xx'*P_1_inv*xx <= 1;
F_2 = xx'*P_2_inv*xx <= 1;
F_3 = xx'*P_3_inv*xx <= 1;
H = hull(F_1,F_2,F_3);
plot(H,[],'w');
alpha(0);
hold on;
grid on;

nLS = 8;
for k = 1:nLS

    F_1 = xx'*(lambda^(-k))*P_1_inv*xx <= 1;
    F_2 = xx'*(lambda^(-k))*P_2_inv*xx <= 1;
    F_3 = xx'*(lambda^(-k))*P_3_inv*xx <= 1;
    H = hull(F_1,F_2,F_3);
    plot(H,[],'w');
    hold on;

end
alpha(0);
plot(xx'*P_1_inv*xx<=1,[],'r');
plot(xx'*P_2_inv*xx<=1,[],'g');
plot(xx'*P_3_inv*xx<=1,[],'c');
alpha(0.5);
plot(xx'*P_LQR_DD_2_opt_max*xx<=1,[],'k');
alpha(0.3);
plot(x0_sim(1),x0_sim(2),'sm','MarkerSize',15,'MarkerFaceColor','m');
hold on;
for ii = 1:nRealization

    plot(x_o(2*ii-1,:),x_o(2*ii,:),'b-','LineWidth',1,'MarkerSize',12);

end
grid on;
% axis equal;
plot(x_o_end(1),x_o_end(2),'pm','MarkerSize',15,'MarkerFaceColor','m');
xlabel('x_1 (t)','FontSize',20);
ylabel('x_2 (t)','FontSize',20);
% title('Optimal Controller');
legend('','','','','','','','','','','','','','Initial point','FontSize',20,'TextColor','k');

subplot(1,3,2);
F_1 = xx'*P_1_inv*xx <= 1;
F_2 = xx'*P_2_inv*xx <= 1;
F_3 = xx'*P_3_inv*xx <= 1;
H = hull(F_1,F_2,F_3);
plot(H,[],'w');
alpha(0);
hold on;
grid on;

nLS = 8;
for k = 1:nLS

    F_1 = xx'*(lambda^(-k))*P_1_inv*xx <= 1;
    F_2 = xx'*(lambda^(-k))*P_2_inv*xx <= 1;
    F_3 = xx'*(lambda^(-k))*P_3_inv*xx <= 1;
    H = hull(F_1,F_2,F_3);
    plot(H,[],'w');
    hold on;

end
alpha(0);
plot(xx'*P_1_inv*xx<=1,[],'r');
plot(xx'*P_2_inv*xx<=1,[],'g');
plot(xx'*P_3_inv*xx<=1,[],'c');
alpha(0.5);
plot(x0_sim(1),x0_sim(2),'sm','MarkerSize',15,'MarkerFaceColor','m');
hold on;
for ii = 1:nRealization

    plot(x_s(2*ii-1,:),x_s(2*ii,:),'b-','LineWidth',1,'MarkerSize',12);

end
grid on;
% axis equal;
plot(x_s_end(1),x_s_end(2),'pm','MarkerSize',15,'MarkerFaceColor','m');
xlabel('x_1 (t)','FontSize',20);
ylabel('x_2 (t)','FontSize',20);
% title('Minimum Variance-based Safe Controller');
legend('','','','','','','','','','','','','Initial point','FontSize',20,'TextColor','k');

subplot(1,3,3);
F_1 = xx'*P_1_inv*xx <= 1;
F_2 = xx'*P_2_inv*xx <= 1;
F_3 = xx'*P_3_inv*xx <= 1;
H = hull(F_1,F_2,F_3);
plot(H,[],'w');
alpha(0);
hold on;
grid on;

nLS = 8;
for k = 1:nLS

    F_1 = xx'*(lambda^(-k))*P_1_inv*xx <= 1;
    F_2 = xx'*(lambda^(-k))*P_2_inv*xx <= 1;
    F_3 = xx'*(lambda^(-k))*P_3_inv*xx <= 1;
    H = hull(F_1,F_2,F_3);
    plot(H,[],'w');
    hold on;

end
alpha(0);
plot(xx'*P_1_inv*xx<=1,[],'r');
plot(xx'*P_2_inv*xx<=1,[],'g');
plot(xx'*P_3_inv*xx<=1,[],'c');
alpha(0.5);
plot(xx'*P_LQR_DD_2_opt_max*xx<=1,[],'k');
alpha(0.3);
plot(x0_sim(1),x0_sim(2),'sm','MarkerSize',15,'MarkerFaceColor','m');
hold on;
for ii = 1:nRealization

    plot(x_so(2*ii-1,:),x_so(2*ii,:),'b-','LineWidth',1,'MarkerSize',12);

end
grid on;
% axis equal;
plot(x_so_end(1),x_so_end(2),'pm','MarkerSize',15,'MarkerFaceColor','m');
xlabel('x_1 (t)','FontSize',20);
ylabel('x_2 (t)','FontSize',20);
% title('Minimum Variance-based Safe Optimal Controller');
legend('','','','','','','','','','','','','','Initial point','FontSize',20,'TextColor','k');
