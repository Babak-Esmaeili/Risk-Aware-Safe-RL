function [x, u, P_1_inv, P_2_inv, P_3_inv, P_LQR_DD_2_opt_max, cost_average, t] = DD_optimal_controller(N,U0,x0_id,x0_sim,lambda,delta,sigma,noise)

    %% System matrices
    A = [ 0.2895 -0.0001
         -1.6012  0.0295 ];
    
    % A = [4/5 1/2;-2/5 6/5];
    % 
    % A = [ 0.2895 -0.0001
    %         7.2768  0.1519 ];
    % A = [-0.5975 -0.0127
    %         7.2768  0.1519 ];
    % A = [-0.5975  -0.0127
    %         96.0568  1.4139 ];
    
    B = [ 0 
          1 ];
    
    nStates = size(A,1);
    nInputs = size(B,2);
    
    %% Safe set
    F = [ 1/3   1/4
          0     1/4
         -4/12 -1/12
         -1/3  -1/4
          0    -1/4
          4/12  1/12 ];
    
    g = [ 1  1  1  1  1  1 ]';
    
    F_vertices = [ 3  0 -4 -3  0  4
                   0  4  4  0 -4 -4 ];
    
    nF = size(F,1);
    nVertices = nF;
    
    %% Noise realizations
%     nRealization = 100;
%     nSamples = 100000;
%     noise = zeros(nSamples,nRealization*nStates);
%     w = zeros(nRealization,nSamples);
%     
%     varNoise = 0.005;
%     coVarNoise = 0;
%     coVarMatrix = [  varNoise   coVarNoise
%         coVarNoise   varNoise  ];
%     sigma = coVarMatrix;
%     mu = zeros(1,nStates);
%     R = chol(sigma);
%     
%     for ii = 1:nRealization
%     
%         noise(:,nStates*ii-(nStates-1):nStates*ii) = repmat(mu,nSamples,1) + randn(nSamples,nStates)*R;
%     
%     end
    
    nRealization = size(noise,2)/nStates;
    
    %% Controller parameters
%     lambda = 0.8;
    
    nEllips = nVertices/2;
    
%     delta = 0.1;
    
    delta_n = nStates + 2*sqrt(nStates*log(1/delta)) + 2*log(1/delta);
    % delta_n = 1;
    
    %% Data collection
%     N = nStates + 10;
    
%     x0 = randn(nStates,1);
    
    X = zeros(nStates,N+1);
    X(:,1) = x0_id;
    
%     U0 = zeros(nInputs,N);
%     % U0(:,1:N) = 5*randn(nInputs,N);
%     
%     for i = 1:N
%     
%         rand_signal = rand(nInputs,1);
%     
%         if(round(rand_signal(1))==1)
%     
%             U0(1,i) = 1;
%     
%         else
%     
%             U0(1,i) = -1;
%     
%         end
%     
%     end
    
    W0 = noise(1:N,1:nStates)';
    
    for kk = 1:N
    
    %     X(:,kk+1) = A*X(:,kk) + B*U0(:,kk);
        X(:,kk+1) = A*X(:,kk) + B*U0(:,kk) + W0(:,kk);
    
    end
    
    X0 = X(:,1:N);
    X1 = X(:,2:N+1);
    
    % data = [U0; X0];
    data = X0;
    
    % Check rank
    if(rank(data)<min(size(data)))
    
        errordlg('Data are not full row rank!');
    
    end
    
    %% Optimization
    d_1 = F_vertices(:,1);
    % d_1 = d_1/norm(d_1);
    
    d_2 = F_vertices(:,2);
    % d_2 = d_2/norm(d_2);
    
    d_3 = F_vertices(:,3);
    % d_3 = d_3/norm(d_3);
    
    P_1 = sdpvar(nStates,nStates);
    P_2 = sdpvar(nStates,nStates);
    P_3 = sdpvar(nStates,nStates);
    
    Z_x_1 = sdpvar(nF,nF);
    Z_x_2 = sdpvar(nF,nF);
    Z_x_3 = sdpvar(nF,nF);
    
    mu_1 = sdpvar(1,1);
    mu_2 = sdpvar(1,1);
    mu_3 = sdpvar(1,1);
    
    Y_1 = sdpvar(N,nStates);
    Y_2 = sdpvar(N,nStates);
    Y_3 = sdpvar(N,nStates);
    
    H_1 = sdpvar(N,N);
    H_2 = sdpvar(N,N);
    H_3 = sdpvar(N,N);
    
    eta_1 = sdpvar(1,1);
    eta_2 = sdpvar(1,1);
    eta_3 = sdpvar(1,1);
    
    alpha_1 = sdpvar(1,1);
    alpha_2 = sdpvar(1,1);
    alpha_3 = sdpvar(1,1);
    
    tau_val = 0.01;
    tau_1 = tau_val;
    tau_2 = tau_val;
    tau_3 = tau_val;
    
    c1 = [         P_1                   X1*Y_3            eta_3*(sigma.^(1/2))
                 Y_3'*X1'         (lambda-tau_3)*P_3         zeros(nStates)
           eta_3'*(sigma.^(1/2))    zeros(nStates)    (tau_3/delta_n)*eye(nStates)  ] >= 0;
    c2 = [         P_2                   X1*Y_1            eta_1*(sigma.^(1/2))
                 Y_1'*X1'         (lambda-tau_1)*P_1         zeros(nStates)
           eta_1'*(sigma.^(1/2))    zeros(nStates)    (tau_1/delta_n)*eye(nStates)  ] >= 0;
    c3 = [         P_3                   X1*Y_2            eta_2*(sigma.^(1/2))
                 Y_2'*X1'         (lambda-tau_2)*P_2         zeros(nStates)
           eta_2'*(sigma.^(1/2))    zeros(nStates)    (tau_2/delta_n)*eye(nStates)  ] >= 0;
    % c1 = [    P_1       X1*Y_3
    %        (X1*Y_3)'  lambda*P_3   ] >= 0;
    % c2 = [    P_2       X1*Y_1
    %        (X1*Y_1)'  lambda*P_1   ] >= 0;
    % c3 = [    P_3       X1*Y_2
    %        (X1*Y_2)'  lambda*P_2   ] >= 0;
    
    c4 = P_1 >= 0.01*eye(2);
    c5 = P_2 >= 0.01*eye(2);
    c6 = P_3 >= 0.01*eye(2);
    % c4 = P_1 >= 0;
    % c5 = P_2 >= 0;
    % c6 = P_3 >= 0;
    
    c7 = [  Z_x_1    F*P_1
            (F*P_1)'   P_1  ] >= 0;
    c8 = [  Z_x_2    F*P_2
            (F*P_2)'   P_2  ] >= 0;
    c9 = [  Z_x_3    F*P_3
            (F*P_3)'   P_3  ] >= 0;
    
    c10 = Z_x_1(1,1) <= g(1)^2;
    c11 = Z_x_1(2,2) <= g(2)^2;
    c12 = Z_x_1(3,3) <= g(3)^2;
    c13 = Z_x_1(4,4) <= g(4)^2;
    c14 = Z_x_1(5,5) <= g(5)^2;
    c15 = Z_x_1(6,6) <= g(6)^2;
    
    c16 = Z_x_2(1,1) <= g(1)^2;
    c17 = Z_x_2(2,2) <= g(2)^2;
    c18 = Z_x_2(3,3) <= g(3)^2;
    c19 = Z_x_2(4,4) <= g(4)^2;
    c20 = Z_x_2(5,5) <= g(5)^2;
    c21 = Z_x_2(6,6) <= g(6)^2;
    
    c22 = Z_x_3(1,1) <= g(1)^2;
    c23 = Z_x_3(2,2) <= g(2)^2;
    c24 = Z_x_3(3,3) <= g(3)^2;
    c25 = Z_x_3(4,4) <= g(4)^2;
    c26 = Z_x_3(5,5) <= g(5)^2;
    c27 = Z_x_3(6,6) <= g(6)^2;
    
    c28 = Z_x_1 >= 0;
    c29 = Z_x_2 >= 0;
    c30 = Z_x_3 >= 0;
    
    c31 = [    1      mu_1*d_1'
            mu_1*d_1     P_1    ] >= 0;
    c32 = [    1      mu_2*d_2'
            mu_2*d_2     P_2    ] >= 0;
    c33 = [    1      mu_3*d_3'
            mu_3*d_3     P_3    ] >= 0;
    
    c34 = mu_1 >= 0;
    c35 = mu_2 >= 0;
    c36 = mu_3 >= 0;
    
    c37 = X0*Y_1 == P_1;
    c38 = X0*Y_2 == P_2;
    c39 = X0*Y_3 == P_3;
    
    c40 = [ H_1   Y_1
            Y_1'  P_1 ] >= 0;
    c41 = [ H_2   Y_2
            Y_2'  P_2 ] >= 0;
    c42 = [ H_3   Y_3
            Y_3'  P_3 ] >= 0;
    
    c43 = [ alpha_1+1  eta_1
              eta_1      1   ] >= 0;
    c44 = [ alpha_2+1  eta_2
              eta_2      1   ] >= 0;
    c45 = [ alpha_3+1  eta_3
              eta_3      1   ] >= 0;
    
    c46 = trace(H_1) <= alpha_1;
    c47 = trace(H_2) <= alpha_2;
    c48 = trace(H_3) <= alpha_3;
    
    c49  = eta_1 >= 0;
    c50  = eta_2 >= 0;
    c51  = eta_3 >= 0;
    
    c52  = alpha_1 >= 0;
    c53  = alpha_2 >= 0;
    c54  = alpha_3 >= 0;
    
    % c55 = Y_1 >= 0;
    % c56 = Y_2 >= 0;
    % c57 = Y_3 >= 0;
    % 
    % c58 = H_1 >= 0;
    % c59 = H_2 >= 0;
    % c60 = H_3 >= 0;
    
    constraints = c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10 + c11 + c12 + ...
                  c13 + c14 + c15 + c16 + c17 + c18 + c19 + c20 + c21 + c22 + ...
                  c23 + c24 + c25 + c26 + c27 + c28 + c29 + c30 + c31 + c32 + c33 + ...
                  c34 + c35 + c36 + c37 + c38 + c39 + c40 + c41 + c42 + c43 + c44 + c45 + ...
                  c46 + c47 + c48 + c49 + c50 + c51 + c52 + c53 + c54;
    %  + c55 + c56 + c57 + c58 + c59 + c60
    
    cost = -1000*(mu_1+mu_2+mu_3) + (eta_1+eta_2+eta_3) + (alpha_1+alpha_2+alpha_3);
    
    ops = sdpsettings('solver','mosek','verbose',0);
    diagnostics = optimize(constraints, cost, ops);
    
    if(diagnostics.problem == 0)
        
        P_1_opt = value(P_1);
        P_2_opt = value(P_2);
        P_3_opt = value(P_3);
    
        Y_1_opt = value(Y_1);
        Y_2_opt = value(Y_2);
        Y_3_opt = value(Y_3);
    
    else
    
        fprintf('Error! Optimization problem is infeasible!\n\n');
    
    end
    
    P_1_inv = P_1_opt\[1 0;0 1];
    P_2_inv = P_2_opt\[1 0;0 1];
    P_3_inv = P_3_opt\[1 0;0 1];
    
    G_K_1_opt = Y_1_opt*P_1_inv;
    G_K_2_opt = Y_2_opt*P_2_inv;
    G_K_3_opt = Y_3_opt*P_3_inv;
    
    K_1_opt = U0*G_K_1_opt;
    K_2_opt = U0*G_K_2_opt;
    K_3_opt = U0*G_K_3_opt;
    
    %% Set partitioning
    P_all = {P_1_opt,P_2_opt,P_3_opt};
    v_h_all = [];
    K_opt_all = {K_1_opt,K_2_opt,K_3_opt};
    
    syms beta_1 beta_2 real;
    beta = [beta_1 beta_2]';
    
    for j_1 = 1:nEllips-1
    
        for j_2 = j_1+1:nEllips
    
            eq1 = beta'*P_all{j_1}*beta == 1;
            eq2 = beta'*P_all{j_2}*beta == 1;
    
            sol = solve([eq1,eq2]);
    
            beta_sol = [double(sol.beta_1) double(sol.beta_2)]';
    
            v_j1 = P_all{j_1}*beta_sol;
            v_j2 = P_all{j_2}*beta_sol;
    
            v_h = [v_j1 v_j2];
            v_h_all = [v_h_all v_h];    %#ok
    
        end
    
    end
    
    indices = convhull(v_h_all(1,:),v_h_all(2,:));
    % indices = indices(1:end-1);
    v_h_all_accepted = v_h_all(:,indices);
    
    nPartitions = size(indices,1) - 1;
    
    %% Compute K for all partitions
    label_v_h_all_accepted = zeros(1,size(v_h_all_accepted,2));
    
    for ii = 1:size(indices,1)
    
        if((indices(ii)-rem(indices(ii),2*(nStates^2)+1))/(2*(nStates^2)+1)==0) % First group (1:8)
    
            j_1 = 1;
            j_2 = 2;
    
            if((indices(ii)-rem(indices(ii),1+(nStates^2)))/(1+(nStates^2))==0)
    
                label_v_h_all_accepted(ii) = j_1;
    
            else
    
                label_v_h_all_accepted(ii) = j_2;
    
            end
    
        elseif((indices(ii)-rem(indices(ii),2*(nStates^2)+1))/(2*(nStates^2)+1)==1) % Second group (9:16)
    
            j_1 = 1;
            j_2 = 3;
    
            if((indices(ii)-rem(indices(ii),1+2*(nStates^2)+(nStates^2)))/(1+2*(nStates^2)+(nStates^2))==0)
    
                label_v_h_all_accepted(ii) = j_1;
    
            else
    
                label_v_h_all_accepted(ii) = j_2;
    
            end
    
        else    % Third group (17:24) 
    
            j_1 = 2;
            j_2 = 3;
    
            if((indices(ii)-rem(indices(ii),1+2*(nStates^2)+2*(nStates^2)+(nStates^2)))/(1+2*(nStates^2)+2*(nStates^2)+(nStates^2))==0)
    
                label_v_h_all_accepted(ii) = j_1;
    
            else
    
                label_v_h_all_accepted(ii) = j_2;
    
            end
    
        end
    
    end
    
    K_opt_all_partiotions = zeros(nPartitions,nStates);
    
    for iii = 1:nPartitions
    
        if(label_v_h_all_accepted(iii) == label_v_h_all_accepted(iii+1))
    
            K_opt_all_partiotions(iii,:) = K_opt_all{label_v_h_all_accepted(iii)};
    
        else
    
            K_opt_all_partiotions(iii,:) = [K_opt_all{label_v_h_all_accepted(iii)}*v_h_all_accepted(:,iii) K_opt_all{label_v_h_all_accepted(iii+1)}*v_h_all_accepted(:,iii+1)]*([v_h_all_accepted(:,iii) v_h_all_accepted(:,iii+1)]^-1);
    
        end
    
    end
    
    %% F of convex hull
    F_CH = find_F(v_h_all_accepted(:,1:end-1));
    g_CH = ones(size(v_h_all_accepted(:,1:end-1),2),1);
    
    nF_CH = size(F_CH,1);
    
    %% Model-based LQR design
    C = [ 10  0
          0   1 ];
    D = 0;
    
    Q = C'*C;
    R = D'*D;
    
    [K_LQR,P_LQR] = dlqr(A,B,Q,R);
    
    Z_x_LQR = sdpvar(nF_CH,nF_CH);
    zeta = sdpvar(1,1);
    
    c1_LQR = [         Z_x_LQR          F_CH*(P_LQR^-1*zeta)
               (F_CH*(P_LQR^-1*zeta))'     P_LQR^-1*zeta     ] >= 0;
    
    c2_LQR = Z_x_LQR(1,1) <= g_CH(1)^2;
    c3_LQR = Z_x_LQR(2,2) <= g_CH(2)^2;
    c4_LQR = Z_x_LQR(3,3) <= g_CH(3)^2;
    c5_LQR = Z_x_LQR(4,4) <= g_CH(4)^2;
    c6_LQR = Z_x_LQR(5,5) <= g_CH(5)^2;
    c7_LQR = Z_x_LQR(6,6) <= g_CH(6)^2;
    c8_LQR = Z_x_LQR(7,7) <= g_CH(7)^2;
    c9_LQR = Z_x_LQR(8,8) <= g_CH(8)^2;
    c10_LQR = Z_x_LQR(9,9) <= g_CH(9)^2;
    c11_LQR = Z_x_LQR(10,10) <= g_CH(10)^2;
    c12_LQR = Z_x_LQR(11,11) <= g_CH(11)^2;
    c13_LQR = Z_x_LQR(12,12) <= g_CH(12)^2;
    
    c14_LQR = zeta >= 0;
    
    constraints_LQR = c1_LQR + c2_LQR + c3_LQR + c4_LQR + c5_LQR + c6_LQR + c7_LQR + ...
                      c8_LQR + c9_LQR + c10_LQR + c11_LQR + c12_LQR + c13_LQR + c14_LQR;
    
    cost_LQR = -zeta;
    
    ops = sdpsettings('solver','mosek','verbose',0);
    diagnostics = optimize(constraints_LQR, cost_LQR, ops);
    
    if(diagnostics.problem == 0)
        
        zeta_opt = value(zeta);
    
    else
    
        fprintf('Error! Optimization problem is infeasible!\n\n');
    
    end
    
    P_LQR_opt = P_LQR/zeta_opt;
    
    %% Data-based LQR design #1
    Z = C*X0 + D*U0;
    
    Y_LQR_DD = sdpvar(nStates,nStates);
    THETA = sdpvar(N,nStates);
    
    c1_LQR_DD = [  Y_LQR_DD      THETA'*X1'       THETA'*Z'
                  X1*THETA       Y_LQR_DD       zeros(nStates)
                   Z*THETA  zeros(nStates)  eye(nStates)   ] >= eps;
    
    c2_LQR_DD = Y_LQR_DD == X0*THETA;
    
    constraints_LQR_DD = c1_LQR_DD + c2_LQR_DD;
    
    cost_LQR_DD = [];
    
    ops = sdpsettings('solver','mosek','verbose',0);
    diagnostics = optimize(constraints_LQR_DD, cost_LQR_DD, ops);
    
    if(diagnostics.problem == 0)
        
        Y_LQR_DD_opt = value(Y_LQR_DD);
    
        THETA_opt = value(THETA);
    
    else
    
        fprintf('Error! Optimization problem is infeasible!\n\n');
    
    end
    
    P_LQR_DD_opt = Y_LQR_DD_opt\eye(nStates);
    
    K_LQR_DD = U0*THETA_opt*P_LQR_DD_opt;
    
    % Maximum level-set of the optimal set (P)
    Z_x_LQR_DD_P = sdpvar(nF_CH,nF_CH);
    zeta_DD_P = sdpvar(1,1);
    
    c1_LQR_DD_P = [              Z_x_LQR_DD_P              F_CH*((P_LQR_DD_opt^-1)*zeta_DD_P)
                    (F_CH*((P_LQR_DD_opt^-1)*zeta_DD_P))'      (P_LQR_DD_opt^-1)*zeta_DD_P    ] >= 0;
    
    c2_LQR_DD_P = Z_x_LQR_DD_P(1,1) <= g_CH(1)^2;
    c3_LQR_DD_P = Z_x_LQR_DD_P(2,2) <= g_CH(2)^2;
    c4_LQR_DD_P = Z_x_LQR_DD_P(3,3) <= g_CH(3)^2;
    c5_LQR_DD_P = Z_x_LQR_DD_P(4,4) <= g_CH(4)^2;
    c6_LQR_DD_P = Z_x_LQR_DD_P(5,5) <= g_CH(5)^2;
    c7_LQR_DD_P = Z_x_LQR_DD_P(6,6) <= g_CH(6)^2;
    c8_LQR_DD_P = Z_x_LQR_DD_P(7,7) <= g_CH(7)^2;
    c9_LQR_DD_P = Z_x_LQR_DD_P(8,8) <= g_CH(8)^2;
    c10_LQR_DD_P = Z_x_LQR_DD_P(9,9) <= g_CH(9)^2;
    c11_LQR_DD_P = Z_x_LQR_DD_P(10,10) <= g_CH(10)^2;
    c12_LQR_DD_P = Z_x_LQR_DD_P(11,11) <= g_CH(11)^2;
    c13_LQR_DD_P = Z_x_LQR_DD_P(12,12) <= g_CH(12)^2;
    
    c14_LQR_DD_P = zeta_DD_P >= 0;
    
    constraints_LQR_DD_P = c1_LQR_DD_P + c2_LQR_DD_P + c3_LQR_DD_P + c4_LQR_DD_P + c5_LQR_DD_P + ...
                           c6_LQR_DD_P + c7_LQR_DD_P + c8_LQR_DD_P + c9_LQR_DD_P + c10_LQR_DD_P + ...
                           c11_LQR_DD_P + c12_LQR_DD_P + c13_LQR_DD_P + c14_LQR_DD_P;
    
    cost_LQR_DD_P = -zeta_DD_P;
    
    ops = sdpsettings('solver','mosek','verbose',0);
    diagnostics = optimize(constraints_LQR_DD_P, cost_LQR_DD_P, ops);
    
    if(diagnostics.problem == 0)
        
        zeta_DD_P_opt = value(zeta_DD_P);
    
    else
    
        fprintf('Error! Optimization problem is infeasible!\n\n');
    
    end
    
    P_LQR_DD_opt_max = P_LQR_DD_opt/zeta_DD_P_opt;
    
    %% Data-based LQR design #2
    QQ_LQR_DD_2 = [ 10000  0
                    0      0.01 ];
    RR_LQR_DD_2 = 50;

    alpha_LQR_DD_2 = 1;
    
    gamma_LQR_DD_2 = sdpvar(1,1);
    P_LQR_DD_2 = sdpvar(nStates,nStates);
    Q_LQR_DD_2 = sdpvar(N,nStates);
    L_LQR_DD_2 = sdpvar(nInputs,nInputs);
    V_LQR_DD_2 = sdpvar(N,N);
    
    c1_LQR_DD_2 = [ P_LQR_DD_2-eye(nStates)  X1*Q_LQR_DD_2
                        (X1*Q_LQR_DD_2)'       P_LQR_DD_2  ] >= 0;
    c2_LQR_DD_2 = P_LQR_DD_2 >= eye(nStates);
    c3_LQR_DD_2 = [    L_LQR_DD_2    U0*Q_LQR_DD_2
                    (U0*Q_LQR_DD_2)'   P_LQR_DD_2  ] >= 0;
    c4_LQR_DD_2 = [ V_LQR_DD_2   Q_LQR_DD_2
                    Q_LQR_DD_2'  P_LQR_DD_2 ] >= 0;
    c5_LQR_DD_2 = X0*Q_LQR_DD_2 == P_LQR_DD_2;
    c6_LQR_DD_2 = trace(QQ_LQR_DD_2*P_LQR_DD_2) + trace(RR_LQR_DD_2*L_LQR_DD_2) + alpha_LQR_DD_2*trace(V_LQR_DD_2) <= gamma_LQR_DD_2;
    
    constraints_LQR_DD_2 = c1_LQR_DD_2 + c2_LQR_DD_2 + c3_LQR_DD_2 + c4_LQR_DD_2 + c5_LQR_DD_2 + c6_LQR_DD_2;
    
    cost_LQR_DD_2 = gamma_LQR_DD_2;
    
    ops = sdpsettings('solver','mosek','verbose',0);
    diagnostics = optimize(constraints_LQR_DD_2, cost_LQR_DD_2, ops);
    
    if(diagnostics.problem == 0)
        
        P_LQR_DD_2_opt = value(P_LQR_DD_2);
        Q_LQR_DD_2_opt = value(Q_LQR_DD_2);
    
    else
    
        fprintf('Error! Optimization problem is infeasible!\n\n');
    
    end
    
    K_LQR_DD_2 = U0*Q_LQR_DD_2_opt*(P_LQR_DD_2_opt\eye(nStates));
    
    % Maximum level-set of the optimal set (P)
    Z_x_LQR_DD_2_P = sdpvar(nF_CH,nF_CH);
    zeta_DD_2_P = sdpvar(1,1);
    
    c1_LQR_DD_2_P = [              Z_x_LQR_DD_2_P                F_CH*((P_LQR_DD_2_opt^-1)*zeta_DD_2_P)
                      (F_CH*((P_LQR_DD_2_opt^-1)*zeta_DD_2_P))'      (P_LQR_DD_2_opt^-1)*zeta_DD_2_P    ] >= 0;
    
    c2_LQR_DD_2_P = Z_x_LQR_DD_2_P(1,1) <= g_CH(1)^2;
    c3_LQR_DD_2_P = Z_x_LQR_DD_2_P(2,2) <= g_CH(2)^2;
    c4_LQR_DD_2_P = Z_x_LQR_DD_2_P(3,3) <= g_CH(3)^2;
    c5_LQR_DD_2_P = Z_x_LQR_DD_2_P(4,4) <= g_CH(4)^2;
    c6_LQR_DD_2_P = Z_x_LQR_DD_2_P(5,5) <= g_CH(5)^2;
    c7_LQR_DD_2_P = Z_x_LQR_DD_2_P(6,6) <= g_CH(6)^2;
    c8_LQR_DD_2_P = Z_x_LQR_DD_2_P(7,7) <= g_CH(7)^2;
    c9_LQR_DD_2_P = Z_x_LQR_DD_2_P(8,8) <= g_CH(8)^2;
    c10_LQR_DD_2_P = Z_x_LQR_DD_2_P(9,9) <= g_CH(9)^2;
    c11_LQR_DD_2_P = Z_x_LQR_DD_2_P(10,10) <= g_CH(10)^2;
    c12_LQR_DD_2_P = Z_x_LQR_DD_2_P(11,11) <= g_CH(11)^2;
    c13_LQR_DD_2_P = Z_x_LQR_DD_2_P(12,12) <= g_CH(12)^2;
    
    c14_LQR_DD_2_P = zeta_DD_2_P >= 0;
    
    constraints_LQR_DD_2_P = c1_LQR_DD_2_P + c2_LQR_DD_2_P + c3_LQR_DD_2_P + c4_LQR_DD_2_P + c5_LQR_DD_2_P + ...
                             c6_LQR_DD_2_P + c7_LQR_DD_2_P + c8_LQR_DD_2_P + c9_LQR_DD_2_P + c10_LQR_DD_2_P + ...
                             c11_LQR_DD_2_P + c12_LQR_DD_2_P + c13_LQR_DD_2_P + c14_LQR_DD_2_P;
    
    cost_LQR_DD_2_P = -zeta_DD_2_P;
    
    ops = sdpsettings('solver','mosek','verbose',0);
    diagnostics = optimize(constraints_LQR_DD_2_P, cost_LQR_DD_2_P, ops);
    
    if(diagnostics.problem == 0)
        
        zeta_DD_2_P_opt = value(zeta_DD_2_P);
    
    else
    
        fprintf('Error! Optimization problem is infeasible!\n\n');
    
    end
    
    P_LQR_DD_2_opt_max = P_LQR_DD_2_opt/zeta_DD_2_P_opt;
    
    %% Time sequence
    T0 = 0;
    Ts = 1;
    Tf = 20;
    
    t = T0:Ts:Tf;
    
    nSteps = numel(t);
    
    %% Main loop
    % x0 = [-2,4]';
    % x0 = [-3.5,2]';
    % x0 = [1.5,2]';
    % x0 = [0,-4]';
    % x0 = [-3,4]';
    % x0 = [0,3.7]';    % For LQR
%     x0 = [-3.5,2.0]';
    
    x = zeros(nRealization,nSteps+1);
    
    norm_K = zeros(nRealization,nSteps);
    
    a = zeros(1,nSteps);
    b = zeros(1,nSteps);
    
    nu = zeros(1,nSteps);

    cost_total = zeros(nRealization,nSteps+1);
    
    noise_flag = 1;
    
    for ii = 1:nRealization
    
        % Simulation
        x(nStates*ii-(nStates-1):nStates*ii,1) = x0_sim;
    
        u = zeros(nInputs,nSteps);
    
        w = noise(1:nSteps,nStates*ii-(nStates-1):nStates*ii)';
    
        for k = 1:nSteps
    
            [~, idx] = dist_calc(x(nStates*ii-(nStates-1):nStates*ii,k),v_h_all_accepted(:,1:end-1));
    
            K_safe = K_opt_all_partiotions(idx,:);
    
            if(x(nStates*ii-(nStates-1):nStates*ii,k)'*P_LQR_DD_2_opt_max*x(nStates*ii-(nStates-1):nStates*ii,k)<=1)
    
                K_interpolation = K_LQR_DD_2;
    
            else
    
                a(k) = sqrt(x(nStates*ii-(nStates-1):nStates*ii,k)'*(P_LQR_DD_2_opt_max\eye(nStates))*x(nStates*ii-(nStates-1):nStates*ii,k));
                b(k) = find_distance_interpolation(x(nStates*ii-(nStates-1):nStates*ii,k),v_h_all_accepted,idx);
    
                nu(k) = b(k)*(a(k)-1)/(a(k)-b(k));
                
%                 K_interpolation = nu(k)*K_safe/b(k) + (1-nu(k))*K_LQR_DD_2/a(k);
                K_interpolation = K_LQR_DD_2;
    
            end
            
            u(:,k) = K_interpolation*x(nStates*ii-(nStates-1):nStates*ii,k);
    
            x(nStates*ii-(nStates-1):nStates*ii,k+1) = A*x(nStates*ii-(nStates-1):nStates*ii,k) + B*u(:,k) + noise_flag*w(:,k);
    
            norm_K(ii,k) = norm(K_interpolation);

            cost_total(ii,k+1) = cost_total(ii,k) + Ts*(x(nStates*ii-(nStates-1):nStates*ii,k)'*QQ_LQR_DD_2*x(nStates*ii-(nStates-1):nStates*ii,k) + u(:,k)'*RR_LQR_DD_2*u(:,k));
    
        end
    
    end
    
    x = x(:,1:end-1);

    cost_average = sum(cost_total(:,1:end-1),1)/nRealization;

end
