
% Discrete Input Delay Solver. This function accepts the following parameters:

% 1. The agent model (A, B, C). Types: matrices.
% 2. The input delays (kappas). Type: array of (positive) integers.
% 3. The adjacency matrix of the graph (A_script). Type: matrix.
% 4. The set of leader agents (leader_set). Type: array.
% 5. The initial conditions (initial_conditions). Type: matrix.

% Note on initial conditions: Define these as a 2*n*N + n by K_max + 1 
% matrix for full-state, and 3*n*N + n by K_max + 1 for partial state. The 
% kth column of your matrix represent the stacked state ([x;chi;x_r] for
% full-state, [x;xhat;xhi;x_r] for partial state) evaluated at time k.
% Currently, there is a placeholder initial conditions defined in the 
% function, so any argument passed into the corresponding argument
% (e.g. 0) will be ignored in favor of the pre-defined initial conditions.
% Delete these if you wish to pass in your own.

% 6. Kmax: The maximum time step to evaluate the solution at, i.e.,
% you have the states for times k = 1, ..., Kmax.

% From this information, input_delay_solver implements the protocol 
% described in the Discrete Time Scale-free Input Delay paper. 
% Specifically, the discrete-time linear system is solved, 
% and synchronization of the agents is achieved. If the 
% given C is the identity, full-state coupling is implemented. Otherwise, 
% the protocol for  partial-state coupling is enacted. The algorithm
% is the obvious one (i.e., for k = 0, ..., Kmax-1, start with x(k), 
% update to x(k+1) according to the linear system).

% The following data is returned.

% 1. x (the state). Type: matrix. Remark: Specifically, an n*N by 
% Kmax matrix, where n is the dimension of the state, and N is the number 
% of agents.
% The columns of this matrix give the state values at the corresponding
% time steps. Furthermore, the individual agent states are stacked, i.e., 
% at a given time step, x = [x1; ...; xN]. You can get the individual 
% states at all times by writing x_1 = x(1:n,:), x_2 = (n+1:2*n,:), and so on.
% 2. x_r (the exosystem). Type: matrix (n by Kmax + kappa_bar + 1)
% 4. u (the input). Type: matrix (same structure as above).
%
% NOTE: Make sure to download the initial_conditions.m and
% protocol_design.m files for this function to compile in its current form,
% since it calls the functions defined in those files.
%
% Some things you can do after solving:
% 1. Let k = 1:K_max (time steps)
% 2. Plot the states: call plot(k,x).
% Note: whenever plotting any of this data, make sure to specify the time 
% steps k  (e.g., plot(k,x), NOT plot(x)). MATLAB will go awry without them.
% 3. Plot xtilde (which should -> 0). Note x = [x_1; x_2; ...; x_N], so
% xtilde = [x_1 - x_r; x_2 - x_r; ...; x_N - x_r]. Then, call
% plot(k,xtilde).

function [x x_r] = discrete_input_delay_solver(A,B,C, kappas, initial_conditions, ...
    A_script, leader_set, K_max)

    kappa_bar = max(kappas);
    
    % extract dimensions from data (i.e., N, n, q)
    sz1 = size(A);
    N = sz1(1);
    sz2 = size(A_script);
    n = sz2(1);
    sz3 = size(C);
    q = sz3(1);

    % initialize solution
    X = zeros(3*N*n + n, K_max);


    % get leader from leader_set
    leader = zeros(N,1);
    for i = 1:N
        if ismember(i, leader_set)
            leader(i) = 1;
        end
    end
    
    % get Laplacian, d_in
    L = zeros(N,N);
    d_in = zeros(N,1);
    for i = 1:N 
        v = 0;
        d_in(i) = 0;
        for j = 1:N
            v = v + A_script(i,j);
            d_in(i) = d_in(i) + A_script(i,j);
        end
        L(i,i) = v;
    end
    
    for i = 1:N
        for j = 1:N
            if i ~= j
                L(i,j) = -A_script(i,j);
            end
        end
    end
    % extended Laplacian
    L_bar = L + diag(leader);
    
    % get D_bar
    D_in = diag(d_in);
    D_bar = eye(N) - inv(2*eye(N) + D_in)*L_bar;
    
    % define helper vector
    d_hat = zeros(N,1);
    for i = 1:N
         d_hat(i) = 0;
        for j = 1:N
            d_hat(i) = d_hat(i) + D_bar(i,j);
        end
    end
    
    % initial conditions are given as an example here. Currently,
    % the code runs by putting any arbitrary argument (e.g. 0) in
    % the initial conditions argument of the solver, as it will use these.
    % When putting initial_conditions in for yourself, please 
    % format them as follows (i.e., as a matrix with the proper
    % dimensions)
    initial_conditions = zeros(3*N*n + n,kappa_bar + 1);
    ic = zeros(3*N*n + n,1);
    for j = 1:3*N*n+n
        ic(j) = j;
    end
    for i = 1:kappa_bar+1
        initial_conditions(:,i) = ic;
    end
    
    
    % design your protocol, using the provided function. Make
    % sure to download the corresponding m-file for this to compile
    % [epsilon, rho, taubar_max, P, K] = protocol_design(A,B,C,tau_bar);
    [epsilon, rho, kappabar_max, K, F] = discrete_protocol_design(A,B,C,kappa_bar);
    
    % initialize protocol matrices
    Dtilde = 0;
    Atilde = 0;

    % frame of algorithm
    if C == eye(n)
    % full-state coupling
        Atilde = [kron(eye(N),A) zeros(N*n, N*n) zeros(N*n, n); ...
            kron(eye(N) - D_bar, A) kron(D_bar, A) -kron(ones(N,1) - d_hat, A); ...
            zeros(n,N*n) zeros(n,N*n) A];
        Dtilde = -rho*[zeros(N*n, N*n) kron(eye(N), B*K) zeros(N*n, n); ...
            zeros(N*n, N*n) kron(eye(N), B*K) zeros(N*n,n); ...
            zeros(n,N*n) zeros(n, N*n) zeros(n,n)];
    else
    % partial-state coupling
        Atilde = [kron(eye(N),A) zeros(N*n, N*n) zeros(N*n, N*n) ...
            zeros(N*n, n); kron(eye(N) - D_bar, F*C) kron(eye(N), A - F*C) ...
            zeros(n*N, n*N) kron(d_hat - ones(N,1), F*C); ...
            zeros(n*N, n*N) eye(n*N) kron(eye(N) - D_bar, A) zeros(N*n,n);
            zeros(n,N*n) zeros(n,N*n) zeros(n,N*n) A];
        Dtilde = -rho*[zeros(N*n, N*n) kron(eye(N), B*K) zeros(N*n,N*n) ...
            zeros(N*n,n); zeros(N*n, N*n) kron(eye(N) - D_bar, B*K) zeros(N*n, N*n) ...
            zeros(N*n, n); zeros(N*n, N*n) kron(eye(N), B*K) zeros(N*n,N*n) ...
            zeros(N*n,n); zeros(n, N*n) zeros(n,N*n) zeros(n,N*n) ...
            zeros(n,n)];
    end
    
    % for each k, get x(k)
    for k = 0:K_max-1
        % for early values of k, use initial conditions
        if k <= kappa_bar
            if k == 0
                X(:,k+1) = Atilde*initial_conditions(:,kappa_bar+1) + ...
                    Dtilde*initial_conditions(:,k+1);
            else
                X(:,k+1) = Atilde*X(:,k) + Dtilde*initial_conditions(:,k+1);
            end
        else
            % delay x
            X_delay = zeros(3*N*n + n, 1);
            for j = 1:N
                % define R_j
                R_j = blkdiag(zeros((j-1)*n, (j-1)*n), eye(n), ...
                    zeros((N-1)*n, (N-1)*n), eye(n), zeros((N-j)*n, (N-j)*n));
                if C == eye(n) % full state
                    R_j = [R_j zeros(2*N*n,n); zeros(n,2*N*n) zeros(n,n)];
                else % partial state
                    R_j = [R_j zeros(2*N*n,n); zeros(n,2*N*n) zeros(n,n)];
                    R_j = [zeros(N*n, 3*N*n +n); zeros(2*N*n + n, N*n) R_j];
                end
                X_delay = X_delay + R_j*X(:, k-kappas(j));
            end
            % update x(k) to get x(k+1) via the linear system
            X(:, k+1) = Atilde*X(:,k) + Dtilde*X_delay;
        end
    end

    % solve the system
    
    if C == eye(n)
        % full state
        x_r = X(2*N*n+1:2*N*n+n,:);
        x = X(1:N*n,:);
    else 
       %  partial state
        x_r = X(3*N*n+1:3*N*n + n, :);
        x = X(1:N*n,:);
    end
    
    % get input
    % u = -rho*kron(eye(N), B*B'*P)*chi;

end
