% Discrete Protocol Design. This function accepts the following parameters:
% 
% 1. The agent model: A, B, C. Types: matrices.
% 2. An upper bound on the delays: kappa_bar. Type: positive integer.
% 
% From this information, discrete_protocol_design picks parameters that will
% go into the protocol described in the Discrete Time Scale-Free input 
% delay paper. Our metric of performance is given by speed of 
% synchronization, i.e., how quickly the agents x_i converge to the 
% exosystem x_r. Typically, larger epsilon lead to better performance. 
% As such, we have designed our algorithm with an eye towards larger 
% epsilon. The algorithm is designed as follows. First, we find the maximum
% upper bound on delays (kappabar_max) through the solvability condition 
% given in the paper. Then, we check that the given tau_bar is strictly less
% than kappabar_max. If the problem is solvable, we design the protocol,
% returning the following parameters:

% 1. epsilon. Type: positive number.
% 2. rho. Type: positive number.
% 3. taubar_max. Type: positive number.
% 4. P. Type: positive definite matrix.
% 5. K. Type: matrix.

% NOTE: This algorithm is meant to generate valid parameters (epsilon, rho) 
% given any agent model and valid upper bound on delays. We make no claim
% that this choice is optimal for convergence. In fact, you may often
% find that choosing epsilon larger than we give here will make performance
% much better. What we DO guarantee is that you need not choose epsilon any
% smaller than what this function returns. We welcome any improvements 
% upon this algorithm to get a larger allowable epsilon.

function [epsilon, rho, kappabar_max, K, F] = discrete_protocol_design(A,B,C,kappa_bar)
    
    % get dimensions
    sizes = size(A);
    n = sizes(1);
    
    % get omega_max from the eigenvalues of A
    eigenvalues = eig(A);
    eigenvalues = eigenvalues';
    omega_max = -realmax;
    for lambda = eigenvalues
        if abs(lambda) == 1
            temp = (lambda + conj(lambda))/2;
            omega = acos(temp);
            if omega > omega_max
                omega_max = omega;
            end
        end
    end

    % if A is Schur, set omega_max = 0
    if omega_max == -realmax
        omega_max = 0;
    end

    % get kappabar_max
    kappabar_max = pi/(2*omega_max);
    
    if kappa_bar >= kappabar_max
        error("kappa_bar is too large for the given A. Try a smaller value")
    end
    
    % check controllability
    if rank(ctrb(A,B)) < 0
        error ("(A,B) are not controllable.")
    end
    
    % check observability
    if rank(obsv(A,C)) < 0
        error("(A,C) are not observable.")
    end
  
    % fix rho
    rho = 1/cos(kappa_bar*omega_max);
    
    
    theta = (acos(1/2*rho)/kappa_bar) - omega_max;
    % algorithm to compute a non-conservative mu
    omega = omega_max + theta;
    omega_bar = pi;
    length = omega_bar - omega;
    step_size = length/100;
    mu = sqrt(min(eig((exp(-j*omega)*eye(n) - A')*(exp(j*omega)*eye(n) - A))));
    while abs(omega - omega_bar) > 10^(-5)*length
        % search the interval (omega_max + theta, omega_bar) to see
        % if the given mu minimizes sv
        for omega = omega_max + theta:step_size:omega_bar
            sv = sqrt(min(eig((exp(-j*omega)*eye(n) - A')*(exp(j*omega)*eye(n) - A))));
            if mu >= sv
                % mu does not minimize sv; update mu
                mu = mu*0.999;
                break
            end
        end
    end
   
    % algorithm to compute epsilon
    % initialize values
    epsilon = 1;
     % get K
    [P,V,K] = dare(A,B,epsilon*eye(n),1);
    size(rho)
    size(B)
    size(K)
    value = rho*norm(B*K);
    % repeat until value <= 0.99*mu
    while value > mu*0.99
        % update epsilon, value
        epsilon = 0.99*epsilon;
         % get K
        [P,V,K] = dare(A,B,epsilon*eye(n),1);
        value = rho*norm(B*K);
    end
    
    % algorithm to get F. Change this if you want, but it will not
    % affect the solution too much.
    p = [-0.5,0.5,0.1];
    F = place(A',C',p);
    F = F';
    
end