function [P_solution, t_nodes] = solve_rmde_msa(A, B, C, R, F, t0, tf, N, epsilon, l)
% SOLVE_RMDE_MSA Solves the Matrix Differential Riccati Equation (RMDE).
%
% Algorithm:
%   Matrix Series Approximation (MSA). This method solves the RMDE by decomposing
%   the solution into an algebraic stabilizing part and a transient part. The
%   algebraic term is approximated using a matrix series expansion with
%   spectral scaling to ensure convergence within a specified tolerance.
%
% Equation:
%   -dP(t)/dt = A'P(t) + P(t)A - P(t)BR^{-1}B'P(t) + C'C
%   Terminal Condition: P(tf) = F
%
%   The solver integrates backward from t = tf to t = t0.
%
% Inputs:
%   A       : (n x n) State transition matrix.
%   B       : (n x s) Input matrix.
%   C       : (p x n) Output matrix.
%   R       : (s x s) Control weighting matrix (Symmetric positive definite).
%   F       : (n x n) Terminal condition matrix at t = tf.
%   t0, tf  : Start and end time of the integration interval (t0 < tf).
%   N       : Number of time discretization steps.
%   epsilon : Precision tolerance for the matrix series approximation.
%   l       : Number of refinement iterations for the algebraic solution.
%
% Outputs:
%   P_solution : (n x n x N) Solution trajectory P(t). P_solution(:,:,j)
%                corresponds to the solution at time t_nodes(j).
%   t_nodes    : (1 x N) Time vector corresponding to the integration steps.

% --- 1. Initialization and Parameter Setup ---
Q = C.' * C;                 % State weighting matrix
na_dim = size(A, 1);
nr_dim = size(R,1);
Ia = eye(na_dim);
Ir = eye(nr_dim);

% Compute control term coefficients
R_inv = Ir / R;
BRB = B * R_inv * B.';

delta_T = (tf - t0) / N;

% --- 2. Compute Stabilizing Algebraic Solution ---
% Solve the Algebraic Riccati Equation (ARE) for the stabilizing solution P_a-.
% Note: The system matrix is negated (-A) to align with the formulation
% required for the backward-time solution properties.
iP_a_minus = care(-A, B, Q, R);
P_a_minus = -iP_a_minus;

% --- 3. Closed-Loop System and Spectral Transformation ---
% Compute the closed-loop matrix A_check
A_check = A - BRB * P_a_minus;

% Perform spectral analysis for scaling parameters
eig_A_check = eig(A_check);
abs_eig_A_check = abs(eig_A_check);
max_eig_val = max(abs_eig_A_check);
min_eig_val = min(abs_eig_A_check);

% Calculate geometric mean of eigenvalues for scaling factor q
q = sqrt(max_eig_val * min_eig_val);
qI = q * Ia;

% Construct the transformed matrix A_hat via Cayley-like transform
% A_hat = (qI + A_check)^{-1} * (qI - A_check)
qI_plus_Acheck_inv = Ia / (qI + A_check);
A_hat = qI_plus_Acheck_inv * (qI - A_check);
A_hat_T = A_hat.';

% Transform the quadratic term Q
Q_hat = 2 * q * qI_plus_Acheck_inv * BRB * qI_plus_Acheck_inv.';

% --- 4. Matrix Series Approximation Parameters ---
% Determine the number of terms (n_guess) required for convergence
% based on the spectral radius (rho) and tolerance (epsilon).
rho_A_hat = max(abs(eig(A_hat)));
xi = (1 - rho_A_hat) / 2;
beta = (rho_A_hat + xi)^2;
norm_Q_hat_F = norm(Q_hat);

% Estimate iteration count
log_num = log(epsilon * (1 - beta) / norm_Q_hat_F);
log_den = log(beta);
n_guess = floor(log_num / log_den) - l;

% Refine n_guess to satisfy the norm condition ||A_hat^n|| < limit
A_n_guess = A_hat^n_guess;
condition_val_guess = (rho_A_hat + xi)^n_guess;
while norm(A_n_guess) >= condition_val_guess
    n_guess = n_guess + 1;
    A_n_guess = A_n_guess * A_hat;
    condition_val_guess = condition_val_guess * (rho_A_hat + xi);
end

% --- 5. Compute Approximate Algebraic Term (S_tilde) ---
% Decompose series summation indices
m_0_prime = floor(sqrt(n_guess)) + 1;
n_0_prime = floor(sqrt(n_guess)) + 1;

% Compute partial sums using the Lyapunov sum helper function
A_n0_prime = A_hat^n_0_prime;
S_q0 = Q_hat;

S_n0_standard = compute_lyap_sum(A_hat, Q_hat, n_0_prime);
S_n0 = A_hat * S_n0_standard * A_hat_T;

S_m0_standard = compute_lyap_sum(A_n0_prime, S_n0, m_0_prime-1);
S_m0 = A_n0_prime * S_m0_standard * A_n0_prime.';

% Combine terms to form the initial approximation
S_tilde_0 = S_q0 + S_n0 + S_m0;

% --- 6. Refine Approximation (Lyapunov Iterations) ---
P_hat_l = S_tilde_0;
for i = 1:l
    P_hat_l = A_hat * P_hat_l * A_hat_T + Q_hat;
end

% --- 7. Backward Recursive Integration ---
t_nodes = t0 + (0:N-1) * delta_T;
P_solution = zeros(na_dim, na_dim, N);

% The calculated P_hat_l serves as the algebraic term E
E_term = P_hat_l;

% Pre-compute matrix exponential for time stepping
X = expm(-A_check * delta_T);
X_T = X.';

% Transform terminal condition F to the initial recursive state
P_next = Ia / (F - P_a_minus);

for j = N:-1:1
    % Update state via the recursive formula
    P_current = X * (P_next - E_term) * X_T + E_term;

    % Recover the solution P(t) in the original coordinates
    P_solution(:, :, j) = Ia / P_current + P_a_minus;

    P_next = P_current;
end
end

function S_sum = compute_lyap_sum(A, Q, L)

A_T = A';
S_sum = Q;
P_term = Q;

for k = 1:L-1
    P_term = A * P_term * A_T;
    S_sum = S_sum + P_term;
end

end
