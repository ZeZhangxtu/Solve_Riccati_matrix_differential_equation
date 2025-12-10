function [P_solution, t_nodes] = solve_rmde_msa(A, B, C, R, F, t0, tf, N, epsilon, l)
%SOLVE_RMDE_MSA Solves the Matrix Differential Riccati Equation using MSA.
%
%   [P_solution, t_nodes] = SOLVE_RMDE_MSA(A, B, C, R, F, t0, tf, N, epsilon, l)
%   computes the solution to the standard Riccati Matrix Differential Equation
%   using the Matrix Series Approximation (MSA) method.
%
%   Inputs:
%       A, B, C    - System matrices defining the state-space model.
%       R          - Control weighting matrix (positive definite).
%       F          - Terminal condition matrix P(tf) = F.
%       t0, tf     - Integration time interval [t0, tf].
%       N          - Number of time steps.
%       epsilon    - Tolerance for the approximation series.
%       l          - Number of refinement iterations (Lyapunov steps).
%
%   Outputs:
%       P_solution - 3D array of size (n, n, N) containing solution matrices.
%       t_nodes    - Vector of time points corresponding to the solution.

% Dimensions and Identity Matrices
Q = C.' * C;
na_dim = size(A, 1);
nr_dim = size(R, 1);
Ia = eye(na_dim);
Ir = eye(nr_dim);

% Control term coefficients
R_inv = Ir / R;
BRB = B * R_inv * B.';
delta_T = (tf - t0) / N;

% Compute Stabilizing Algebraic Solution (P_a-)
% Solve ARE for P_a-. System matrix is negated (-A).
iP_a_minus = care(-A, B, Q, R);
P_a_minus = -iP_a_minus;

% Closed-Loop System Construction
A_check = A - BRB * P_a_minus;

% Spectral Analysis and Scaling
eig_A_check = eig(A_check);
abs_eig_A_check = abs(eig_A_check);
max_eig_val = max(abs_eig_A_check);
min_eig_val = min(abs_eig_A_check);

% Geometric mean scaling factor q
q = sqrt(max_eig_val * min_eig_val);
qI = q * Ia;

% Cayley-like Transformation
% A_hat = (qI + A_check)^{-1} * (qI - A_check)
qI_plus_Acheck_inv = Ia / (qI + A_check);
A_hat = qI_plus_Acheck_inv * (qI - A_check);
A_hat_T = A_hat.';

% Transform Quadratic Term Q
Q_hat = 2 * q * qI_plus_Acheck_inv * BRB * qI_plus_Acheck_inv.';

% Estimation of Series Terms (n_guess)
rho_A_hat = max(abs(eig(A_hat)));
xi = (1 - rho_A_hat) / 2;
beta = (rho_A_hat + xi)^2;
norm_Q_hat = norm(Q_hat);

log_num = log(epsilon * (1 - beta) / norm_Q_hat);
log_den = log(beta);
n_guess = floor(log_num / log_den) - l;

% Refine n_guess using Norm Condition
if ~issymmetric(A_check)
    A_hat_schur = schur(A_hat);
    T_pow = A_hat_schur^n_guess;

    condition_val_guess = (rho_A_hat + xi)^n_guess;

    while norm(T_pow) >= condition_val_guess
        n_guess = n_guess + 1;
        T_pow = T_pow * A_hat_schur;
        condition_val_guess = condition_val_guess * (rho_A_hat + xi);
    end
end

% Compute Approximate Algebraic Term (S_tilde)
m_0_prime = floor(sqrt(n_guess)) + 1;
n_0_prime = floor(sqrt(n_guess)) + 1;

% Partial sums calculation
A_n0_prime = A_hat^n_0_prime;
S_q0 = Q_hat;
S_n0_standard = compute_lyap_sum(A_hat, Q_hat, n_0_prime);
S_n0 = A_hat * S_n0_standard * A_hat_T;
S_m0_standard = compute_lyap_sum(A_n0_prime, S_n0, m_0_prime-1);
S_m0 = A_n0_prime * S_m0_standard * A_n0_prime.';

S_tilde_0 = S_q0 + S_n0 + S_m0;

% Refinement (Lyapunov Iterations)
P_hat_l = S_tilde_0;
for i = 1:l
    P_hat_l = A_hat * P_hat_l * A_hat_T + Q_hat;
end

% Backward Recursive Integration
t_nodes = t0 + (0:N-1) * delta_T;
P_solution = zeros(na_dim, na_dim, N);

E_term = P_hat_l;
X = expm(-A_check * delta_T);
X_T = X.';

% Pre-compute constant offset
Constant_Offset = E_term - X * E_term * X_T;
P_next = Ia / (F - P_a_minus);

for j = N:-1:1
    P_current = X * P_next * X_T + Constant_Offset;

    % Recover solution
    P_solution(:, :, j) = Ia / P_current + P_a_minus;

    P_next = P_current;
end
end

function S_sum = compute_lyap_sum(A, Q, L)
%COMPUTE_LYAP_SUM Helper function to compute linear Lyapunov sum.
A_T = A';
S_sum = Q;
P_term = Q;
for k = 1:L-1
    P_term = A * P_term * A_T;
    S_sum = S_sum + P_term;
end
end
