function [P_solution, t_nodes] = MSA_RMDE_solver(A, B, C, R, F, t0, tf, N, epsilon, l)
% [P_solution, t_nodes] = MSA_RMDE_solver(A, B, C, R, F, t0, tf, N, epsilon, l)
%
% Solves the matrix Riccati differential equation (RMDE)
% with a terminal condition using the Mmatrix Series Approximation (MSA).
%
% The equation being solved is:
%    -dP/dt = A'P + PA - P(B*R_inv*B')P + C'C
%    P(tf) = F
%
% The function integrates backward from t_f to t_0.
%
% Input:
%   A         : (n x n) State matrix
%   B         : (n x s) Input matrix
%   C         : (p x n) Output matrix
%   R         : (s x s) Control weighting matrix (symmetric positive definite)
%   F         : (n x n) Terminal condition matrix (at t=tf)
%   t0, tf    : Integration start and end time (t0 < tf)
%   N         : Number of time steps
%   epsilon   : Precision tolerance for the algorithm
%   l         : Number of refinement iterations for the algebraic solution
%
% Output:
%   P_solution : (n x n x N) Solution matrix P(t) at N time nodes.
%                P_solution(:,:,j) is the solution at t_nodes(j).
%   t_nodes    : (1 x N) Time nodes corresponding to the solution (from t0 to tf - delta_T).
%
% Remarks:
%   This implementation is based on "Algorithm 4" and related equations (e.g., (23), (24)),
%   solving via backward recurrence.
%
% Example:
%   % (User should provide specific values for A, B, C, R, F, t0, tf, N, epsilon, l)
%   % [P, t] = MSA_RMDE_solver(A, B, C, R, F, 0, 10, 1000, 1e-9, 2);
%
%   Author                : [Ze Zhang]
%   E-mail                : [zezhang@smail.xtu.edu.cn]
%   Last modification     : [2025/11/01]

Q = C.' * C;
na_dim = size(A, 1);
nr_dim = size(R,1);
Ia = eye(na_dim);
Ir = eye(nr_dim);
R_inv = Ir / R;

BRB = B * R_inv * B.';
delta_T = (tf - t0) / N;

iP_a_minus = care(-A, B, Q, R);
P_a_minus = -iP_a_minus; % P_a_minus is the negative definite solution

A_check = A - BRB * P_a_minus;

eig_A_check = eig(A_check);
abs_eig_A_check = abs(eig_A_check);

max_eig_val = max(abs_eig_A_check);
min_eig_val = min(abs_eig_A_check);
q = sqrt(max_eig_val * min_eig_val);

qI = q * Ia;
qI_plus_Acheck_inv = Ia / (qI + A_check);

A_hat = qI_plus_Acheck_inv * (qI - A_check);
A_hat_T = A_hat.';
Q_hat = 2 * q * qI_plus_Acheck_inv * BRB * qI_plus_Acheck_inv.';

rho_A_hat = max(abs(eig(A_hat)));
xi = (1 - rho_A_hat) / 2;
beta = (rho_A_hat + xi)^2;

norm_Q_hat_F = norm(Q_hat, 'fro');

log_num = log(epsilon * (1 - beta) / norm_Q_hat_F);
log_den = log(beta);
n_guess = floor(log_num / log_den) - l;
n =  n_guess;
n_0_guess = floor(sqrt(n_guess)) + 1;
A_n0_guess = A_hat^n_0_guess;
A_n0_prime = A_n0_guess;
condition_val_guess = (rho_A_hat + xi)^n_0_guess;

max_iter = 1000;
iter_count = 0;

while norm(A_n0_guess, 'fro') >= condition_val_guess
    n_0_guess = n_0_guess + 1;
    A_n0_guess = A_n0_guess * A_hat;
    condition_val_guess = condition_val_guess * (rho_A_hat + xi);
    if n_0_guess > n_guess
        n = n_0_guess;
    else
        n = n_guess;
    end
    iter_count = iter_count + 1;
    if iter_count > max_iter
        error("MSA_RMDE_solver:LoopError", "The while loop for determining 'n' exceeded %d iterations.", n_guess);
    end

end

m_0_prime = floor(sqrt(n)) + 1;
n_0_prime = floor(sqrt(n)) + 1;
A_n0_prime = A_n0_prime * A^(floor(sqrt(n)) + 1 - floor(sqrt(n_guess)) - 1);

S_q0 = Q_hat; % q_0_prime = 0; % As defined in the pseudo-code q_0' = 0

S_n0_standard = compute_lyap_sum(A_hat, Q_hat, n_0_prime);
S_n0 = A_hat * S_n0_standard * A_hat_T;

S_m0_standard = compute_lyap_sum(A_n0_prime, S_n0, m_0_prime-1);
S_m0 = A_n0_prime * S_m0_standard * A_n0_prime.';

S_tilde_0 = S_q0 + S_n0 + S_m0;

P_hat_l = S_tilde_0;
for i = 1:l
    P_hat_l = A_hat * P_hat_l * A_hat_T + Q_hat;
end

t_nodes = t0 + (0:N-1) * delta_T;
P_solution = zeros(na_dim, na_dim, N);

E_term = P_hat_l;

X = expm(-A_check * delta_T);
X_T = X.';

P_next = Ia / (F - P_a_minus);

for j = N:-1:1
    P_current = X * (P_next - E_term) * X_T + E_term;
    P_solution(:, :, j) = Ia / P_current + P_a_minus;
    P_next = P_current;
end

end


