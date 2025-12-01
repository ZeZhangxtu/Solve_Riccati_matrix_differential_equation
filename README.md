# MSA Solver for Riccati Matrix Differential Equation

A MATLAB implementation of the **Multi-Step Simplified Algorithm (MSA)** for solving the Riccati Matrix Differential Equation (RMDE) with a terminal condition. 

This solver is designed for high-precision numerical integration of matrix differential equations arising in optimal control, specifically for Finite-Time Linear Quadratic Regulator (LQR) problems.

## 📝 Mathematical Formulation

This solver addresses the following **Riccati Matrix Differential Equation (RMDE)**:

$$
-\dot{P}(t) = A^T P(t) + P(t) A - P(t) B R^{-1} B^T P(t) + C^T C, \quad t \in [t_0, t_f]
$$

Subject to the terminal condition:
$$
P(t_f) = F
$$

Where:
* $P(t) \in \mathbb{R}^{n \times n}$ is the unknown symmetric solution matrix.
* $A, B, C$ are system matrices ($Q = C^T C$).
* $R$ is the positive definite control weighting matrix.
* $F$ is the terminal state weighting matrix.

The algorithm integrates **backward** from $t_f$ to $t_0$.

## 🚀 Features

* **Algorithm**: Implements the **Multi-Step Simplified Algorithm (MSA)** (referenced as "Algorithm 4" in associated literature).
* **Optimization**: Utilizes algebraic refinements and matrix shifts (via the Continuous Algebraic Riccati Equation solution $P^-$) to ensure numerical stability and convergence.
* **Efficiency**: Computes solutions over a time grid $N$ using optimized matrix operations (Lyapunov sum computation).

## 📦 Installation & Requirements

### Prerequisites
* **MATLAB** (Tested on R2020a and later)
* **Control System Toolbox** (Required for the `care` function)

### Setup
1.  Clone this repository:
    ```bash
    git clone [https://github.com/ZeZhangxtu/Solve_Riccati_matrix_differential_equation.git](https://github.com/ZeZhangxtu/Solve_Riccati_matrix_differential_equation.git)
    ```
2.  Add the folder to your MATLAB path or navigate to the folder in MATLAB.

## 💻 Usage Example

Here is a simple script to demonstrate how to use the solver:

```matlab
% 1. Define System Matrices (Example: 2D System)
A = [0 1; -2 -3];
B = [0; 1];
C = [1 0];       % Note: The solver computes Q = C'*C internally
R = 1;           % Control weighting (scalar or matrix)
F = eye(2);      % Terminal condition P(tf) = F

% 2. Define Time and Solver Parameters
t0 = 0;          % Start time
tf = 5;          % End time
N  = 1000;       % Number of time steps
epsilon = 1e-9;  % Precision tolerance
l = 2;           % Number of refinement iterations

% 3. Call the MSA Solver
[P_solution, t_nodes] = MSA_RMDE_solver(A, B, C, R, F, t0, tf, N, epsilon, l);

% 4. Visualize Results (Example: Plotting element P_11 over time)
figure;
P11 = squeeze(P_solution(1, 1, :));
plot(t_nodes, P11, 'LineWidth', 2);
xlabel('Time (t)');
ylabel('P_{11}(t)');
title('Solution of RMDE using MSA');
grid on;
