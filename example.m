clear all

% Time Vector
t0 = 0; dt = 0.011/2; tf = 0.23;
t = t0:dt:tf;

% Mass Matrix
M = [100 0 0;
    0 100 0;
    0 0 50];

% Stiffness Matrix
K = 1e7 * [2 -1 0;
        -1 2.5 -0.5;
        -0 -0.5 0.5];

% Damping Matrix
C = [5000 0 0;
    0 2500 -1000;
    0 -1000 1000];

% Initial Position Vector
X0 = [0 0 0];

% Initial Velocity Vector
Xdot0 = [0 0 0];

% Force Matrix
F(1, :) = 0 * t; F(2, :) = 0 * t; F(3, :) = get_force(t);

% Newmark Integration Constants
beta_ = 0.5; alfa = 0.25 * (0.5 + beta_)^2;

% Integrating
[X, Xdot, X2dot] = newmark_integration(t, F, X0, Xdot0, M, K, C, alfa, beta_);

% Converting displacement result from column to row for better display
X = X';
disp(X);

% Plotting
x1 = X(:, 1); x2 = X(:, 2); x3 = X(:, 3);
figure("name", "Displacement versus Time for 3DOF system solved using Newmark method");
plot(t, x1,
t, x2,
t, x3);
xlabel("Time [s]");
ylabel("Displacement [mm]");
