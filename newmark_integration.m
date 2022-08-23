function [X, Xdot, X2dot] = newmark_integration(t, F, X0, Xdot0, M, K, C, alfa, beta_)
    % Newmark Integration Method
    %--------------------------------------------------------------------------
    % Integrates a N-DOF system with a mass matrix "M", stiffness matrix "K" and
    % damping matrix "C" subjected to an external force matrix F.
    % Returns the displacement, velocity and acceleration of the system with
    % respect to an inertial frame of reference.
    %
    % Input
    % ----------
    %       [t] :       Time Vector                 [n, 1]
    %       [F] :       External Force Matrix       [n, DOF]
    %       [X0]:       Initial Position Vector     [1, DOF]
    %       [Xdot0]:    Initial Velocity Vector     [1, DOF]
    %       [M]:        Mass Matrix                 [DOF, DOF]
    %       [K]:        Stiffness Matrix            [DOF, DOF]
    %       [C]:        Damping Matrix              [DOF, DOF]
    %       [alfa]:     Integration Parameter       scalar
    %       [beta_]:    Integration Parameter       scalar
    %
    % Output
    % ----------
    %       [X]:        Displacement Response   [n, DOF]
    %       [Xdot]:     Velocity                [n, DOF]
    %       [X2dot]:    Acceleration            [n, DOF]

    % Initializing Variables
    dt = t(2) - t(1); DOF = length(K); n = length(t);
    X = zeros(DOF, n); Xdot = X; X2dot = X;

    % Initial Conditions
    X(:, 1) = X0; Xdot(:, 1) = Xdot0;
    X2dot(:, 1) = M \ (F(:, 1) - K * X(:, 1) - C * Xdot(:, 1));

    % Integration Constants
    a0 = 1 / (alfa * dt^2);
    a1 = beta_ / (alfa * dt);
    a2 = 1 / (alfa * dt);
    a3 = ((1 / (2 * alfa)) - 1);
    a4 = (beta_ / alfa) - 1;
    a5 = ((alfa * dt) / 2) * ((beta_ / alfa) - 2);
    a6 = beta_ * dt;
    a7 = (1 - beta_) * dt;

    % Form Effective Stiffness Matrix
    Keff = K + a0 * M + a1 * C;
    inv_Keff = inv(Keff);

    for i = 1:(n - 1)
        Feff = F(:, i + 1) + ...
            M * (a0 * X(:, i) + a2 * Xdot(:, i) + a3 * X2dot(:, i)) + ...
            C * (a1 * X(:, i) + a4 * Xdot(:, i) + a5 * X2dot(:, i));

        % ALERT: DO NOT CHANGE THE ORDER OF THE VARIABLES BELOW
        % Xdot(:, i + 1), for example, accesses the row X2dot(:, i + 1),
        % which is defined in the previous line
        X(:, i + 1) = inv_Keff * Feff;
        X2dot(:, i + 1) = a0 * (X(:, i + 1) - X(:, i)) - a2 * Xdot(:, i) - a3 * X2dot(:, i);
        Xdot(:, i + 1) = Xdot(:, i) + a7 * X2dot(:, i) + a6 * X2dot(:, i + 1);
    end

end
