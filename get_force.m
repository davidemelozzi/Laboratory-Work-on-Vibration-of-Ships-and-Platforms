function [F] = get_force(t)
    % Create vector of pulse forces

    T_impulse = 0.011;
    T_half = T_impulse / 2;
    Fmax = 10000;
    F = zeros(length(t), 1);

    for i = 1:length(t)

        if t(i) < T_half
            F(i) = Fmax * t(i) / T_half;
        elseif t(i) > T_impulse
            F(i) = 0;
        else
            F(i) = (Fmax * (t(i) - T_impulse)) / T_half;
        end

    end

end
