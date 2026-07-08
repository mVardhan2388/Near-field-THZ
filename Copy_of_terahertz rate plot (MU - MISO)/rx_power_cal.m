function Pr = rx_power_cal(Ak, Wk)
% Ak : Mr x Mt channel for user k
% Wk : Mt x Qk precoder for user k

    Pr = real(trace( Ak * Wk * Wk' * Ak' ));
end