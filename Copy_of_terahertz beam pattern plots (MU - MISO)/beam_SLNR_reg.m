function W = beam_SLNR_reg(para, A)
% Regularized SLNR precoder for THz MU-MIMO
% A : (Mr*K) x Mt
% W : Mt x (Qk*K)

Mt = para.M_t;
Mr = para.M_r;
K  = para.K;
Qk = para.Qk;

sigma2 = para.noise;

% regularization strength (THz-safe default)
eta = 0.05;   % 0.01–0.1 typical

W = zeros(Mt, Qk*K);

for k = 1:K

    % Desired channel
    Hk = A((k-1)*Mr+1 : k*Mr, :);    % Mr x Mt
    R_sig = Hk' * Hk;               % Mt x Mt

    % Interference covariance
    R_int = zeros(Mt, Mt);
    for j = 1:K
        if j == k, continue; end
        Hj = A((j-1)*Mr+1 : j*Mr, :);
        R_int = R_int + Hj' * Hj;
    end

    % Regularization term
    alpha = eta * trace(R_int) / Mt;

    % Denominator matrix
    R_den = R_int + (sigma2 + alpha) * eye(Mt);

    % Generalized eigenvalue problem
    [V, D] = eig(R_sig, R_den, 'vector');
    [~, idx] = sort(real(D), 'descend');
    V = V(:, idx);

    % Select streams
    Wk = V(:, 1:Qk);

    % Normalize per-user
    Wk = Wk / norm(Wk, 'fro');

    % Store
    W(:, (k-1)*Qk+1 : k*Qk) = Wk;
end

% Global normalization
W = W / norm(W, 'fro');
end
