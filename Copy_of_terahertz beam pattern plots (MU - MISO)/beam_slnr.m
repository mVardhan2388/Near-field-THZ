function W = beam_slnr(para, A)
% Multi-user near-field SLNR precoder
% A : (Mr*K) x Mt  stacked channel
% W : Mt x (Qk*K)

Mt = para.M_t;
Mr = para.M_r;
K  = para.K;
Qk = para.Qk;

sigma2 = para.noise;     % thermal noise power

W = zeros(Mt, Qk*K);

for k = 1:K

    % -------- Desired user channel ----------
    Hk = A((k-1)*Mr + 1 : k*Mr, :);     % Mr x Mt

    % -------- Interference covariance -------
    R_int = zeros(Mt, Mt);
    for j = 1:K
        if j == k, continue; end
        Hj = A((j-1)*Mr + 1 : j*Mr, :);
        R_int = R_int + Hj' * Hj;
    end

    % -------- SLNR matrices -----------------
    R_sig = Hk' * Hk;
    R_den = R_int + sigma2 * eye(Mt);

    % -------- Generalized eigenvalue problem
    % R_sig v = lambda R_den v
    [V, D] = eig(R_sig, R_den, 'vector');

    % sort descending SLNR
    [~, idx] = sort(real(D), 'descend');
    V = V(:, idx);

    % take Qk strongest streams
    Wk = V(:, 1:Qk);

    % per-user normalization
    Wk = Wk / norm(Wk, 'fro');

    % store
    W(:, (k-1)*Qk + 1 : k*Qk) = Wk;
end

% -------- Global power normalization --------
W = W / norm(W, 'fro');
end
