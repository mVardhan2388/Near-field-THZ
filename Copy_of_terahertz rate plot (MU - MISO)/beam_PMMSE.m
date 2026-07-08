function W = beam_PMMSE(para, A)
% -------------------------------------------------------
% Stacked MMSE precoder for near-field THz MU-MIMO
% - Single-stream MMSE per user
% - Replicated (stacked) to K*Qk columns for fair comparison
% -------------------------------------------------------

Mt = para.M_t;
Mr = para.M_r;
K  = para.K;
Qk = para.Qk;

Pt    = para.Pt;
noise = para.noise;

%% ===== Build effective MU channel =====
% Each row = one effective user channel
H_eff = zeros(K, Mt);   % K x Mt

for k = 1:K
    Hk = A((k-1)*Mr+1 : k*Mr, :);    % Mr x Mt

    % Near-field stable combining (array-average)
    hk = mean(Hk, 1);               % 1 x Mt

    H_eff(k,:) = hk;
end

%% ===== MMSE precoder (single-stream per user) =====
Reg = (noise / Pt) * eye(K);

t_W_mmse = H_eff' / (H_eff*H_eff' + Reg);   % Mt x K
W_mmse = sqrt(Mt/trace(t_W_mmse*t_W_mmse')) * t_W_mmse ;

%% ===== Stack to match K*Qk streams =====
W = zeros(Mt, K*Qk);

for k = 1:K
    wk = W_mmse(:,k);     % Mt x 1

    % Replicate across Qk streams
    W(:, (k-1)*Qk+1 : k*Qk) = repmat(wk, 1, Qk);
end

%% ===== Power normalization =====
W = W / (sqrt(trace(W'*W)) * sqrt(Pt));

end
