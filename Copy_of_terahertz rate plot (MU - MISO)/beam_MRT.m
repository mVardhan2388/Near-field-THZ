function [W] = beam_MRT(para, A)
% MRT precoder for near-field MU-MIMO
% A: (Mr*K) x Mt

Mt = para.M_t;
Mr = para.M_r;
K  = para.K;
Qk = para.Qk;

W = zeros(Mt, Qk*K);

for k = 1:K
    Hk = A((k-1)*Mr+1 : k*Mr, :);   % Mr x Mt

    % MRT direction
    Vk = Hk';                       % Mt x Mr

    % Take strongest Qk directions
    [U,~,~] = svd(Vk,'econ');
    Wk = U(:,1:Qk);

    % Normalize
    Wk = Wk / norm(Wk,'fro');
    % Wk = proR

    W(:,(k-1)*Qk+1:k*Qk) = Wk;
end
end
