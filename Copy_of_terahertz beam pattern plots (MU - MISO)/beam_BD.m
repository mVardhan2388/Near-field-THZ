function W = beam_BD(para, A)
Mt = para.M_t; Mr = para.M_r; K = para.K; Qk = para.Qk;
W = zeros(Mt, Qk*K);

for k = 1:K
    rows = setdiff(1:K, k);
    H_int = zeros((K-1)*Mr, Mt);
    ptr = 1;
    for ii = 1:length(rows)
        r = rows(ii);
        H_int(ptr:ptr+Mr-1, :) = A((r-1)*Mr + 1 : r*Mr, :);
        ptr = ptr + Mr;
    end

    Vk0 = null(H_int);      % Mt x d
    d = size(Vk0,2);
    if d == 0
        Vk0 = eye(Mt); d = Mt;
    end

    Hk = A((k-1)*Mr+1 : k*Mr, :);    % Mr x Mt
    H_tilde = Hk * Vk0;              % Mr x d

    [~, ~, V2full] = svd(H_tilde, 'econ');
    Quse = min(Qk, d);
    if Quse >= 1
        V2 = V2full(:,1:Quse);
        Wk = Vk0 * V2;               % Mt x Quse
    else
        Wk = zeros(Mt, 0);
    end

    if Quse < Qk
        need = Qk - Quse;
        extra = [];
        if d > Quse
            avail = min(d - Quse, need);
            extra = Vk0(:, Quse+1:Quse+avail);
            need = need - avail;
        end
        if need > 0
            Rnd = randn(Mt, need) + 1i*randn(Mt, need);
            if ~isempty(Wk)
                Qfull = orth([Wk, Rnd]);
                extra = [extra, Qfull(:, size(Wk,2)+1 : size(Wk,2)+need)];
            else
                [Qfull,~] = qr(Rnd,0);
                extra = [extra, Qfull(:,1:need)];
            end
        end
        Wk = [Wk, extra];
    end

    if ~isempty(Wk)
        Wk = Wk / norm(Wk, 'fro');
    end

    col1 = (k-1)*Qk + 1; col2 = k*Qk;
    W(:, col1:col2) = Wk;
end
end
