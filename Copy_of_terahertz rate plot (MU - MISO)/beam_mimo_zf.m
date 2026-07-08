function F = beam_mimo_zf(para, H)
% compute per-user effective row by taking first left singular vector as combiner
Heq = zeros(para.K, para.M_t);
for k=1:para.K
    Hk = squeeze(H(:,:,k));          % Mr x Mt
    [U,~,~] = svd(Hk,'econ');
    u1 = U(:,1);                     % Mr x 1
    Heq(k,:) = (u1' * Hk);           % 1 x Mt effective channel
end
% now do ZF on Heq (K x Mt)
Fzf = pinv(Heq);   % Mt x K (pseudo-inverse)
F = Fzf;
% scale to total power
scale = sqrt( real(para.Pt / (trace(F*F') + 1e-15)) );
F = scale * F;
end
