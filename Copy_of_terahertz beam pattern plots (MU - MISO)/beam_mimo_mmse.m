function F = beam_mimo_mmse(para, H)
% build effective channels Heq as in ZF
Heq = zeros(para.K, para.M_t);
for k=1:para.K
    Hk = squeeze(H(:,:,k));
    [U,~,~] = svd(Hk,'econ');
    u1 = U(:,1);
    Heq(k,:) = (u1' * Hk);
end
sigma2 = para.noise;
Mt = para.M_t;
la = sigma2/para.Pt;   % or use sigma2
F0 = (Heq' * Heq + la * eye(Mt)) \ Heq';  % Mt x K
scale = sqrt( real(para.Pt / (trace(F0*F0') + 1e-15)) );
F = scale * F0;
end
