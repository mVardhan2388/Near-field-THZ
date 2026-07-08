function F = beam_mimo_slnr(para, H)
F = zeros(para.M_t, para.K);
sigma2 = para.noise;
for k=1:para.K
    Hk = squeeze(H(:,:,k));
    [U,~,~] = svd(Hk,'econ');
    u1 = U(:,1);                     % combiner
    a_eff = (u1' * Hk).';            % Mt x 1 effective transmit vector
    idx = setdiff(1:para.K, k);
    Interf = zeros(para.M_t);
    for j = idx
        Hj = squeeze(H(:,:,j));
        uj = svd(Hj,'econ'); uj = uj(:,1);
        a_j = (uj' * Hj).';
        Interf = Interf + (a_j * a_j');
    end
    M = Interf + sigma2 * eye(para.M_t);
    [v,~] = eigs( M \ (a_eff * a_eff'), 1);
    F(:,k) = v;
end
scale = sqrt( real(para.Pt / (trace(F*F') + 1e-15)) );
F = scale * F;
end
