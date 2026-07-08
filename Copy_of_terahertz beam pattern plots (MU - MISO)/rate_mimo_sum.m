function R = rate_mimo_sum(para, H, F)
R = 0;
for k=1:para.K
    Hk = squeeze(H(:,:,k));           % Mr x Mt
    % assume combiner is left singular vector u1
    [U,~,~] = svd(Hk,'econ');
    u1 = U(:,1);                       % Mr x 1
    hk_eq = u1' * Hk * F(:,k);         % scalar desired
    interf = 0;
    for j=1:para.K
        if j==k, continue; end
        interf = interf + abs(u1' * Hk * F(:,j))^2;
    end
    sig = abs(hk_eq)^2;
    snr = sig / (interf + para.noise);
    R = R + log2(1 + real(snr));
end
end
