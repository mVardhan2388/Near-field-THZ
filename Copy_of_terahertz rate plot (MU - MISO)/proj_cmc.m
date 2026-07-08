function W_cmc = proj_cmc(W)
% Element-wise constant modulus projection (paper definition)
% |w_n| = 1/sqrt(Mt)

Mt = size(W,1);
W_cmc = exp(1j * angle(W)) / sqrt(Mt);

end
