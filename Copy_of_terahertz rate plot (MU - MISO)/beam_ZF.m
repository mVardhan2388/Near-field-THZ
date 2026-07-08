function W = beam_ZF(para, A)
% ZF: W = A^H (A A^H)^(-1)
Mt = para.M_t;
K  = para.K;

W = A' * pinv(A*A');      % MT × K
W = W / norm(W,'fro');    % normalize

end
