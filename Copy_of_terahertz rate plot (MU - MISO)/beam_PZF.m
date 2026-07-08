function W = beam_PZF(para, A)
% ZF: W = A^H (A A^H)^(-1)
Mt = para.M_t;
K  = para.K;

beta = sqrt(Mt/trace(inv(A)*(inv(A)')));

W = A' * pinv(A*A');      % MT × K

W = beta * W / norm(W,'fro');    % normalize

end
