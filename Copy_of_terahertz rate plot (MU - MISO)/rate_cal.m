function R = rate_cal(para, A, W)

Mr = size(A,1);
sigma2 = para.noise;
P = para.Pt;
M_t = para.M_t;


Heff = A * W;                  % Mr x Q
C = eye(Mr) + (P/(M_t*sigma2)) * (Heff * Heff');


C = (C + C')/2;

% Rate
R = real(log2(det(C)));

end
