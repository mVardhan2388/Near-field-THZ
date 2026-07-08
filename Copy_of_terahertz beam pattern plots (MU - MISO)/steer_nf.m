function A = steer_nf(para, r, theta)

M_r = para.M_r;
lambda = para.lambda;
y_t = para.y_bs_t;
y_r = para.y_bs_r(1:M_r);

A = zeros(M_r, para.M_t);

for mr = 1:M_r
    for mt = 1:para.M_t
        d = sqrt( (r*sin(theta))^2 + (r*cos(theta)-y_t(mt))^2 + (y_r(mr))^2 );
        A(mr, mt) = exp(-1j*2*pi*d/lambda);
    end
end
end
