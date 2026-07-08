function Htest = gen_test_channel(para, r, theta)
Mt = para.M_t; 
Mr = para.M_r;
lambda = para.lambda;

y_t = para.y_bs_t(:).'; 
x_t = zeros(1,Mt);

x0 = r*cos(theta);
y0 = r*sin(theta);

y_r = para.y_bs_r(:).'; 
x_r = x0*ones(1,Mr);
y_r = y0 + y_r;

Htest = zeros(Mr, Mt);
for rr=1:Mr
    for tt=1:Mt
        d = sqrt( (x_r(rr)-x_t(tt))^2 + (y_r(rr)-y_t(tt))^2 );
        Htest(rr,tt) = 1/d * exp(-1j*2*pi*d/lambda);
    end
end
end
