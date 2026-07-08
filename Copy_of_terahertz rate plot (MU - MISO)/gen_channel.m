function [H, A_tx] = gen_channel(para)
H = zeros(para.M_r, para.M_t, para.K);
A_tx = zeros(para.K, para.M_t);

k_abs = get_absorption_coef(para.f);

for k = 1:para.K
    r_k = para.r(k);
    th_k = para.theta_users(k);

    a_t = steer_nf_tx(para, r_k, th_k);      % Mt × 1
    a_r = steer_nf_rx(para, r_k, th_k, k);   % Mr × 1

    PL = (4*pi*r_k/para.lambda)^2 * exp(2 * k_abs * r_k);
    beta = 1/sqrt(PL);

    H(:,:,k) = beta * (a_r * a_t.');         % Mr × Mt
    A_tx(k,:) = a_t.';                       % 1 × Mt
end
end
