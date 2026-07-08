function F = beam_mimo_mrt(para)

Mt = para.M_t;
K  = para.K;
lambda = para.lambda;

% --- enforce column vectors ---
x_t = para.x_bs_t(:);   % Mt x 1
y_t = para.y_bs_t(:);   % Mt x 1

F = zeros(Mt, K);

for k = 1:K
    r     = para.r_user(k);
    theta = para.theta_users(k);

    x0 = r * cos(theta);
    y0 = r * sin(theta);

    % --- distance from each Tx element to user ---
    d = sqrt( (x_t - x0).^2 + (y_t - y0).^2 );   % Mt x 1

    % --- near-field matched MRT (phase-only) ---
    a = exp(-1j * 2*pi/lambda * d);              % Mt x 1
    F(:,k) = a / norm(a);
end

% --- total power normalization ---
scale = sqrt(1 / trace(F*F'));
F = scale * F;

end
