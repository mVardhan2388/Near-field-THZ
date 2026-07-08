function A = gen_channel_thz_mimo(para)
% returns A stacked as (Mr*K) x Mt, LOS per-element spherical waves

Mt = para.M_t;
Mr = para.M_r;
K  = para.K;
lambda = para.lambda;

% absorption coefficient (1/m). set to realistic or 0
if isfield(para,'k_abs')
    k_abs = para.k_abs;
else
    k_abs = 0; % set 0 for qualitative sims, or e.g. 0.02 for THz
end

A = zeros(Mr*K, Mt);

% Tx element positions (assume ULA on y axis, x=0)
y_t = para.y_bs_t(:).';         % 1 x Mt
x_t = zeros(1, Mt);             % x coordinate = 0

% iterate users
for k = 1:K
    r_k = para.r(k);
    th = para.theta_users(k);

    % user centroid coordinates (x along range, y lateral)
    x0 = r_k * cos(th);
    y0 = r_k * sin(th);

    % Rx element coordinates for this user (assume small ULA around (x0,y0))
    % Rx positions: x = x0, y = y0 + para.y_bs_r
    y_r = para.y_bs_r(:).';     % 1 x Mr offsets about center
    x_r = x0 * ones(1, Mr);
    y_r = y0 + y_r;             % actual y positions of Rx elements

    Hk = zeros(Mr, Mt);
    for rr = 1:Mr
        for tt = 1:Mt
            d_rt = sqrt( (x_r(rr) - x_t(tt))^2 + (y_r(rr) - y_t(tt))^2 );
            % amplitude (use 1/d and absorption) and phase (spherical)
            amp = (1 / d_rt) * exp(- k_abs * d_rt);            % Thz extension to previous model
            phase = exp(-1j * 2*pi * d_rt / lambda);
            Hk(rr, tt) = amp * phase;
        end
    end

    Q(k) = min([Mr,Mt,rank(Hk)]);

    A((k-1)*Mr + 1 : k*Mr, :) = Hk;
end

para.Qk = min(Q);
end
