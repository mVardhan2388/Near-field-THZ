function atx = steer_tx(para, r, theta)
% returns 1 x Mt transmit steering (near-field spherical wave phases only)
lambda = para.lambda;
y_t = para.y_bs_t(:).';      % 1 x Mt
x_t = zeros(1, para.M_t);

% target point coordinates
x0 = r * cos(theta);
y0 = r * sin(theta);

% distances Tx element -> target point (point not Rx array center)
d = sqrt( (x0 - x_t).^2 + (y0 - y_t).^2 );
atx = exp(-1j * 2*pi * d / lambda);   % 1 x Mt
end
