function a_rx = steer_nf_rx(para, r, theta, k)
lambda = para.lambda;

x_rx = squeeze(para.u(1,:,k)).';   % ensure column (Mr×1)
y_rx = squeeze(para.u(2,:,k)).';   % ensure column (Mr×1)

ux = r * cos(theta);
uy = r * sin(theta);

d_n = sqrt( (ux - x_rx).^2 + (uy - y_rx).^2 );
a_rx = exp(-1j * 2*pi * d_n / lambda);

a_rx = a_rx(:);                    % FORCE column (Mr×1)
end
