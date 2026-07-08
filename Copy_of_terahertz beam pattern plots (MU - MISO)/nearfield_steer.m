function a_nf = nearfield_steer(para, r, theta)
    Mt = para.M_t;
    lambda = para.lambda;

    % Tx element coordinates
    x_t = zeros(1,Mt);
    y_t = para.y_bs_t(:).';

    % User location
    x0 = r*cos(theta);
    y0 = r*sin(theta);

    % Distance from each element to user focal point
    d = sqrt( (x_t - x0).^2 + (y_t - y0).^2 );

    % Near-field steering (phase only)
    a_nf = exp(-1j * 2*pi/lambda * d);
    a_nf = a_nf(:) / norm(a_nf);   % normalize
end
