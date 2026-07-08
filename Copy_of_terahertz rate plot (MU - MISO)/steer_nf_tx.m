function a_tx = steer_nf_tx(para, r, theta)

lambda = para.lambda;
y = para.y_bs_t;                      % Mt positions

dist_elem = sqrt( (r*sin(theta)).^2 + (r*cos(theta) - y).^2 );

a_tx = exp(-1j*2*pi * dist_elem/lambda).';   % 1×Mt

end
