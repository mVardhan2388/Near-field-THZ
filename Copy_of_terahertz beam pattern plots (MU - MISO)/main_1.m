clc;
clear;
close all;

% Load / build system -----------------------------
para = para_init();           % user-provided
Mt   = para.M_t;
Mr   = para.M_r;
K    = para.K;
Qk   = para.Qk;

user_labels = arrayfun(@(i) sprintf('User%d (r=%dm, \\theta=%d°)', ...
    i, round(para.r(i)), round(rad2deg(para.theta_users(i))) ), ...
    1:para.K, 'UniformOutput', false);

A = gen_channel_thz_mimo(para);    % (Mr*K) x Mt

W_BD = beam_BD(para, A);

x = linspace(1,40,200);   % range direction (m)
y = linspace(-20,20,200); % lateral direction (m)
[X,Y] = meshgrid(x,y);

R = sqrt(X.^2 + Y.^2);
Theta = atan2(Y, X);



spectrum = zeros(size(X));

for ix = 1:size(X,1)
    for iy = 1:size(X,2)

        r_test = R(ix,iy);
        theta_test = Theta(ix,iy);

        if r_test < 0.5
            spectrum(ix,iy) = eps;
            continue;
        end

        % Near-field channel to a virtual RX point
        Htest = gen_test_channel(para, r_test, theta_test); % Mr x Mt

        % Choose precoder (BD or SLNR)
        W = W_BD;     % <-- change to W_SLNR for comparison

        Heff = Htest * W; 

        if para.k_abs == 0
            PL_scalar = (1 / r_test)^2;    
        else
            PL_scalar = (1 / r_test)^2 * exp(-para.k_abs * r_test);  % include absorption
        end

        Heff_nopl = Heff / sqrt(PL_scalar);   
        s = sum((Heff_nopl));
        g = (s(1))^2;

        % Received power
        spectrum(ix,iy) = norm(g, 'fro');
    end
end

% Normalize
spectrum = spectrum / max(spectrum(:));

figure; colormap jet;
mesh(X, Y, 10*log10(spectrum));
view([60,60]);
xlim([0,40]);
ylim([-20,20]);
xlabel('x (m)');
ylabel('y (m)');
zlabel('Received Power (dB)');
title('Near-Field Power Focusing (THz MU-MIMO)');

