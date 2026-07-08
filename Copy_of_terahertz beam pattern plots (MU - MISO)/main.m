% % main_BD.m
% clc; clear; close all;
% 
% para = para_init();
% Mr = para.M_r; Mt = para.M_t; K = para.K; Qk = para.Qk;
% 
% user_labels = arrayfun(@(i) sprintf('User%d (r=%dm, \\theta=%d°)', ...
%     i, round(para.r(i)), round(rad2deg(para.theta_users(i))) ), ...
%     1:para.K, 'UniformOutput', false);
% 
% % Channel
% A = gen_channel_thz_mimo(para);   % (Mr*K) × Mt
% 
% % BD Precoder
% W_BD = beam_BD(para, A);          % Mt × (Qk*K)
% W_BD = W_BD% / norm(W_BD,'fro');   % normalize
% 
% %% ---------------------- BEAM PATTERN VS DISTANCE ------------------------
% dist = linspace(1,25,1000);
% pattern_BD = zeros(length(dist), K);
% 
% for i = 1:length(dist)
%     for k = 1:K
% 
%         % Test channel (Mr × Mt)
%         Htest = gen_test_channel(para, dist(i), para.theta_users(k));
% 
%         % Extract user-k BD precoder (Mt × Qk)
%         Wk = W_BD(:, (k-1)*Qk + 1 : k*Qk);
% 
%         % Beampattern
%         pattern_BD(i,k) = norm( Htest * Wk )^2;
%     end
% end
% 
% pattern_BD_db = 20*log10(pattern_BD ./ max(pattern_BD(:)));
% 
% figure; hold on; grid on;
% plot(dist, pattern_BD_db, 'LineWidth', 1.4);
% xlabel('Distance (m)'); ylabel('Beam Pattern (dB)');
% title('THz MU-MIMO Beam Pattern vs Distance (BD)');
% legend(user_labels, 'Location','best');
% 
% %% ---------------------- BEAM PATTERN VS ANGLE ---------------------------
% ang = linspace(0,180,720);
% pattern_BD_ang = zeros(length(ang), K);
% 
% for i = 1:length(ang)
%     theta = deg2rad(ang(i));
%     for k = 1:K
% 
%         % Test Channel (Mr × Mt)
%         Htest = gen_test_channel(para, para.r(k), theta);
% 
%         % BD precoder for user k
%         Wk = W_BD(:, (k-1)*Qk + 1 : k*Qk);
% 
%         pattern_BD_ang(i,k) = norm( Htest * Wk )^2;
%     end
% end
% 
% pattern_BD_ang_db = 20*log10(pattern_BD_ang ./ max(pattern_BD_ang(:)));
% 
% figure; hold on; grid on;
% plot(ang, pattern_BD_ang_db, 'LineWidth', 1.4);
% xlabel('Angle (deg)'); ylabel('Beam Pattern (dB)');
% title('THz MU-MIMO Beam Pattern vs Angle (BD)');
% legend(user_labels, 'Location','best');



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

% Pattern vs Angle (fixed per-user distance) -----
ang_deg = linspace(-90,90,721);    
ang_rad = deg2rad(ang_deg);
pattern_ang = zeros(length(ang_deg), K);

for ia = 1:length(ang_rad)
    theta = ang_rad(ia);
    for k = 1:K
        
        Htest = gen_test_channel(para, para.r(k), theta);   % Mr x Mt

       
        Wk = W_BD(:, (k-1)*Qk + 1 : k*Qk);                 % Mt x Qk

       
        Heff = Htest * Wk;                                   % Mr x Qk

        
        s = svd(Heff);
        if isempty(s)
            g = 0;
        else
            g = (s(1))^2;
        end

        pattern_ang(ia, k) = g;
    end
end

% Normalize per-user and convert to dB
pattern_ang_db = zeros(size(pattern_ang));
eps_db = 1e-12;
for k = 1:K
    maxval = max(pattern_ang(:,k));
    if maxval <= 0
        pattern_ang_db(:,k) = -140 * ones(size(pattern_ang(:,k)));
    else
        pattern_ang_db(:,k) = 10*log10( (pattern_ang(:,k) + eps_db) / maxval );
    end
end

% Plot angle patterns
figure('Name','MIMO Perceived Beam Pattern vs Angle','NumberTitle','off');
hold on; grid on;
colors = lines(K);
for k = 1:K
    plot(ang_deg, pattern_ang_db(:,k), 'Color', colors(k,:));
end
xlabel('Angle (deg)');
ylabel('Normalized Gain (dB)');
title('MIMO Perceived Beam Pattern vs Angle (per-user normalized)');
legend(arrayfun(@(i) sprintf('User %d',i), 1:K, 'UniformOutput', false), 'Location','best');
xlim([-90 90]);
%ylim([-40 0]);    
% Pattern vs Distance (fixed per-user angle) ------
dist_vec = linspace(0, 100, 1000);   
pattern_dist = zeros(length(dist_vec), K);

for id = 1:length(dist_vec)
    dtest = dist_vec(id);
    for k = 1:K
        theta_k = para.theta_users(k);   % fixed angle for user k

        
        Htest = gen_test_channel(para, dtest, theta_k);   % Mr x Mt

        
        Wk = W_BD(:, (k-1)*Qk + 1 : k*Qk);               % Mt x Qk

        % Effective MIMO matrix
        Heff = Htest * Wk;

        % % MIMO perceived gain (max singular value squared)
        % s = (svd(Heff));
        % if isempty(s)
        %     g = 0;
        % else
        %     g = (s(1))^2;
        % end

        if para.k_abs == 0
            PL_scalar = (1 / dtest)^2;    
        else
            PL_scalar = (1 / dtest)^2 * exp(-para.k_abs * dtest);  % include absorption
        end

        Heff_nopl = Heff / sqrt(PL_scalar);   
        s = sum((Heff_nopl));
        g = (s(1))^2;
        pattern_dist(id, k) = g;
    end
end

%Normalize per-user and convert to dB

pattern_dist_db = zeros(size(pattern_dist));

for k = 1:K
    maxval = max(pattern_dist(:,k));               % strongest focusing distance
    pattern_dist_db(:,k) = 10*log10(pattern_dist(:,k) / maxval);
end


% Plot distance patterns
figure('Name','MIMO Perceived Beam Pattern vs Distance','NumberTitle','off');
hold on; grid on;
for k = 1:K
    plot(dist_vec, pattern_dist_db(:,k), 'Color', colors(k,:));
end
xlabel('Distance (m)');
ylabel('Normalized Gain (dB)');
title('MIMO Perceived Beam Pattern vs Distance (per-user normalized)');
legend(arrayfun(@(i) sprintf('User %d',i), 1:K, 'UniformOutput', false), 'Location','best');
xlim([min(dist_vec) max(dist_vec)]);
% ylim([-60 0]);   % larger dynamic range for distance plots

% Optional: save figures --------------------------
% saveas(gcf,'MIMO_Beampattern_Distance.png');
% figure(1); saveas(gcf,'MIMO_Beampattern_Angle.png');

fprintf('Done. Plots generated: Angle and Distance (MIMO perceived patterns).\n');

para.r
rad2deg(para.theta_users)
para.k_abs
