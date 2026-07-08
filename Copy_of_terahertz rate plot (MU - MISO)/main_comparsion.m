clc; 
clear; 
close all;


para = para_init();

Mt = para.M_t;
Mr = para.M_r;
K  = para.K;
Qk = para.Qk;

user_labels = arrayfun(@(i) ...
    sprintf('User %d (r=%dm, \\theta=%d°)', ...
    i, round(para.r(i)), round(rad2deg(para.theta_users(i))) ), ...
    1:K, 'UniformOutput', false);
colors = lines(K);


A = gen_channel_thz_mimo(para);    % (Mr*K) x Mt


W_BD   = beam_BD(para, A);               % Mt x (Qk*K)
W_SLNR = beam_SLNR_reg(para, A);          % Mt x (Qk*K)

%  ANGLE BEAM PATTERN 
ang_deg = linspace(-90,90,721);
ang_rad = deg2rad(ang_deg);

pattern_BD_ang   = zeros(length(ang_deg), K);
pattern_SLNR_ang = zeros(length(ang_deg), K);

for ia = 1:length(ang_rad)
    theta = ang_rad(ia);
    for k = 1:K

        Htest = gen_test_channel(para, para.r(k), theta);   % Mr x Mt

        Wk_BD   = W_BD(:,   (k-1)*Qk+1 : k*Qk);
        Wk_SLNR = W_SLNR(:, (k-1)*Qk+1 : k*Qk);

        gBD   = svds(Htest * Wk_BD,   1)^2;
        gSLNR = svds(Htest * Wk_SLNR, 1)^2;

        pattern_BD_ang(ia,k)   = gBD;
        pattern_SLNR_ang(ia,k) = gSLNR;
    end
end

% pattern_BD_ang_db   = 10*log10(pattern_BD_ang ./ max(pattern_BD_ang(:)));
% pattern_SLNR_ang_db = 10*log10(pattern_SLNR_ang ./ max(pattern_SLNR_ang(:)));

% Normalize per-user and convert to dB
pattern_BD_ang_db = zeros(size(pattern_BD_ang));
eps_db = 1e-12;
for k = 1:K
    maxval = max(pattern_BD_ang(:,k));
    if maxval <= 0
        pattern_BD_ang_db(:,k) = -140 * ones(size(pattern_BD_ang(:,k)));
    else
        pattern_BD_ang_db(:,k) = 10*log10( (pattern_BD_ang(:,k) + eps_db) / maxval );
    end
end

pattern_SLNR_ang_db = zeros(size(pattern_SLNR_ang));
eps_db = 1e-12;
for k = 1:K
    maxval = max(pattern_SLNR_ang(:,k));
    if maxval <= 0
        pattern_SLNR_ang_db(:,k) = -140 * ones(size(pattern_SLNR_ang(:,k)));
    else
        pattern_SLNR_ang_db(:,k) = 10*log10( (pattern_SLNR_ang(:,k) + eps_db) / maxval );
    end
end

% PLOT ANGLE 
figure('Name','Angle – BD'); 
hold on; 
grid on;
for k = 1:K
    plot(ang_deg, pattern_BD_ang_db(:,k), 'Color', colors(k,:));
end
xlabel('Angle (deg)'); ylabel('Gain (dB)');
title('THz MU-MIMO Beam Pattern vs Angle (BD)');
legend(user_labels,'Location','best'); 
xlim([-90 90]);

figure('Name','Angle – SLNR'); hold on; grid on;
for k = 1:K
    plot(ang_deg, pattern_SLNR_ang_db(:,k), 'Color', colors(k,:));
end
xlabel('Angle (deg)'); 
ylabel('Gain (dB)');
title('THz MU-MIMO Beam Pattern vs Angle (Regularized SLNR)');
legend(user_labels,'Location','best'); 
xlim([-90 90]);

%  DISTANCE BEAM PATTERN 
dist = linspace(0.01,60,800);

pattern_BD_dist   = zeros(length(dist), K);
pattern_SLNR_dist = zeros(length(dist), K);

for id = 1:length(dist)
    d = dist(id);
    for k = 1:K

        Htest = gen_test_channel(para, d, para.theta_users(k));

        Wk_BD   = W_BD(:,   (k-1)*Qk+1 : k*Qk);
        Wk_SLNR = W_SLNR(:, (k-1)*Qk+1 : k*Qk);

        % gBD   = svds(Htest * Wk_BD,   1)^2;
        % gSLNR = svds(Htest * Wk_SLNR, 1)^2;

        Heff_BD = Htest * Wk_BD;
        Heff_SLNR = Htest * Wk_SLNR;

        if para.k_abs == 0
            PL_scalar = (1 / d)^2;    
        else
            PL_scalar = (1 / d)^2 * exp(-2*para.k_abs * d);  % include absorption
        end

        Heff_nopl = Heff_BD / sqrt(PL_scalar);   
        s = sum((Heff_nopl));
        gBD = (s(1))^2;

        Heff_nopl = Heff_SLNR / sqrt(PL_scalar);   
        s = sum((Heff_nopl));
        gSLNR = (s(1))^2;
        

        pattern_BD_dist(id,k)   = gBD;
        pattern_SLNR_dist(id,k) = gSLNR;
    end
end

% pattern_BD_dist_db   = 10*log10(pattern_BD_dist / max(pattern_BD_dist(:)));
% pattern_SLNR_dist_db = 10*log10(pattern_SLNR_dist / max(pattern_SLNR_dist(:)));

pattern_BD_dist_db = zeros(size(pattern_BD_dist));

for k = 1:K
    maxval = max(pattern_BD_dist(:,k));               % strongest focusing distance
    pattern_BD_dist_db(:,k) = 10*log10(pattern_BD_dist(:,k) / maxval);
end

pattern_SLNR_dist_db = zeros(size(pattern_SLNR_dist));

for k = 1:K
    maxval = max(pattern_SLNR_dist(:,k));               % strongest focusing distance
    pattern_SLNR_dist_db(:,k) = 10*log10(pattern_SLNR_dist(:,k) / maxval);
end
% PLOT DISTANCE 
figure('Name','Distance – BD'); 
hold on; 
grid on;
for k = 1:K
    plot(dist, pattern_BD_dist_db(:,k), 'Color', colors(k,:));
end
xlabel('Distance (m)');
ylabel('Gain (dB)');
title('THz MU-MIMO Beam Pattern vs Distance (BD)');
legend(user_labels,'Location','best');

figure('Name','Distance – SLNR'); 
hold on; 
grid on;
for k = 1:K
    plot(dist, pattern_SLNR_dist_db(:,k), 'Color', colors(k,:));
end
xlabel('Distance (m)'); 
ylabel('Gain (dB)');
title('THz MU-MIMO Beam Pattern vs Distance (Regularized SLNR)');
legend(user_labels,'Location','best');


disp('✔ BD vs Regularized SLNR comparison completed');
disp('User distances (m):'); 
disp(para.r);
disp('User angles (deg):'); 
disp(rad2deg(para.theta_users));


% comparsion SLNR vs BD

figure('Name', 'Distance - BD vs SLNR');
hold on;
grid on;
plot(dist, pattern_BD_dist_db(:,1), 'Color', colors(1,:));
plot(dist, pattern_SLNR_dist_db(:,1), 'Color', colors(2,:));
title('THz MU-MIMO Beam Pattern vs Distance (SLNR vs BD)');


figure('Name', 'Angle - BD vs SLNR');
hold on;
grid on;
plot(ang_deg, pattern_BD_ang_db(:,1), 'Color', colors(1,:));
plot(ang_deg, pattern_SLNR_ang_db(:,1), 'Color', colors(2,:));
title('THz MU-MIMO Beam Pattern vs Angle (SLNR vs BD)');



