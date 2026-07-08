clc;
clear;
close all;

%% ================== PARAMETERS ==================
para = para_init();

% para.k = 2;

Mt = para.M_t;
Mr = para.M_r;
K  = para.K;
Qk = para.Qk;

colors = lines(4);   % BD, SLNR, MRT, MMSE

user_labels = arrayfun(@(i) ...
    sprintf('User %d (r=%dm, \\theta=%d°)', ...
    i, round(para.r(i)), round(rad2deg(para.theta_users(i))) ), ...
    1:K, 'UniformOutput', false);

%% ================== CHANNEL ==================
A = gen_channel_thz_mimo(para);     % (Mr*K) x Mt

%% ================== PRECODERS ==================
disp('Computing precoders...');

W.BD   = beam_BD(para, A);
W.SLNR = beam_SLNR_reg(para, A);
W.MRT  = beam_MRT(para, A);
W.MMSE = beam_MMSE(para, A);
W.BD_CMC   = proj_cmc(W.BD);
W.SLNR_CMC = proj_cmc(W.SLNR);
W.MRT_CMC  = proj_cmc(W.MRT);
W.MMSE_CMC = proj_cmc(W.MMSE);


precoder_names = fieldnames(W);

%% ================== ANGLE GRID ==================
ang_deg = linspace(-90,90,721);
ang_rad = deg2rad(ang_deg);

pattern_ang = struct();

for p = 1:length(precoder_names)
    pname = precoder_names{p};
    pattern_ang.(pname) = zeros(length(ang_rad), K);
end

for ia = 1:length(ang_rad)
    theta = ang_rad(ia);
    for k = 1:K
        Htest = gen_test_channel(para, para.r(k), theta);

        for p = 1:length(precoder_names)
            pname = precoder_names{p};
            Wk = W.(pname)(:, (k-1)*Qk+1:k*Qk);
            g = svds(Htest * Wk, 1)^2;
            pattern_ang.(pname)(ia,k) = g;
        end
    end
end

%% ================== ANGLE NORMALIZATION ==================
pattern_ang_db = struct();
eps_db = 1e-12;

for p = 1:length(precoder_names)
    pname = precoder_names{p};
    X = pattern_ang.(pname);
    Xdb = zeros(size(X));
    for k = 1:K
        Xdb(:,k) = 10*log10((X(:,k)+eps_db)/max(X(:,k)));
    end
    pattern_ang_db.(pname) = Xdb;
end

figure;
hold on; grid on;
for k = para.K:-1:1
    plot(ang_deg, pattern_ang_db.BD(:,k), ...
        'LineWidth',1.2, 'Color', colors(k,:));
end
xlabel('Angle (deg)');
ylabel('Normalized Power (dB)');
xlim([-90 90]);


figure('Name','Angle Beam Comparison');
hold on; grid on;

% ----- Unconstrained -----
plot(ang_deg, pattern_ang_db.BD(:,k),   ...
    'LineWidth',1.2,'Color',colors(1,:));

plot(ang_deg, pattern_ang_db.SLNR(:,k), ...
    'LineWidth',1.2,'Color',colors(2,:));

plot(ang_deg, pattern_ang_db.MRT(:,k),  ...
    'LineWidth',1.2,'Color',colors(3,:));

plot(ang_deg, pattern_ang_db.MMSE(:,k), ...
    'LineWidth',1.2,'Color',colors(4,:));

% ----- CMC (same colors, dashed) -----
plot(ang_deg, pattern_ang_db.BD_CMC(:,k),   ...
    '--','LineWidth',1.2,'Color',colors(1,:));

plot(ang_deg, pattern_ang_db.SLNR_CMC(:,k), ...
    '--','LineWidth',1.2,'Color',colors(2,:));

plot(ang_deg, pattern_ang_db.MRT_CMC(:,k),  ...
    '--','LineWidth',1.2,'Color',colors(3,:));

plot(ang_deg, pattern_ang_db.MMSE_CMC(:,k), ...
    '--','LineWidth',1.2,'Color',colors(4,:));

xlabel('Angle (deg)');
ylabel('Normalized Power (dB)');
xlim([-90 90]);

legend('BD','SLNR','MRT','MMSE', ...
       'BD CMC','SLNR CMC','MRT CMC','MMSE CMC', ...
       'Location','best');


%% ================== DISTANCE GRID ==================
dist = linspace(5,70,1000);
pattern_dist = struct();

for p = 1:length(precoder_names)
    pname = precoder_names{p};
    pattern_dist.(pname) = zeros(length(dist), K);
end

for id = 1:length(dist)
    d = dist(id);
    for k = 1:K
        Htest = gen_test_channel(para, d, para.theta_users(k));

        % if para.k_abs == 0
        %     PL = (1/d)^2;
        % else
        %     PL = (1/d)^2 * exp(-2*para.k_abs*d);
        % end

        for p = 1:length(precoder_names)
            pname = precoder_names{p};
            Wk = W.(pname)(:, (k-1)*Qk+1:k*Qk);
            Heff = Htest * Wk;
            Heff_nopl = Heff;% / sqrt(PL);
            s = sum(Heff_nopl);
            pattern_dist.(pname)(id,k) = abs(s(1))^2;
        end
    end
end

%% ================== DISTANCE NORMALIZATION ==================
pattern_dist_db = struct();

for p = 1:length(precoder_names)
    pname = precoder_names{p};
    X = pattern_dist.(pname);
    Xdb = zeros(size(X));
    for k = 1:K
        Xdb(:,k) = 10*log10(X(:,k)/max(X(:,k)));
    end
    pattern_dist_db.(pname) = Xdb;
end

% DISTANCE PLOT 

figure;
hold on;
grid on;
for k = 1:para.K
    plot(dist, pattern_dist_db.BD(:,k),'LineWidth',1.2,'Color',colors(k,:));
end
xlabel('Distance (m)');
ylabel('Normalized Power (dB)');
xlim([5 70]);

figure('Name','Distance Beam Comparison');
hold on; grid on;

% ----- Unconstrained -----
plot(dist, pattern_dist_db.BD(:,k),   ...
    'LineWidth',1.2,'Color',colors(1,:));

plot(dist, pattern_dist_db.SLNR(:,k), ...
    'LineWidth',1.2,'Color',colors(2,:));

plot(dist, pattern_dist_db.MRT(:,k),  ...
    'LineWidth',1.2,'Color',colors(3,:));

plot(dist, pattern_dist_db.MMSE(:,k), ...
    'LineWidth',1.2,'Color',colors(4,:));

% ----- CMC -----
plot(dist, pattern_dist_db.BD_CMC(:,k),   ...
    '--','LineWidth',1.2,'Color',colors(1,:));

plot(dist, pattern_dist_db.SLNR_CMC(:,k), ...
    '--','LineWidth',1.2,'Color',colors(2,:));

plot(dist, pattern_dist_db.MRT_CMC(:,k),  ...
    '--','LineWidth',1.2,'Color',colors(3,:));

plot(dist, pattern_dist_db.MMSE_CMC(:,k), ...
    '--','LineWidth',1.2,'Color',colors(4,:));

xlabel('Distance (m)');
ylabel('Normalized Power (dB)');
xlim([5 70]);

legend('BD','SLNR','MRT','MMSE', ...
       'BD CMC','SLNR CMC','MRT CMC','MMSE CMC', ...
       'Location','best');


%% ================== SUMMARY ==================
disp('✔ Beam comparator completed');
disp(['User index used for comparison: ', num2str(k)]);
disp('Precorder set: BD, SLNR, MRT, MMSE');
