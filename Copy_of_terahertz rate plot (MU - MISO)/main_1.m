clc;
clear;
% close all;


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


% A : (Mr*K) x Mt



SNR_dB = -40:8:50;
SNR    = 10.^(SNR_dB/10);

R_BD   = zeros(size(SNR));
R_SLNR = zeros(size(SNR));
R_MRT = zeros(size(SNR));
R_MMSE = zeros(size(SNR));
R_PMMSE = zeros(size(SNR));
R_BD_CMC   = zeros(size(SNR));
R_SLNR_CMC = zeros(size(SNR));
R_MRT_CMC = zeros(size(SNR));
R_MMSE_CMC = zeros(size(SNR));
R_PMMSE_CMC = zeros(size(SNR));

noise0 = para.noise;

for i = 1:length(SNR)

    para.noise = 1 / SNR(i);

    A = gen_channel_thz_mimo(para);

    % Precoding
    W_BD   = beam_BD(para,A);
    W_SLNR = beam_SLNR_reg(para,A);
    W_MRT  = beam_MRT(para,A);
    W_MMSE = beam_MMSE(para,A);
    W_PMMSE = beam_PMMSE(para,A);
    W_BD_CMC   = proj_cmc(W_BD);
    W_SLNR_CMC = proj_cmc(W_SLNR);
    W_MRT_CMC  = proj_cmc(W_MRT);
    W_MMSE_CMC = proj_cmc(W_MMSE);
    W_PMMSE_CMC = proj_cmc(W_PMMSE);

    sumR_BD = 0; sumR_SLNR = 0;
    sumR_MRT = 0; sumR_MMSE = 0;
    sumR_PMMSE = 0; 
    sumR_BD_CMC = 0; sumR_SLNR_CMC = 0;
    sumR_MRT_CMC = 0; sumR_MMSE_CMC = 0;
    sumR_PMMSE_CMC = 0; 

    for k = 1:K

        Hk = A((k-1)*para.M_r+1 : k*para.M_r, :);  % Mr x Mt

        % Desired precoder blocks
        Wk_BD   = W_BD(:,   (k-1)*Qk+1 : k*Qk);
        Wk_SLNR = W_SLNR(:, (k-1)*Qk+1 : k*Qk);
        Wk_MRT  = W_MRT(:,  (k-1)*Qk+1 : k*Qk);
        Wk_MMSE = W_MMSE(:, (k-1)*Qk+1 : k*Qk);
        Wk_PMMSE = W_PMMSE(:, (k-1)*Qk+1 : k*Qk);
        Wk_BD_CMC   = W_BD_CMC(:,   (k-1)*Qk+1 : k*Qk);
        Wk_SLNR_CMC = W_SLNR_CMC(:, (k-1)*Qk+1 : k*Qk);
        Wk_MRT_CMC  = W_MRT_CMC(:,  (k-1)*Qk+1 : k*Qk);
        Wk_MMSE_CMC = W_MMSE_CMC(:, (k-1)*Qk+1 : k*Qk);
        Wk_PMMSE_CMC = W_PMMSE_CMC(:, (k-1)*Qk+1 : k*Qk);

        % ===== Interference covariance =====
        Int_BD   = para.noise * eye(para.M_r);
        Int_SLNR = para.noise * eye(para.M_r);
        Int_MRT  = para.noise * eye(para.M_r);
        Int_MMSE = para.noise * eye(para.M_r);
        Int_PMMSE = para.noise * eye(para.M_r);
        Int_BD_CMC   = para.noise * eye(para.M_r);
        Int_SLNR_CMC = para.noise * eye(para.M_r);
        Int_MRT_CMC  = para.noise * eye(para.M_r);
        Int_MMSE_CMC = para.noise * eye(para.M_r);
        Int_PMMSE_CMC = para.noise * eye(para.M_r);

        for j = 1:K
            if j ~= k
                Hj = Hk;   % same user channel for interference

                Wj_BD   = W_BD(:,   (j-1)*Qk+1 : j*Qk);
                Wj_SLNR = W_SLNR(:, (j-1)*Qk+1 : j*Qk);
                Wj_MRT  = W_MRT(:,  (j-1)*Qk+1 : j*Qk);
                Wj_MMSE = W_MMSE(:, (j-1)*Qk+1 : j*Qk);
                Wj_PMMSE = W_PMMSE(:, (j-1)*Qk+1 : j*Qk);
                Wj_BD_CMC   = W_BD_CMC(:,   (j-1)*Qk+1 : j*Qk);
                Wj_SLNR_CMC = W_SLNR_CMC(:, (j-1)*Qk+1 : j*Qk);
                Wj_MRT_CMC  = W_MRT_CMC(:,  (j-1)*Qk+1 : j*Qk);
                Wj_MMSE_CMC = W_MMSE_CMC(:, (j-1)*Qk+1 : j*Qk);
                Wj_PMMSE_CMC = W_PMMSE_CMC(:, (j-1)*Qk+1 : j*Qk);

                Int_BD   = Int_BD   + Hj*Wj_BD*Wj_BD'*Hj';
                Int_SLNR = Int_SLNR + Hj*Wj_SLNR*Wj_SLNR'*Hj';
                Int_MRT  = Int_MRT  + Hj*Wj_MRT*Wj_MRT'*Hj';
                Int_MMSE = Int_MMSE + Hj*Wj_MMSE*Wj_MMSE'*Hj';
                Int_PMMSE = Int_PMMSE + Hj*Wj_PMMSE*Wj_PMMSE'*Hj';
                Int_BD_CMC   = Int_BD_CMC   + Hj*Wj_BD_CMC*Wj_BD_CMC'*Hj';
                Int_SLNR_CMC = Int_SLNR_CMC + Hj*Wj_SLNR_CMC*Wj_SLNR_CMC'*Hj';
                Int_MRT_CMC  = Int_MRT_CMC  + Hj*Wj_MRT_CMC*Wj_MRT_CMC'*Hj';
                Int_MMSE_CMC = Int_MMSE_CMC + Hj*Wj_MMSE_CMC*Wj_MMSE_CMC'*Hj';
                Int_PMMSE_CMC = Int_PMMSE_CMC + Hj*Wj_PMMSE_CMC*Wj_PMMSE_CMC'*Hj';
            end
        end

        % ===== Rate computation (log-det) =====
        sumR_BD = sumR_BD + real(log2(det( ...
            eye(Qk) + (Hk*Wk_BD)'/Int_BD*(Hk*Wk_BD) )));

        sumR_SLNR = sumR_SLNR + real(log2(det( ...
            eye(Qk) + (Hk*Wk_SLNR)'/Int_SLNR*(Hk*Wk_SLNR) )));

        sumR_MRT = sumR_MRT + real(log2(det( ...
            eye(Qk) + (Hk*Wk_MRT)'/Int_MRT*(Hk*Wk_MRT) )));

        sumR_MMSE = sumR_MMSE + real(log2(det( ...
            eye(Qk) + (Hk*Wk_MMSE)'/Int_MMSE*(Hk*Wk_MMSE) )));

        sumR_PMMSE = sumR_PMMSE + real(log2(det( ...
            eye(Qk) + (Hk*Wk_PMMSE)'/Int_PMMSE*(Hk*Wk_PMMSE) )));
        % sumR_DPC = sumR_DPC + real(log2(det(eye(Qk) + norm(Hk)^2/para.noise)));
        % ===== Rate computation (log-det) =====
        sumR_BD_CMC = sumR_BD_CMC + real(log2(det( ...
            eye(Qk) + (Hk*Wk_BD_CMC)'/Int_BD_CMC*(Hk*Wk_BD_CMC) )));

        sumR_SLNR_CMC = sumR_SLNR_CMC + real(log2(det( ...
            eye(Qk) + (Hk*Wk_SLNR_CMC)'/Int_SLNR_CMC*(Hk*Wk_SLNR_CMC) )));

        sumR_MRT_CMC = sumR_MRT_CMC + real(log2(det( ...
            eye(Qk) + (Hk*Wk_MRT_CMC)'/Int_MRT_CMC*(Hk*Wk_MRT_CMC) )));

        sumR_MMSE_CMC = sumR_MMSE_CMC + real(log2(det( ...
            eye(Qk) + (Hk*Wk_MMSE_CMC)'/Int_MMSE_CMC*(Hk*Wk_MMSE_CMC) )));

        sumR_PMMSE_CMC = sumR_PMMSE_CMC + real(log2(det( ...
            eye(Qk) + (Hk*Wk_PMMSE_CMC)'/Int_PMMSE_CMC*(Hk*Wk_PMMSE_CMC) )));
        % sumR_DPC = sumR_DPC + real(log2(det(eye(Qk) + norm(Hk)^2/para.noise)));
    end

    R_BD(i)   = sumR_BD;
    R_SLNR(i) = sumR_SLNR;
    R_MRT(i)  = sumR_MRT;
    R_MMSE(i) = sumR_MMSE;
    R_PMMSE(i) = sumR_PMMSE;
    % R_DPC(i) = sumR_DPC;
    R_BD_CMC(i)   = sumR_BD_CMC;
    R_SLNR_CMC(i) = sumR_SLNR_CMC;
    R_MRT_CMC(i)  = sumR_MRT_CMC;
    R_MMSE_CMC(i) = sumR_MMSE_CMC;
    R_PMMSE_CMC(i) = sumR_PMMSE_CMC;
    % R_DPC_CMC(i) = sumR_DPC_CMC;
end

para.noise = noise0;

figure;
plot(SNR_dB, real(R_BD),   '-o', 'LineWidth',1.2,Color='#000000');
hold on;
plot(SNR_dB, real(R_SLNR), '-s', 'LineWidth',1.2,Color='r');
plot(SNR_dB, real(R_MRT), '-d', 'LineWidth',1.2,Color='b');
plot(SNR_dB, real(R_MMSE), '-x', 'LineWidth',1,Color='g');
% plot(SNR_dB, real(R_PMMSE), '-+', 'LineWidth',1);
% Rplot(SNR_dB, real(R_DPC), ':', 'LineWidth',1);
plot(SNR_dB, real(R_BD_CMC),   '--', 'LineWidth',1.2,Color='#000000');
plot(SNR_dB, real(R_SLNR_CMC), '--', 'LineWidth',1.2,Color='r');
plot(SNR_dB, real(R_MRT_CMC), '--', 'LineWidth',1.2,Color='b');
plot(SNR_dB, real(R_MMSE_CMC), '--', 'LineWidth',1.2,Color='g');
% plot(SNR_dB, real(R_PMMSE_CMC), '--', 'LineWidth',1);
% plot(SNR_dB, real(R_DPC), '--', 'LineWidth',1);

grid on;
xlabel('SNR (dB)');
ylabel('Sum Rate (bps/Hz)');
legend('BD','Reg-SLNR','MRT','MMSE','BD CMC','Reg-SLNR CMC','MRT CMC','MMSE CMC','Location','northwest');
title('THz MU-MIMO: Sum Rate vs SNR');





