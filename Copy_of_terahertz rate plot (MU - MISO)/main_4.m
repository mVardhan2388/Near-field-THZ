clc; clear; 
% close all;

para = para_init();

para.Qk = 2;                 
Qk      = para.Qk;
K       = para.K;

dist_vec = 1:4:41;           % meters

R_BD   = zeros(size(dist_vec));
R_SLNR = zeros(size(dist_vec));
R_MRT  = zeros(size(dist_vec));
R_MMSE = zeros(size(dist_vec));
R_BD_CMC   = zeros(size(dist_vec));
R_SLNR_CMC = zeros(size(dist_vec));
R_MRT_CMC  = zeros(size(dist_vec));
R_MMSE_CMC = zeros(size(dist_vec));

for i = 1:length(dist_vec)

    para.r = dist_vec(i) * ones(1,K);

    A = gen_channel_thz_mimo(para);

    % Precoding
    W_BD   = beam_BD(para,A);
    W_SLNR = beam_SLNR_reg(para,A);
    W_MRT  = beam_MRT(para,A);
    W_MMSE = beam_MMSE(para,A);
    W_BD_CMC   = proj_cmc(W_BD);
    W_SLNR_CMC = proj_cmc(W_SLNR);
    W_MRT_CMC  = proj_cmc(W_MRT);
    W_MMSE_CMC = proj_cmc(W_MMSE);

    sumR_BD = 0; sumR_SLNR = 0;
    sumR_MRT = 0; sumR_MMSE = 0; 
    sumR_BD_CMC = 0; sumR_SLNR_CMC = 0;
    sumR_MRT_CMC = 0; sumR_MMSE_CMC = 0; 

    for k = 1:K

        Hk = A((k-1)*para.M_r+1 : k*para.M_r, :);  % Mr x Mt

        % Desired precoder blocks
        Wk_BD   = W_BD(:,   (k-1)*Qk+1 : k*Qk);
        Wk_SLNR = W_SLNR(:, (k-1)*Qk+1 : k*Qk);
        Wk_MRT  = W_MRT(:,  (k-1)*Qk+1 : k*Qk);
        Wk_MMSE = W_MMSE(:, (k-1)*Qk+1 : k*Qk);
        Wk_BD_CMC   = W_BD_CMC(:,   (k-1)*Qk+1 : k*Qk);
        Wk_SLNR_CMC = W_SLNR_CMC(:, (k-1)*Qk+1 : k*Qk);
        Wk_MRT_CMC  = W_MRT_CMC(:,  (k-1)*Qk+1 : k*Qk);
        Wk_MMSE_CMC = W_MMSE_CMC(:, (k-1)*Qk+1 : k*Qk);

        % ===== Interference covariance =====
        Int_BD   = para.noise * eye(para.M_r);
        Int_SLNR = para.noise * eye(para.M_r);
        Int_MRT  = para.noise * eye(para.M_r);
        Int_MMSE = para.noise * eye(para.M_r);
        Int_BD_CMC   = para.noise * eye(para.M_r);
        Int_SLNR_CMC = para.noise * eye(para.M_r);
        Int_MRT_CMC  = para.noise * eye(para.M_r);
        Int_MMSE_CMC = para.noise * eye(para.M_r);

        for j = 1:K
            if j ~= k
                Hj = Hk;   % same user channel for interference

                Wj_BD   = W_BD(:,   (j-1)*Qk+1 : j*Qk);
                Wj_SLNR = W_SLNR(:, (j-1)*Qk+1 : j*Qk);
                Wj_MRT  = W_MRT(:,  (j-1)*Qk+1 : j*Qk);
                Wj_MMSE = W_MMSE(:, (j-1)*Qk+1 : j*Qk);
                Wj_BD_CMC   = W_BD_CMC(:,   (j-1)*Qk+1 : j*Qk);
                Wj_SLNR_CMC = W_SLNR_CMC(:, (j-1)*Qk+1 : j*Qk);
                Wj_MRT_CMC  = W_MRT_CMC(:,  (j-1)*Qk+1 : j*Qk);
                Wj_MMSE_CMC = W_MMSE_CMC(:, (j-1)*Qk+1 : j*Qk);


                Int_BD   = Int_BD   + Hj*Wj_BD*Wj_BD'*Hj';
                Int_SLNR = Int_SLNR + Hj*Wj_SLNR*Wj_SLNR'*Hj';
                Int_MRT  = Int_MRT  + Hj*Wj_MRT*Wj_MRT'*Hj';
                Int_MMSE = Int_MMSE + Hj*Wj_MMSE*Wj_MMSE'*Hj';
                Int_BD_CMC   = Int_BD_CMC   + Hj*Wj_BD_CMC*Wj_BD_CMC'*Hj';
                Int_SLNR_CMC = Int_SLNR_CMC + Hj*Wj_SLNR_CMC*Wj_SLNR_CMC'*Hj';
                Int_MRT_CMC  = Int_MRT_CMC  + Hj*Wj_MRT_CMC*Wj_MRT_CMC'*Hj';
                Int_MMSE_CMC = Int_MMSE_CMC + Hj*Wj_MMSE_CMC*Wj_MMSE_CMC'*Hj';
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
        sumR_BD_CMC = sumR_BD_CMC + real(log2(det( ...
            eye(Qk) + (Hk*Wk_BD_CMC)'/Int_BD_CMC*(Hk*Wk_BD_CMC) )));

        sumR_SLNR_CMC = sumR_SLNR_CMC + real(log2(det( ...
            eye(Qk) + (Hk*Wk_SLNR_CMC)'/Int_SLNR_CMC*(Hk*Wk_SLNR_CMC) )));

        sumR_MRT_CMC = sumR_MRT_CMC + real(log2(det( ...
            eye(Qk) + (Hk*Wk_MRT_CMC)'/Int_MRT_CMC*(Hk*Wk_MRT_CMC) )));

        sumR_MMSE_CMC = sumR_MMSE_CMC + real(log2(det( ...
            eye(Qk) + (Hk*Wk_MMSE_CMC)'/Int_MMSE_CMC*(Hk*Wk_MMSE_CMC) )));
    end

    R_BD(i)   = sumR_BD;
    R_SLNR(i) = sumR_SLNR;
    R_MRT(i)  = sumR_MRT;
    R_MMSE(i) = sumR_MMSE;
    R_BD_CMC(i)   = sumR_BD_CMC;
    R_SLNR_CMC(i) = sumR_SLNR_CMC;
    R_MRT_CMC(i)  = sumR_MRT_CMC;
    R_MMSE_CMC(i) = sumR_MMSE_CMC;
end

%% ===================== PLOT =====================
figure;
plot(dist_vec,R_BD/4,'-o','LineWidth',1.2,Color='#000000'); hold on;
plot(dist_vec,R_SLNR/4,'-s','LineWidth',1.2,Color='r');
plot(dist_vec,R_MRT/4,'-d','LineWidth',1.2,Color='b');
plot(dist_vec,R_MMSE/4,'-x','LineWidth',1.2,Color='g');
plot(dist_vec,R_BD_CMC/4,'--','LineWidth',1.2,Color='#000000'); 
plot(dist_vec,R_SLNR_CMC/4,'--','LineWidth',1.2,Color='r');
plot(dist_vec,R_MRT_CMC/4,'--','LineWidth',1.2,Color='b');
plot(dist_vec,R_MMSE_CMC/4,'--','LineWidth',1.2,Color='g');
grid on;
xlabel('User distance (m)');
ylabel('Sum Rate (bps/Hz)');
legend('BD','Reg-SLNR','MRT','MMSE','BD CMC','Reg-SLNR CMC','MRT CMC','MMSE CMC','Location','northeast');
title('THz Near-field MU-MIMO: Sum Rate vs Distance (Multi-stream)');
xlim([1,41]);
