clc; 
clear; 
%close all;

% PARAMETERS 
para = para_init();

para.Qk = 1;                 % SINGLE STREAM PER USER (IMPORTANT)
Qk      = 1;

Mt = para.M_t;
K  = para.K;


%          RATE & EE vs Rx ARRAY SIZE

Mr_list = 8:16:256;

R_BD   = zeros(size(Mr_list));
R_SLNR = zeros(size(Mr_list));
R_MRT  = zeros(size(Mr_list));
R_MMSE = zeros(size(Mr_list));
R_PMMSE = zeros(size(Mr_list));

EE_BD   = zeros(size(Mr_list));
EE_SLNR = zeros(size(Mr_list));
EE_MRT  = zeros(size(Mr_list));
EE_MMSE = zeros(size(Mr_list));
EE_PMMSE = zeros(size(Mr_list));

for m = 1:length(Mr_list)

    para.M_r = Mr_list(m);
    para.y_bs_r = ((0:para.M_r-1) - (para.M_r-1)/2) * para.d;

    A = gen_channel_thz_mimo(para);

    % Precoding matrices
    W_BD   = beam_BD(para, A);          
    W_SLNR = beam_SLNR_reg(para, A);    
    W_MRT  = beam_MRT(para, A);
    W_MMSE = beam_MMSE(para, A);
    W_PMMSE = beam_PMMSE(para, A);

    sumR_BD = 0; sumR_SLNR = 0;
    sumR_MRT = 0; sumR_MMSE = 0;
    sumR_PMMSE = 0;

    for k = 1:K
        Ak = A((k-1)*para.M_r+1 : k*para.M_r, :);

        wk_BD   = W_BD(:,k);
        wk_SLNR = W_SLNR(:,k);
        wk_MRT  = W_MRT(:,k);
        wk_MMSE = W_MMSE(:,k);
        wk_PMMSE = W_PMMSE(:,k);

        sumR_BD   = sumR_BD   + rate_cal(para, Ak, wk_BD);
        sumR_SLNR = sumR_SLNR + rate_cal(para, Ak, wk_SLNR);
        sumR_MRT  = sumR_MRT  + rate_cal(para, Ak, wk_MRT);
        sumR_MMSE = sumR_MMSE + rate_cal(para, Ak, wk_MMSE);
        sumR_PMMSE = sumR_PMMSE + rate_cal(para, Ak, wk_PMMSE);
    end

    R_BD(m)   = real(sumR_BD);
    R_SLNR(m) = real(sumR_SLNR);
    R_MRT(m)  = real(sumR_MRT);
    R_MMSE(m) = real(sumR_MMSE);
    R_PMMSE(m) = real(sumR_PMMSE);

    EE_BD(m)   = R_BD(m)   / para.Pt;
    EE_SLNR(m) = R_SLNR(m) / para.Pt;
    EE_MRT(m)  = R_MRT(m)  / para.Pt;
    EE_MMSE(m) = R_MMSE(m) / para.Pt;
    EE_PMMSE(m) = R_PMMSE(m) / para.Pt;
end

figure;
plot(Mr_list,R_BD,'-o','LineWidth',1); hold on;
plot(Mr_list,R_SLNR,'-s','LineWidth',1);
plot(Mr_list,R_MRT,'-d','LineWidth',1);
plot(Mr_list,R_MMSE,'-x','LineWidth',1);
plot(Mr_list,R_PMMSE,'-+','LineWidth',1);
grid on;
xlabel('Number of Rx Antennas (M_r)');
ylabel('Sum Rate (bps/Hz)');
legend('BD','Reg-SLNR','MRT','MMSE','PMMSE','Location','northwest');
title('THz MU-MIMO: Sum Rate vs Rx Array Size');

figure;
plot(Mr_list,EE_BD,'-o','LineWidth',1); hold on;
plot(Mr_list,EE_SLNR,'-s','LineWidth',1);
plot(Mr_list,EE_MRT,'-d','LineWidth',1);
plot(Mr_list,EE_MMSE,'-x','LineWidth',1);
plot(Mr_list,EE_PMMSE,'-+','LineWidth',1);
grid on;
xlabel('Number of Rx Antennas (M_r)');
ylabel('Energy Efficiency (bits/J/Hz)');
legend('BD','Reg-SLNR','MRT','MMSE','PMMSE','Location','northwest');
title('THz MU-MIMO: Energy Efficiency vs Rx Array Size');
