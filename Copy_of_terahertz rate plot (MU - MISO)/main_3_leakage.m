

clc; clear; close all;

para = para_init();

% para.Qk = 2;                 
Qk = para.Qk;

Mt = para.M_t;
K  = para.K;

%% =========================================================
% 1) LEAKAGE vs Rx ARRAY SIZE
%% =========================================================

Mr_list = 8:32:256;

Leak_BD   = zeros(size(Mr_list));
Leak_SLNR = zeros(size(Mr_list));
Leak_MRT  = zeros(size(Mr_list));
Leak_MMSE = zeros(size(Mr_list));
Leak_PMMSE = zeros(size(Mr_list));
Leak_BD_CMC   = zeros(size(Mr_list));
Leak_SLNR_CMC = zeros(size(Mr_list));
Leak_MRT_CMC  = zeros(size(Mr_list));
Leak_MMSE_CMC = zeros(size(Mr_list));
Leak_PMMSE_CMC = zeros(size(Mr_list));

for m = 1:length(Mr_list)

    para.M_r = Mr_list(m);
    para.Qk = Mr_list(m);
    Qk = para.Qk;
    para.y_bs_r = ((0:para.M_r-1) - (para.M_r-1)/2) * para.d;

    A = gen_channel_thz_mimo(para);

    W_BD   = beam_BD(para, A);          
    W_SLNR = beam_SLNR_reg(para, A);    
    W_MRT  = beam_MRT(para, A);
    W_MMSE = beam_MMSE(para, A);
    W_PMMSE = beam_PMMSE(para, A);
    W_BD_CMC   = proj_cmc(W_BD);          
    W_SLNR_CMC = proj_cmc(W_SLNR);    
    W_MRT_CMC  = proj_cmc(W_MRT);
    W_MMSE_CMC = proj_cmc(W_MMSE);
    W_PMMSE_CMC = proj_cmc(W_PMMSE);

    leak_bd = 0; leak_slnr = 0;
    leak_mrt = 0; leak_mmse = 0;
    leak_pmmse = 0;
    leak_bd_CMC = 0; leak_slnr_CMC = 0;
    leak_mrt_CMC = 0; leak_mmse_CMC = 0;
    leak_pmmse_CMC = 0;

    for k = 1:K

        
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

        for j = 1:K
            if j ~= k
                Hj = A((j-1)*para.M_r+1 : j*para.M_r, :);

                
                leak_bd   = leak_bd   + norm(Hj*Wk_BD,   'fro')^2;
                leak_slnr = leak_slnr + norm(Hj*Wk_SLNR, 'fro')^2;
                leak_mrt  = leak_mrt  + norm(Hj*Wk_MRT,  'fro')^2;
                leak_mmse = leak_mmse + norm(Hj*Wk_MMSE, 'fro')^2;
                leak_pmmse = leak_pmmse + norm(Hj*Wk_PMMSE, 'fro')^2;
                leak_bd_CMC   = leak_bd_CMC   + norm(Hj*Wk_BD_CMC,   'fro')^2;
                leak_slnr_CMC = leak_slnr_CMC + norm(Hj*Wk_SLNR_CMC, 'fro')^2;
                leak_mrt_CMC  = leak_mrt_CMC  + norm(Hj*Wk_MRT_CMC,  'fro')^2;
                leak_mmse_CMC = leak_mmse_CMC + norm(Hj*Wk_MMSE_CMC, 'fro')^2;
                leak_pmmse_CMC = leak_pmmse_CMC + norm(Hj*Wk_PMMSE_CMC, 'fro')^2;
            end
        end
    end

    Leak_BD(m)   = leak_bd;
    Leak_SLNR(m) = leak_slnr;
    Leak_MRT(m)  = leak_mrt;
    Leak_MMSE(m) = leak_mmse;
    Leak_PMMSE(m) = leak_pmmse;
    Leak_BD_CMC(m)   = leak_bd_CMC;
    Leak_SLNR_CMC(m) = leak_slnr_CMC;
    Leak_MRT_CMC(m)  = leak_mrt_CMC;
    Leak_MMSE_CMC(m) = leak_mmse_CMC;
    Leak_PMMSE_CMC(m) = leak_pmmse_CMC;
end

figure;
plot(Mr_list,10*log10(Leak_BD),'-o','LineWidth',1.2,Color='#000000'); hold on;
plot(Mr_list,10*log10(Leak_SLNR),'-s','LineWidth',1.2,Color='r');
plot(Mr_list,10*log10(Leak_MRT),'-d','LineWidth',1.2,Color='b');
plot(Mr_list,10*log10(Leak_MMSE),'-x','LineWidth',1.2,Color='g');
% plot(Mr_list,10*log10(Leak_PMMSE),'-+','LineWidth',1);
plot(Mr_list,10*log10(Leak_BD_CMC),'--','LineWidth',1.2,Color='#000000'); 
plot(Mr_list,10*log10(Leak_SLNR_CMC),'--','LineWidth',1.2,Color='r');
plot(Mr_list,10*log10(Leak_MRT_CMC),'--','LineWidth',1.2,Color='b');
plot(Mr_list,10*log10(Leak_MMSE_CMC),'--','LineWidth',1.2,Color='g');
% plot(Mr_list,10*log10(Leak_PMMSE_CMC),'--','LineWidth',1);
grid on;
xlabel('Number of Rx Antennas (M_r)');
ylabel('Leakage Power (dB)');
legend('BD','Reg-SLNR','MRT','MMSE','BD_CMC','Reg-SLNR_CMC','MRT_CMC','MMSE_CMC','Location','best');
title('Near-Field Leakage vs Rx Array Size');

%% =========================================================
% 2) LEAKAGE vs DISTANCE
%% =========================================================

dist_vec = 2:2:40;
para.M_r = 32;
para.Qk = 32;
Qk = para.Qk;

Leak_BD_d   = zeros(size(dist_vec));
Leak_SLNR_d = zeros(size(dist_vec));
Leak_MRT_d  = zeros(size(dist_vec));
Leak_MMSE_d = zeros(size(dist_vec));
Leak_PMMSE_d = zeros(size(dist_vec));
Leak_BD_d_CMC   = zeros(size(dist_vec));
Leak_SLNR_d_CMC = zeros(size(dist_vec));
Leak_MRT_d_CMC  = zeros(size(dist_vec));
Leak_MMSE_d_CMC = zeros(size(dist_vec));
Leak_PMMSE_d_CMC = zeros(size(dist_vec));

for i = 1:length(dist_vec)

    para.r = dist_vec(i)*ones(1,K);
    A = gen_channel_thz_mimo(para);

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

    leak_bd = 0; leak_slnr = 0;
    leak_mrt = 0; leak_mmse = 0;
    leak_pmmse = 0;
    leak_bd_CMC = 0; leak_slnr_CMC = 0;
    leak_mrt_CMC = 0; leak_mmse_CMC = 0;
    leak_pmmse_CMC = 0;

    for k = 1:K

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

        for j = 1:K
            if j ~= k
                Hj = A((j-1)*para.M_r+1 : j*para.M_r, :);

                leak_bd   = leak_bd   + norm(Hj*Wk_BD,   'fro')^2;
                leak_slnr = leak_slnr + norm(Hj*Wk_SLNR, 'fro')^2;
                leak_mrt  = leak_mrt  + norm(Hj*Wk_MRT,  'fro')^2;
                leak_mmse = leak_mmse + norm(Hj*Wk_MMSE, 'fro')^2;
                leak_pmmse = leak_pmmse + norm(Hj*Wk_PMMSE, 'fro')^2;
                leak_bd_CMC   = leak_bd_CMC   + norm(Hj*Wk_BD_CMC,   'fro')^2;
                leak_slnr_CMC = leak_slnr_CMC + norm(Hj*Wk_SLNR_CMC, 'fro')^2;
                leak_mrt_CMC  = leak_mrt_CMC  + norm(Hj*Wk_MRT_CMC,  'fro')^2;
                leak_mmse_CMC = leak_mmse_CMC + norm(Hj*Wk_MMSE_CMC, 'fro')^2;
                leak_pmmse_CMC = leak_pmmse_CMC + norm(Hj*Wk_PMMSE_CMC, 'fro')^2;

            end
        end
    end

    Leak_BD_d(i)   = leak_bd;
    Leak_SLNR_d(i) = leak_slnr;
    Leak_MRT_d(i)  = leak_mrt;
    Leak_MMSE_d(i) = leak_mmse;
    Leak_PMMSE_d(i) = leak_pmmse;
    Leak_BD_d_CMC(i)   = leak_bd_CMC;
    Leak_SLNR_d_CMC(i) = leak_slnr_CMC;
    Leak_MRT_d_CMC(i)  = leak_mrt_CMC;
    Leak_MMSE_d_CMC(i) = leak_mmse_CMC;
    Leak_PMMSE_d_CMC(i) = leak_pmmse_CMC;

end

figure;
plot(dist_vec,10*log10(Leak_BD_d),'-o','LineWidth',1.2,Color='#000000'); hold on;
plot(dist_vec,10*log10(Leak_SLNR_d),'-s','LineWidth',1.2,Color='r');
plot(dist_vec,10*log10(Leak_MRT_d),'-d','LineWidth',1.2,Color='b');
plot(dist_vec,10*log10(Leak_MMSE_d),'-x','LineWidth',1.2,Color='g');
% plot(dist_vec,10*log10(Leak_PMMSE_d),'-+','LineWidth',1);
plot(dist_vec,10*log10(Leak_BD_d_CMC),'--','LineWidth',1.2,Color='#000000'); 
plot(dist_vec,10*log10(Leak_SLNR_d_CMC),'--','LineWidth',1.2,Color='r');
plot(dist_vec,10*log10(Leak_MRT_d_CMC),'--','LineWidth',1.2,Color='b');
plot(dist_vec,10*log10(Leak_MMSE_d_CMC),'--','LineWidth',1.2,Color='g');
% plot(dist_vec,10*log10(Leak_PMMSE_d),'--','LineWidth',1);
grid on;
xlabel('User Distance (m)');
ylabel('Leakage Power (dB)');
legend('BD','Reg-SLNR','MRT','MMSE','BD CMC','Reg-SLNR CMC','MRT CMC','MMSE CMC','Location','best');
title('Near-Field Leakage vs Distance');
