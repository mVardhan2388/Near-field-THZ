clc; clear; close all;

para = para_init();
para.K = 1;
K  = para.K;
Qk = para.Qk;          % MULTI-STREAM (DO NOT SET TO 1)

%% =========================================================
% 1) RECEIVED POWER vs SNR
%% =========================================================

SNR_dB = -10:5:30;
Prx_BD   = zeros(size(SNR_dB));
Prx_SLNR = zeros(size(SNR_dB));
Prx_MRT  = zeros(size(SNR_dB));
Prx_MMSE = zeros(size(SNR_dB));

A = gen_channel_thz_mimo(para);   % FIX channel for fair SNR comparison

for i = 1:length(SNR_dB)

    para.Pt    = 10^(SNR_dB(i)/10);
    para.noise = 1;

    W_BD   = beam_BD(para,A);
    W_SLNR = beam_SLNR_reg(para,A);
    W_MRT  = beam_MRT(para,A);
    W_MMSE = beam_MMSE(para,A);

    for k = 1:K
        Hk = A((k-1)*para.M_r+1 : k*para.M_r, :);

        Wk_BD   = W_BD(:,   (k-1)*Qk+1 : k*Qk);
        Wk_SLNR = W_SLNR(:, (k-1)*Qk+1 : k*Qk);
        Wk_MRT  = W_MRT(:,  (k-1)*Qk+1 : k*Qk);
        Wk_MMSE = W_MMSE(:, (k-1)*Qk+1 : k*Qk);

        Prx_BD(i)   = Prx_BD(i)   + real(trace(Hk*Wk_BD*Wk_BD'*Hk'));
        Prx_SLNR(i) = Prx_SLNR(i) + real(trace(Hk*Wk_SLNR*Wk_SLNR'*Hk'));
        Prx_MRT(i)  = Prx_MRT(i)  + real(trace(Hk*Wk_MRT*Wk_MRT'*Hk'));
        Prx_MMSE(i) = Prx_MMSE(i) + real(trace(Hk*Wk_MMSE*Wk_MMSE'*Hk'));
    end
end

figure;
plot(SNR_dB,10*log10(Prx_BD),'-o'); hold on;
plot(SNR_dB,10*log10(Prx_SLNR),'-s');
plot(SNR_dB,10*log10(Prx_MRT),'-d');
% plot(SNR_dB,10*log10(Prx_MMSE),'-x');
grid on;
xlabel('SNR (dB)');
ylabel('Received Power (dB)');
legend('BD','Reg-SLNR','MRT','MMSE','Location','northwest');
title('Received Power vs SNR (THz Near-field MU-MIMO)');

%% =========================================================
% 2) RECEIVED POWER vs Rx ARRAY SIZE
%% =========================================================

Mr_list = 8:16:128;

Prx_BD   = zeros(size(Mr_list));
Prx_SLNR = zeros(size(Mr_list));
Prx_MRT  = zeros(size(Mr_list));
Prx_MMSE = zeros(size(Mr_list));

for m = 1:length(Mr_list)

    para.M_r = Mr_list(m);
    para.Qk = Mr_list(m);
    Qk = para.Qk;
    para.y_bs_r = ((0:para.M_r-1)-(para.M_r-1)/2)*para.d;

    A = gen_channel_thz_mimo(para);

    W_BD   = beam_BD(para,A);
    W_SLNR = beam_SLNR_reg(para,A);
    W_MRT  = beam_MRT(para,A);
    W_MMSE = beam_MMSE(para,A);

    for k = 1:K
        Hk = A((k-1)*para.M_r+1 : k*para.M_r, :);

        Wk_BD   = W_BD(:,   (k-1)*Qk+1 : k*Qk);
        Wk_SLNR = W_SLNR(:, (k-1)*Qk+1 : k*Qk);
        Wk_MRT  = W_MRT(:,  (k-1)*Qk+1 : k*Qk);
        Wk_MMSE = W_MMSE(:, (k-1)*Qk+1 : k*Qk);

        Prx_BD(m)   = Prx_BD(m)   + real(trace(Hk*Wk_BD*Wk_BD'*Hk'));
        Prx_SLNR(m) = Prx_SLNR(m) + real(trace(Hk*Wk_SLNR*Wk_SLNR'*Hk'));
        Prx_MRT(m)  = Prx_MRT(m)  + real(trace(Hk*Wk_MRT*Wk_MRT'*Hk'));
        Prx_MMSE(m) = Prx_MMSE(m) + real(trace(Hk*Wk_MMSE*Wk_MMSE'*Hk'));
    end
end

figure;
plot(Mr_list,10*log10(Prx_BD),'-o'); hold on;
plot(Mr_list,10*log10(Prx_SLNR),'-s');
plot(Mr_list,10*log10(Prx_MRT),'-d');
% plot(Mr_list,10*log10(Prx_MMSE),'-x');
grid on;
xlabel('Number of Rx Antennas (M_r)');
ylabel('Received Power (dB)');
legend('BD','Reg-SLNR','MRT','MMSE','Location','northwest');
title('Received Power vs Rx Array Size');

%% =========================================================
% 3) RECEIVED POWER vs USER DISTANCE
%% =========================================================

dist_vec = 2:2:40;

Prx_BD   = zeros(size(dist_vec));
Prx_SLNR = zeros(size(dist_vec));
Prx_MRT  = zeros(size(dist_vec));
Prx_MMSE = zeros(size(dist_vec));

para.M_r = 32;
para.Qk = 32;
Qk = para.Qk;

for i = 1:length(dist_vec)

    para.r = dist_vec(i)*ones(1,K);
    A = gen_channel_thz_mimo(para);

    W_BD   = beam_BD(para,A);
    W_SLNR = beam_SLNR_reg(para,A);
    W_MRT  = beam_MRT(para,A);
    W_MMSE = beam_MMSE(para,A);

    for k = 1:K
        Hk = A((k-1)*para.M_r+1 : k*para.M_r, :);

        Wk_BD   = W_BD(:,   (k-1)*Qk+1 : k*Qk);
        Wk_SLNR = W_SLNR(:, (k-1)*Qk+1 : k*Qk);
        Wk_MRT  = W_MRT(:,  (k-1)*Qk+1 : k*Qk);
        Wk_MMSE = W_MMSE(:, (k-1)*Qk+1 : k*Qk);

        Prx_BD(i)   = Prx_BD(i)   + real(trace(Hk*Wk_BD*Wk_BD'*Hk'));
        Prx_SLNR(i) = Prx_SLNR(i) + real(trace(Hk*Wk_SLNR*Wk_SLNR'*Hk'));
        Prx_MRT(i)  = Prx_MRT(i)  + real(trace(Hk*Wk_MRT*Wk_MRT'*Hk'));
        Prx_MMSE(i) = Prx_MMSE(i) + real(trace(Hk*Wk_MMSE*Wk_MMSE'*Hk'));
    end
end

figure;
plot(dist_vec,10*log10(Prx_BD),'-o'); hold on;
plot(dist_vec,10*log10(Prx_SLNR),'-s');
plot(dist_vec,10*log10(Prx_MRT),'-d');
% plot(dist_vec,10*log10(Prx_MMSE),'-x');
grid on;
xlabel('User Distance (m)');
ylabel('Received Power (dB)');
legend('BD','Reg-SLNR','MRT','MMSE','Location','northeast');
title('Received Power vs User Distance (Near-field THz)');
