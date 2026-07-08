function para = para_init()

% ----------------------------------------------------
% Parameter initialization for Near-field THz MU-MIMO
% Includes HITRAN-based molecular absorption
% ----------------------------------------------------

% Physical constants
para.c = 3e8;                     % Speed of light (m/s)
para.f = 300e9;                   % Carrier frequency (300 GHz)
para.lambda = para.c / para.f;


% hitran_path = ...
% "C:\Users\Mallena\IIITG research\Near_field_beamforming\HITRAN.par";
% 
% H = importhitran(hitran_path);
% 
% % Convert carrier frequency to wavenumber (cm^-1)
% fc0 = para.f / para.c;         % Hz / (m/s) = 1/m
% fc0 = fc0 / 100;               % convert 1/m → cm^-1
% 
% BW_cm = 1;                     % ±1 cm^-1 window
% 
% idx = abs(H.transitionWavenumber - fc0) <= BW_cm;
% 
% hitran.fc          = H.transitionWavenumber(idx);
% hitran.S           = H.lineIntensity(idx);
% hitran.alpha_air   = H.airBroadenedWidth(idx);
% hitran.alpha_self  = H.selfBroadenedWidth(idx);
% hitran.gamma       = H.temperatureDependence(idx);
% hitran.delta       = H.pressureShift(idx);



para.c = 3e8;
para.f = 300e9;
para.lambda = para.c/para.f;
%para.k_abs = 0.03;
tbl = [1e11, 0.01; 2e11, 0.1; 3e11, 1.5]; 
%para.k_abs = get_absorption_coef(para.f, 'table', tbl, 'units', 'dB_km');

para.M_t = 2048;
para.M_r = 32;
para.Qk  = 32;
para.K   = 4;

para.d = para.lambda/2;

n = 0:para.M_t-1;
para.y_bs_t = (n - (para.M_t-1)/2) * para.d;

m = 0:para.M_r-1;
para.y_bs_r = (m - (para.M_r-1)/2) * para.d;

% para.r = 0.01 + (1-0.01)*rand(1,para.K);
% para.r = (randi([10 60],1,para.K));
 para.r = [18 25 35 50];
% para.theta_users = deg2rad(randi([-20 20],1,para.K));
para.theta_users = deg2rad([-15 -5 5 15]);

para.Pt = 10^(10/10);
para.Bw = 50e9;
para.tem = 300;
para.noise = 1.38e-23 * para.tem * para.Bw;
% para.k_abs = get_absorption_coef(para.f, para, hitran);
para.k_abs = 4.1073e-22;


end
