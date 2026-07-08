function k_abs = get_absorption_coef(f, para, hitran)
% ---------------------------------------------------------
% Molecular absorption coefficient using HITRAN data
% f      : carrier frequency (Hz)
% para   : system parameters
% hitran : extracted HITRAN lines
% Output : k_abs (1/m)
% ---------------------------------------------------------

c  = para.c;
T  = para.tem;          % temperature (K)
P  = 1;                 % pressure (atm) – standard assumption
T0 = 296;               % HITRAN reference temperature (K)

% Frequency → wavenumber
nu = f / c / 100;       % cm^-1

k_abs_cm = 0;

for i = 1:length(hitran.fc)

    nu0 = hitran.fc(i);
    S   = hitran.S(i);
    g_a = hitran.alpha_air(i);
    n   = hitran.gamma(i);

    % Temperature-scaled linewidth
    gamma_L = g_a * (T0 / T)^n * P;

    % Lorentzian line shape
    L = (gamma_L / pi) ./ ...
        ( (nu - nu0).^2 + gamma_L^2 );

    % Accumulate absorption
    k_abs_cm = k_abs_cm + S * L;
end

% Convert cm^-1 → m^-1
k_abs = k_abs_cm * 100;

end
