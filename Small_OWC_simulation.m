% ==========================================
% Simple OWC model
% ==========================================

clear; clc; close all;

% --- INPUT PARAMETERS ---
rho_w = 1025;           % Water density [kg/m^3]
rho_a = 1.225;          % Air density [kg/m^3]
g = 9.81;               % Gravity [m/s^2]

%Waves
H = 1.5;                % Wave height [m]
T = 4;                  % Wave period [s]

%Chamber and turbine parameters
Ac_p = 10;                % Chamber cross-section area [m]normalized by 1 meter of 
At_p = 1;                 % Turbine area [m]
Cp = 0.4;               % Power coefficient (typical Wells turbine)


t = linspace(0, 5*T, 1000);  % Simulate 3 wave periods

% --- WAVE CHARACTERISTICS ---
a = H/2;                % Amplitude [m]
omega = 2*pi/T;         % Angular frequency [rad/s]

% --- WATER SURFACE VELOCITY (sinusoidal) ---
v_w = abs(a * omega * cos(omega * t));

% --- AIR VELOCITY through turbine (From mass balance) ---
v_a = v_w * (Ac_p/At_p);

% --- INSTANTANEOUS TURBINE POWER ---
P_turb = 0.5 * rho_a * At_p * Cp .* (v_a).^3;   % [W]

% --- MEAN POWER over one period ---
P_mean = mean(P_turb);
fprintf('Mean turbine power = %.2f W\n', P_mean);

% --- WAVE POWER per meter of crest width (theoretical) ---
P_wave = (rho_w * g^2 * H^2 * T) / (32 * pi);   % [W/m]
fprintf('Theoritical wave power = %.2f W\n', P_wave)

% --- ESTIMATED EFFICIENCY (dimensionless) ---
eta = (P_mean / (P_wave * At_p))*100;   % efficiency in percentage
fprintf('Approx. OWC efficiency = %.3f\n', eta);

% --- PLOTS ---
figure;
subplot(2,1,1)
plot(t, v_a)
xlabel('Time [s]'); ylabel('Air velocity [m/s]');
title('Air velocity through turbine');

subplot(2,1,2)
plot(t, P_turb/1000)
xlabel('Time [s]'); ylabel('Power [kW]');
title('Instantaneous turbine power');


%% #################################################
% Simple OWC efficiency sweeps (per meter of crest)
clear; clc; close all;

%% Constants (edit as needed)
rho_w = 1025;      % water density [kg/m^3]
rho_a = 1.225;     % air density [kg/m^3]
g     = 9.81;      % gravity [m/s^2]
Cp    = 0.4;       % Wells turbine power coefficient (simplified constant)
Ac_p  = 1.5;       % chamber area per meter of crest [m]  (A_c')

%% Efficiency function per meter (dimensionless)
% eta_pm(H,T,AR) with AR = A_c'/A_t'
eta_pm = @(H,T,AR) ...
    ( 0.5*rho_a*Cp*Ac_p*(4/(3*pi)) .* ((H/2).*2*pi./T).^3 .* (AR.^2) ) ...
  ./ ( (rho_w*g^2.*H.^2.*T)/(32*pi) );

%% 1D sweeps (pick fixed values)
T_fix  = 6;                 % s
H_fix  = 1.5;               % m
AR_fix = 10;                % A_c'/A_t' (dimensionless)

% Ranges
H_vec  = linspace(0.5, 3.0, 60);
T_vec  = linspace(4, 12, 60);
AR_vec = linspace(2, 20, 60);

% Compute
eta_H  = eta_pm(H_vec, T_fix, AR_fix);
eta_T  = eta_pm(H_fix, T_vec, AR_fix);
eta_AR = eta_pm(H_fix, T_fix, AR_vec);

%% Plots: 1D sweeps
figure; plot(H_vec, 100*eta_H, 'LineWidth', 2); grid on;
xlabel('Wave height H [m]'); ylabel('Efficiency \eta [%]');
title(sprintf('\\eta vs H  (T=%.1fs,  AR=A_c''/A_t''=%.1f)', T_fix, AR_fix));

figure; plot(T_vec, 100*eta_T, 'LineWidth', 2); grid on;
xlabel('Wave period T [s]'); ylabel('Efficiency \eta [%]');
title(sprintf('\\eta vs T  (H=%.1fm,  AR=A_c''/A_t''=%.1f)', H_fix, AR_fix));

figure; plot(AR_vec, 100*eta_AR, 'LineWidth', 2); grid on;
xlabel('Area ratio  AR = A_c''/A_t'' [-]'); ylabel('Efficiency \eta [%]');
title(sprintf('\\eta vs AR  (H=%.1fm,  T=%.1fs)', H_fix, T_fix));