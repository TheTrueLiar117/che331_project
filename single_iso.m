%% Single Isothermal CSTR Model (Exact Equations S1-S10 & S33-S39)
% Reference: Bremer & Sundmacher (2021), Supplementary Material Tables S1 & S2
clear; clc; close all;

%% 1. Exact Parameters from Table S1 & Table S2
p.R_gas      = 8.314462618;   % Universal gas constant  [J/(mol*K)]

% Table S1: Reactor Reference Parameters
p.T_in       = 300.0;         % Inlet temperature [K]
p.T_op       = 450.0;         % Operating temperature [K]
p.p_in       = 5.0;           % Inlet pressure [bar]
p.rho_cat    = 2355.2;        % Catalyst density (porous material) [kg/m^3]
p.eps        = 0.4;           % Bed void fraction [-]   
p.eta_meth   = 1.0;           % Effectiveness factor [-]
p.tau        = 0.2;           % Residence time [s]
p.c_CO2_in   = 40.0;          % Inlet CO2 concentration [mol/m^3]
p.St         = 3.4;           % Stanton number [-]
p.dT_ad      = 1074.0;        % Adiabatic temperature rise (ATR) [K]

% Table S2: Kinetic Parameters (Koschany et al., 2016)
p.T_ref      = 550.0;         % Reference temperature [K]
p.k_ref      = 3.46e-4;       % Rate coefficient at T_ref [mol/(bar*s*g_cat)]
p.E_A        = 77500.0;       % Activation energy [J/mol]
p.K_OH_ref   = 0.50;          p.dH_OH  =  22400.0; % [bar^-0.5], [J/mol]
p.K_H2_ref   = 0.44;          p.dH_H2  =  -6200.0; % [bar^-0.5], [J/mol]
p.K_mix_ref  = 0.88;          p.dH_mix = -10000.0; % [bar^-0.5], [J/mol]

% Stoichiometry & Molar Masses: CO2 + 4 H2 <=> CH4 + 2 H2O (Eq. S33)
% Order of species alpha: [CO2, H2, CH4, H2O]
p.nu         = [-1, -4, +1, +2];                 % Stoichiometric coefficients [-]
p.M          = [44.01, 2.016, 16.04, 18.015] * 1e-3; % Molar masses [kg   / mol]

% Stoichiometric feed (H2:CO2 = 4:1 -> x_in = [0.2, 0.8, 0, 0])
x_in         = [0.2, 0.8,  0.0, 0.0];
M_mix_in     = sum(x_in .* p.M);
p.w_in       = (x_in .* p.M) / M_mix_in;         % Inlet mass fractions w_alpha,in [-]

% Single CSTR Residence-Time Constant C_1 = tau / (eps * c_CO2_in) (Eq. S7 / Eq. 5)
p.C_1        = p.tau / (p.eps * p.c_CO2_in);     % [m^3*s/mol] (= 0.0125 for Tab. S1)
slope_CSTR   = 1.0 / p.C_1;                      % Operating line slope [mol/(m^3*s)]

%% 2. Solve Single Isothermal CSTR Across Temperature Range (Eq. S10)
T_vec        = linspace(400, 1000, 250);
X_cstr       = zeros(size(T_vec));
X_eq         = zeros(size(T_vec));

for i = 1:length(T_vec) 
    T = T_vec(i);
    
    % Equilibrium conversion X_eq(T) where R_meth = 0
    X_eq(i)   = fzero(@(X) R_meth_S3(X, T, p), [1e-9, 1.0 - 1e-9]);
    
    % Implicit nonlinear mass balance Eq. (S10): 0 = X_CO2 - C_1 * R_meth(X_CO2)
    X_cstr(i) = fzero(@(X) X - p.C_1 * R_meth_S3(X, T, p), [0.0, X_eq(i)]);
end

%% 3. Plot Results for Single Isothermal CSTR
figure('Color', 'w', 'Position', [100, 100, 1180, 490]);

% SUBPLOT 1: R_CO2 vs. X_CO2 (R-X Plane, Fig. 6 Left)
subplot(1, 2, 1); hold on; box on; grid on;
T_isotherms = [500, 510, 520, 530, 540, 550, 560, 570, 580, 600, 650, 700];

for T_iso =  T_isotherms
    X_eq_iso = fzero(@(X) R_meth_S3(X, T_iso, p), [1e-9, 1.0  - 1e-9]);
    X_grid   = linspace(0, X_eq_iso, 150);
    R_grid   = arrayfun(@(x) R_meth_S3(x, T_iso, p), X_grid);
    
    if T_iso == 550
        plot(X_grid, R_grid, 'Color', [0.25 0.25 0.25], 'LineWidth', 2.0, 'DisplayName', 'Isotherm T = 550 K');
    else
        plot(X_grid, R_grid, 'Color', [0.70 0.70 0.70], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end

% Plot Single CSTR Operating Line: R_meth = (eps * c_CO2_in / tau) * X_CO2
X_line = linspace(0, 1, 100);
plot(X_line, slope_CSTR * X_line, 'Color', [0.13 0.70 0.67], 'LineWidth', 2.2, 'DisplayName', sprintf('Single CSTR Line (slope = %.1f)', slope_CSTR));

% Operating Point at T = 550 K
X_eq_550 = fzero(@(X) R_meth_S3(X, 550, p), [1e-9, 1.0 - 1e-9]);
X_op_550 = fzero(@(X) X - p.C_1 * R_meth_S3(X, 550, p), [0.0, X_eq_550]);
plot(X_op_550, slope_CSTR * X_op_550, 'o', 'Color', [0.13 0.70 0.67], ...
    'MarkerFaceColor', [0.13 0.70 0.67], 'MarkerSize', 8, ...
    'DisplayName', sprintf('OP at 550 K (X = %.3f)', X_op_550));

xlim([0, 1.0]); ylim([0, 110]);
xlabel('Conversion X_{CO_2} / -', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Reaction rate R_{CO_2}/mol m^{-3} s^{-1}', 'FontSize', 11, 'FontWeight', 'bold');
title('Single Isothermal CSTR in R-X Plane', 'FontSize', 12);
legend('Location', 'northwest');

% SUBPLOT 2: X_CO2 vs. T (X-T Plane, Eq. 5 in Figs. 4, 5 & 6 Right)
subplot(1, 2, 2); hold on; box on; grid on;

plot(T_vec, X_eq, 'k:', 'LineWidth', 2.0, 'DisplayName', 'Equilibrium');
plot(T_vec, X_cstr, 'Color', [0.13 0.70 0.67], 'LineWidth', 2.5, 'DisplayName', 'Single CSTR Mass Balance (Eq. 5 / S6)');
plot(550, X_op_550, 'o', 'Color', [0.13 0.70 0.67],  'MarkerFaceColor', [0.13 0.70 0.67], 'MarkerSize', 8,  'DisplayName', 'OP at 550 K');

xlim([400, 1000]); ylim([0, 1.05]);
xlabel('Reaction temp., T / K', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Conversion X_{CO_2} / -', 'FontSize', 11, 'FontWeight', 'bold');
title('Single Isothermal CSTR in X-T Plane (Eq. S6)', 'FontSize', 12);
legend('Location', 'southwest');


%  EFFECTIVE REACTION RATE FUNCTION: Eqs. (S3), (S8), and (S34)-(S39)
function R_meth = R_meth_S3(X_CO2, T, p)
    X_CO2 = max(0.0, min(X_CO2, 1.0 - 1e-12));
    
    % Eq. (S8): Mass fractions w_alpha as a function of X_CO2
    % alpha = 1:CO2, 2:H2, 3:CH4, 4:H2O
    w = -X_CO2 .* (p.nu ./ p.nu(1)) .* (p.M ./ p.M(1)) .* p.w_in(1) + p.w_in;
    
    % Convert mass fractions w_alpha to mole fractions x_alpha & partial pressures p_alpha
    moles_rel = w ./ p.M;
    x_alpha   = moles_rel ./ sum(moles_rel);
    p_alpha   = x_alpha .* p.p_in; % [p_CO2, p_H2, p_CH4, p_H2O] in [bar]
    
    p_CO2 = p_alpha(1);
    p_H2  = p_alpha(2);
    p_CH4 = p_alpha(3);
    p_H2O = p_alpha(4);
    
    % Eqs. (S36)-(S38): Rate, adsorption, and equilibrium constants
    inv_T = (1.0 / p.T_ref) - (1.0 / T);
    k     = p.k_ref     * exp((p.E_A    / p.R_gas) * inv_T);
    K_OH  = p.K_OH_ref  * exp((p.dH_OH  / p.R_gas) * inv_T);
    K_H2  = p.K_H2_ref  * exp((p.dH_H2  / p.R_gas) * inv_T);
    K_mix = p.K_mix_ref * exp((p.dH_mix / p.R_gas) * inv_T);
    K_eq  = 137.0 * (T ^ (-3.998)) * exp(158700.0 / (p.R_gas * T));
    
    % Eqs. (S34)-(S35): Koschany et al. (2016) rate equation [mol/(g_cat*s)]
    DEN    = 1.0 + K_OH * (p_H2O / sqrt(p_H2)) + K_H2 * sqrt(p_H2) + K_mix * sqrt(p_CO2);
    r_meth = k * (p_CO2^0.5) * (p_H2^0.5) * (1.0 - (p_CH4 * p_H2O^2) / (K_eq * p_CO2 * p_H2^4)) / (DEN^2);
    
    % Eq. (S39) & Eq. (S3): Unit conversion & effective volumetric reaction rate [mol/(m^3*s)]
    r_tilde_meth = r_meth * p.rho_cat * 1000.0;
    R_meth       = (1.0 - p.eps) * p.eta_meth * r_tilde_meth;
end
