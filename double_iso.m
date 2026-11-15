%% CO2 methanation in two isothermal CSTRs
clear; clc; close all;

%% Parameters from the supplementary tables
p.R_gas      = 8.314462618;   % Universal gas constant [J/(mol*K)]

% Reactor and catalyst data (Table S1)
p.p_in       = 5.0;           % Inlet pressure [bar]
p.rho_cat    = 2355.2;        % Catalyst density [kg/m^3]
p.eps        = 0.4;           % Bed void fraction [-]
p.eta_meth   = 1.0;           % Catalyst effectiveness factor [-]
p.tau        = 0.2;           % Total cascade residence time [s]

% Inlet CO2 concentration at 550 K
p.c_CO2_in   = 26.365;        % [mol/m^3]
p.C_1        = p.tau / (p.eps * p.c_CO2_in); % Single CSTR constant C_1 [m^3*s/mol]
p.C_2        = p.C_1 / 2.0;                  % 2-Stage constant C_2 = C_1 / 2 [m^3*s/mol]

% Reaction kinetics (Table S2)
p.T_ref      = 550.0;         % Reference temperature [K]
p.k_ref      = 3.46e-4;       % Rate coefficient at T_ref [mol/(bar*s*g_cat)]
p.E_A        = 77500.0;       % Activation energy [J/mol]
p.K_OH_ref   = 0.50;          p.dH_OH  =  22400.0; % [bar^-0.5], [J/mol]
p.K_H2_ref   = 0.44;          p.dH_H2  =  -6200.0; % [bar^-0.5], [J/mol]
p.K_mix_ref  = 0.88;          p.dH_mix = -10000.0; % [bar^-0.5], [J/mol]

% Methanation reaction: CO2 + 4H2 <=> CH4 + 2H2O
% Species order: [CO2, H2, CH4, H2O]
p.nu         = [-1, -4, +1, +2];
p.M          = [44.01, 2.016, 16.04, 18.015] * 1e-3; % [kg/mol]
x_in         = [0.2, 0.8, 0.0, 0.0];                 % Stoichiometric feed (H2/CO2 = 4)
p.w_in       = (x_in .* p.M) / sum(x_in .* p.M);     % Inlet mass fractions [-]

%% Find the conversion and reaction rate at 550 K
T_highlight = 550.0;
[X_tot_550, X_stage_550] = solve_isothermal_cascade(T_highlight, 2, p);
X1_550 = X_stage_550(1);      % Conversion added by stage 1, relative to the original feed
X2_550 = X_stage_550(2);      % Conversion added by stage 2, relative to the original feed
R1_550 = X1_550 / p.C_2;      % Stage 1 rate [mol/(m^3*s)]
R2_550 = X2_550 / p.C_2;      % Stage 2 rate [mol/(m^3*s)]

%% Compare the reactor performance
figure('Color', 'w', 'Position', [80, 100, 1250, 540]);

% Left plot: how the reaction rate changes with conversion
subplot(1, 2, 1); hold on; box on; set(gca, 'FontSize', 11, 'LineWidth', 1.0);

T_isotherms = [500, 510, 520, 530, 540, 550, 560, 570, 600, 700, 750];
for T_iso = T_isotherms
    X_eq_iso = solve_X_eq(T_iso, p);
    X_grid   = linspace(0, X_eq_iso, 250);
    R_grid   = arrayfun(@(x) R_meth_S3(x, T_iso, p), X_grid);

    if T_iso == 550
        plot(X_grid, R_grid, 'Color', [0.45 0.45 0.45], 'LineWidth', 2.5);
    else
        plot(X_grid, R_grid, 'Color', [0.55 0.55 0.55], 'LineWidth', 0.9);
    end
end

slope_1 = 1.0 / p.C_1;
X_line1 = linspace(0, 1.0, 100);
plot(X_line1, slope_1 * X_line1, 'Color', [0.00 0.68 0.84], 'LineWidth', 2.0);
text(0.81, 46, 'single CSTR', 'Color', [0.00 0.68 0.84], 'FontSize', 10, 'Rotation', 23);

% Splitting the residence time between two reactors doubles the line slope.
slope_2 = 1.0 / p.C_2;
plot(X_line1, slope_2 * X_line1, 'Color', [0.96 0.45 0.42], 'LineWidth', 3.0);
text(0.88, 88, 'CSTR^{(1)}', 'Color', [0.85 0.20 0.15], 'FontSize', 10, 'Rotation', 40);

% The second reactor picks up where the first one ends.
X_line2 = linspace(X1_550, 1.0, 100);
plot(X_line2, slope_2 * (X_line2 - X1_550), 'Color', [0.88 0.18 0.08], 'LineWidth', 2.0);
text(0.88, 28, 'CSTR^{(2)}', 'Color', [0.88 0.18 0.08], 'FontSize', 10, 'Rotation', 40);

plot([X1_550, X1_550], [0, R1_550], 'k-', 'LineWidth', 1.0);
plot([X_tot_550, X_tot_550], [0, R2_550], 'k-', 'LineWidth', 1.0);

plot(X1_550, R1_550, 'o', 'Color', [0.55 0.05 0.25], 'MarkerSize', 9, 'LineWidth', 1.5);

plot(X_tot_550, R2_550, 'o', 'Color', [0.45 0.02 0.22], ...
    'MarkerFaceColor', [0.45 0.02 0.22], 'MarkerSize', 10);
text(X_tot_550 - 0.16, R2_550, 'OP 1', 'FontSize', 12, 'FontWeight', 'bold');

text(0.10, 96, 'T = T_{op} = 550 K', 'Color', [0.4 0.4 0.4], 'FontSize', 9, 'Rotation', -28);
text(0.04, 71, '540 K', 'Color', [0.45 0.45 0.45], 'FontSize', 8, 'Rotation', -18);
text(0.04, 50, '530 K', 'Color', [0.45 0.45 0.45], 'FontSize', 8, 'Rotation', -12);
text(0.04, 34, '520 K', 'Color', [0.45 0.45 0.45], 'FontSize', 8, 'Rotation', -8);
text(0.04, 23, '510 K', 'Color', [0.45 0.45 0.45], 'FontSize', 8, 'Rotation', -5);
text(0.04, 15, '500 K', 'Color', [0.45 0.45 0.45], 'FontSize', 8, 'Rotation', -3);

xlim([0, 1.0]); ylim([0, 112]);
xlabel(sprintf('conversion X_{CO_2} / -\n(At 550 K:  X_{CO_2}^{(1)} = %.3f,   X_{CO_2}^{(2)} = %.3f,   Total X_{CO_2}^{(2)} = %.3f)', ...
    X1_550, X2_550, X_tot_550), 'FontSize', 11, 'FontWeight', 'bold');
ylabel('reaction rate R_{CO_2} / mol m^{-3} s^{-1}', 'FontSize', 12, 'FontWeight', 'bold');
title('Figure 6 (Left): R-X Plane for Two Isothermal CSTRs', 'FontSize', 12);

% Right plot: compare conversion for different numbers of CSTRs
subplot(1, 2, 2); hold on; box on; set(gca, 'FontSize', 11, 'LineWidth', 1.0);

T_vec      = linspace(450, 700, 150);
X_eq_vec   = zeros(size(T_vec));
X_n1_vec   = zeros(size(T_vec));
X_n2_tot   = zeros(size(T_vec));
X_n2_stg1  = zeros(size(T_vec));
X_n2_stg2  = zeros(size(T_vec));
X_n10_vec  = zeros(size(T_vec));
X_n100_vec = zeros(size(T_vec));

for i = 1:length(T_vec)
    T = T_vec(i);
    X_eq_vec(i)   = solve_X_eq(T, p);
    X_n1_vec(i)   = solve_isothermal_cascade(T, 1, p);
    [X_n2_tot(i), stg2] = solve_isothermal_cascade(T, 2, p);
    X_n2_stg1(i)  = stg2(1);
    X_n2_stg2(i)  = stg2(2);
    X_n10_vec(i)  = solve_isothermal_cascade(T, 10, p);
    X_n100_vec(i) = solve_isothermal_cascade(T, 100, p);
end

plot(T_vec, X_eq_vec,   'k:', 'LineWidth', 2.0, 'DisplayName', 'equilibrium');
plot(T_vec, X_n100_vec, 'Color', [0.55 0.82 0.88], 'LineWidth', 1.5, 'DisplayName', 'n = 100');
plot(T_vec, X_n10_vec,  'Color', [0.00 0.55 0.65], 'LineWidth', 1.5, 'DisplayName', 'n = 10');
plot(T_vec, X_n2_tot,   'Color', [0.45 0.02 0.22], 'LineWidth', 2.2, 'DisplayName', 'n = 2 (total X_{CO_2}^{(2)})');
plot(T_vec, X_n1_vec,   'Color', [0.96 0.45 0.42], 'LineWidth', 1.5, 'DisplayName', 'n = 1');

plot(T_vec, X_n2_stg1, 'Color', [0.96 0.45 0.42], 'LineWidth', 2.5, 'DisplayName', 'Stage 1: X_{CO_2}^{(1)}');
plot(T_vec, X_n2_stg2, 'Color', [0.88 0.18 0.08], 'LineWidth', 1.8, 'DisplayName', 'Stage 2: X_{CO_2}^{(2)}');

plot([550, 550], [0, X_tot_550], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot(550, X_tot_550, 'o', 'Color', [0.45 0.02 0.22], 'MarkerFaceColor', [0.45 0.02 0.22], ...
    'MarkerSize', 9, 'DisplayName', 'OP 1 (550 K)');
plot(550, X1_550, 'o', 'Color', [0.85 0.20 0.15], 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(550, X2_550, 'o', 'Color', [0.55 0.05 0.25], 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');

text(505, X_tot_550, 'OP 1', 'FontSize', 11, 'FontWeight', 'bold');
text(560, X1_550 - 0.02, 'X_{CO_2}^{(1)}', 'Color', [0.96 0.45 0.42], 'FontSize', 11, 'FontWeight', 'bold');
text(560, X2_550 + 0.01, 'X_{CO_2}^{(2)}', 'Color', [0.55 0.05 0.25], 'FontSize', 11, 'FontWeight', 'bold');

xlim([450, 700]); ylim([0, 1.05]);
xlabel('reaction temp., T / K', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('conversion X_{CO_2} / -', 'FontSize', 12, 'FontWeight', 'bold');
title('Figure 6 (Right): X-T Plane for Two Isothermal CSTRs', 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 9);

%% Functions used for the calculations

function [X_total, X_stage] = solve_isothermal_cascade(T, n_stages, p)
    % Solve the material balance one reactor at a time.
    % Stage conversions are measured relative to the original feed.
    C_n      = p.C_1 / n_stages;
    X_eq     = solve_X_eq(T, p);
    X_prev   = 0.0; % Conversion entering the current reactor
    X_stage  = zeros(1, n_stages);

    for i = 1:n_stages
        % Stop adding conversion once the gas is essentially at equilibrium.
        if (X_eq - X_prev) <= 1e-11
            X_stage(i) = 0.0;
            continue;
        end

        % Material balance for the current reactor (Eq. S26)
        res = @(X_tot) (X_tot - X_prev) - C_n * R_meth_S3(X_tot, T, p);

        % Check both ends before solving with fzero.
        f_low  = res(X_prev); % <= 0 because R_meth(X_prev) >= 0
        f_high = res(X_eq);   % >= 0 because R_meth(X_eq) == 0

        if f_low >= 0.0
            X_tot_sol = X_prev;
        elseif f_high <= 0.0
            X_tot_sol = X_eq;
        else
            X_tot_sol = fzero(res, [X_prev, X_eq]);
        end

        X_stage(i) = X_tot_sol - X_prev;
        X_prev     = X_tot_sol;
    end
    X_total = X_prev;
end

function X_eq = solve_X_eq(T, p)
    X_eq = fzero(@(X) R_meth_S3(X, T, p), [1e-12, 1.0 - 1e-12]);
end

function R_meth = R_meth_S3(X_CO2_cum, T, p)
    X_CO2_cum = max(0.0, min(X_CO2_cum, 1.0 - 1e-14));

    % Update the gas composition using the CO2 conversion (Eq. S8).
    w = -X_CO2_cum .* (p.nu ./ p.nu(1)) .* (p.M ./ p.M(1)) .* p.w_in(1) + p.w_in;

    % Convert mass fractions to mole fractions and partial pressures.
    moles_rel = w ./ p.M;
    x_alpha   = moles_rel ./ sum(moles_rel);
    p_alpha   = x_alpha .* p.p_in;

    p_CO2 = p_alpha(1); p_H2 = p_alpha(2); p_CH4 = p_alpha(3); p_H2O = p_alpha(4);

    % Temperature dependence of the rate and adsorption constants.
    inv_T = (1.0 / p.T_ref) - (1.0 / T);
    k     = p.k_ref     * exp((p.E_A    / p.R_gas) * inv_T);
    K_OH  = p.K_OH_ref  * exp((p.dH_OH  / p.R_gas) * inv_T);
    K_H2  = p.K_H2_ref  * exp((p.dH_H2  / p.R_gas) * inv_T);
    K_mix = p.K_mix_ref * exp((p.dH_mix / p.R_gas) * inv_T);
    K_eq  = 137.0 * (T ^ (-3.998)) * exp(158700.0 / (p.R_gas * T));

    % Methanation rate from the LHHW expression (Eqs. S34-S35).
    DEN    = 1.0 + K_OH * (p_H2O / sqrt(p_H2)) + K_H2 * sqrt(p_H2) + K_mix * sqrt(p_CO2);
    r_meth = k * sqrt(p_CO2 * p_H2) * (1.0 - (p_CH4 * p_H2O^2) / (K_eq * p_CO2 * p_H2^4)) / (DEN^2);

    % Convert the rate to a reactor-bed volume basis.
    % Multiply by 1000 to change catalyst mass from kg to g.
    R_meth = (1.0 - p.eps) * p.eta_meth * (r_meth * p.rho_cat * 1000.0);
end
