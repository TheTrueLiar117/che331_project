%% Single Unified Validation of the Dispersion / Cascade Criterion (Eqs. 16-19)
% Reference: Bremer & Sundmacher (2021), Front. Energy Res. 8:549298
% Validates Eqs. (17)-(18) by finding the critical threshold Bo_crit = 2*n_crit
% and verifying that hysteresis in the cascade vanishes for Bo > Bo_crit.
clear; clc; close all;

%% 1. Reference Parameters (Supplementary Tables S1 & S2)
p.R_gas      = 8.314462618;
p.p_in       = 5.0;           % [bar]
p.rho_cat    = 2355.2;        % [kg/m^3]
p.eps        = 0.4;           % [-]
p.tau        = 0.2;           % [s]
p.c_CO2_in   = 40.0;          % [mol/m^3]
p.dT_ad      = 1074.0;        % [K]
p.C_1        = p.tau / (p.eps * p.c_CO2_in);

% Kinetic parameters (Table S2)
p.T_ref      = 550.0;         p.k_ref      = 3.46e-4;       p.E_A    = 77500.0;
p.K_OH_ref   = 0.50;          p.dH_OH      =  22400.0;
p.K_H2_ref   = 0.44;          p.dH_H2      =  -6200.0;
p.K_mix_ref  = 0.88;          p.dH_mix     = -10000.0;

p.nu         = [-1, -4, +1, +2];
p.M          = [44.01, 2.016, 16.04, 18.015] * 1e-3;
x_in         = [0.2, 0.8, 0.0, 0.0];
p.w_in       = (x_in .* p.M) / sum(x_in .* p.M);

% Operating setting from Table S1 range (St = 15, eta = 0.05) so transition occurs
% at a clean, fast-to-simulate stage number (n_crit ~ 11, Bo_crit ~ 22)
p.St         = 15.0;
p.eta_meth   = 0.05;

%% 2. Evaluate Unified DM/CM Criterion (Eqs. 17-19) to Find Bo_crit (= 2*n_crit)
n_vec     = 1:25;
Bo_vec    = 2 * n_vec; % Eq. (3): Bo = 2*n
mass_sens = zeros(size(n_vec));
eng_sens  = (1.0 + (2.0 * p.St) ./ Bo_vec) ./ p.dT_ad; % RHS of Eqs. (16)-(18)

T_search  = linspace(460, 880, 250);
for idx = 1:length(n_vec)
    C_Bo = (2.0 / Bo_vec(idx)) * p.C_1; % Eq. (15): C_Bo = (2/Bo)*(tau/(eps*c_in))
    dXdT_vals = arrayfun(@(T) calc_surrogate_dXdT(T, C_Bo, p), T_search);
    mass_sens(idx) = max(dXdT_vals);    % LHS of Eqs. (17)-(18)
end

% Find exact crossover where max(dX_tilde/dT) == (1 + 2*St/Bo)/dT_ad
diff_sens = mass_sens - eng_sens;
idx_cross = find(diff_sens(1:end-1) .* diff_sens(2:end) <= 0, 1, 'last');
n_crit    = interp1(diff_sens(idx_cross:idx_cross+1), n_vec(idx_cross:idx_cross+1), 0);
Bo_crit   = 2.0 * n_crit;

fprintf('=========================================================================\n');
fprintf('  UNIFIED DM / CM CRITERION VALIDATION (Eqs. 17-18)\n');
fprintf('  Critical Threshold: n_crit = %.1f stages (Bo_crit = %.1f)\n', n_crit, Bo_crit);
fprintf('=========================================================================\n');

%% 3. Validate Against Full Cascade Simulations Below, At, and Above n_crit
n_test_cases = [4, 11, 20]; % Below (Multiplicity), At (Critical), Above (Unique)
T_cool_vec   = linspace(420, 620, 130);

figure('Name', 'Validation of Unified Uniqueness Criterion (Eqs. 17-18)', ...
    'Color', 'w', 'Position', [100, 100, 1150, 480]);

% --- LEFT SUBPLOT: Analytical Criterion (Eqs. 17-18) ---
subplot(1, 2, 1); hold on; box on; grid on; set(gca, 'FontSize', 11, 'LineWidth', 1.0);
plot(Bo_vec, mass_sens, '-o', 'Color', [0.00 0.68 0.78], 'LineWidth', 2.2, 'MarkerSize', 4.5, ...
    'DisplayName', 'Mass Sensitivity: max(d\tilde{X}/dT) [Eq. 19]');
plot(Bo_vec, eng_sens,  '-s', 'Color', [0.55 0.02 0.25], 'LineWidth', 2.2, 'MarkerSize', 4.5, ...
    'DisplayName', 'Energy Sensitivity: (1 + 2St/Bo)/\DeltaT_{ad} [Eq. 17]');
xline(Bo_crit, 'k--', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('Critical Threshold Bo_{crit} = %.1f (n_{crit} = %.1f)', Bo_crit, n_crit));

xlim([2, 50]); ylim([0, 0.0045]);
xlabel('Bodenstein Number Bo (= 2n) / -', 'FontSize', 11.5, 'FontWeight', 'bold');
ylabel('Thermal sensitivity / K^{-1}', 'FontSize', 11.5, 'FontWeight', 'bold');
title('A Priori Criterion Prediction (Eqs. 17-18)', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 9.0);

% --- RIGHT SUBPLOT: Full Cascade Verification ---
subplot(1, 2, 2); hold on; box on; grid on; set(gca, 'FontSize', 11, 'LineWidth', 1.0);
colors = {[0.85 0.15 0.15], [0.20 0.20 0.20], [0.00 0.60 0.45]};

for k = 1:length(n_test_cases)
    n_k  = n_test_cases(k);
    Bo_k = 2 * n_k;
    [has_mult, T_sim, X_sim] = simulate_cascade_branches(n_k, T_cool_vec, p);
    
    if has_mult
        status_str = 'Multiplicity / Hysteresis (Eq. 18)';
    else
        status_str = 'Unique Monotonic Curve (Eq. 17)';
    end
    plot(T_sim, X_sim, '.', 'Color', colors{k}, 'MarkerSize', 7, ...
        'DisplayName', sprintf('Bo = %2d (n = %2d): %s', Bo_k, n_k, status_str));
    fprintf('  * Bo = %2d (n = %2d): %s\n', Bo_k, n_k, status_str);
end
fprintf('=========================================================================\n');

xlim([430, 610]); ylim([0, 1.02]);
xlabel('Coolant temperature T_{cool} / K', 'FontSize', 11.5, 'FontWeight', 'bold');
ylabel('Total reactor conversion X_{CO_2}^{(n)} / -', 'FontSize', 11.5, 'FontWeight', 'bold');
title('Full Reactor Simulation Verification', 'FontSize', 12);
legend('Location', 'southeast', 'FontSize', 9.0);

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function sens = calc_surrogate_dXdT(T, C_Bo, p)
    % Solves Eq. (15) for surrogate conversion X_tilde and evaluates Eq. (19)
    X_eq    = fzero(@(X) R_meth_signed(X, T, p), [1e-12, 1.0 - 1e-12]);
    X_tilde = fzero(@(X) X - C_Bo * R_meth_signed(X, T, p), [0.0, X_eq]);
    
    dT = 1e-3; dX = 1e-6;
    dRdT = (R_meth_signed(X_tilde, T + dT, p) - R_meth_signed(X_tilde, T - dT, p)) / (2*dT);
    X_hi = min(X_tilde + dX, 1-1e-12); X_lo = max(X_tilde - dX, 0.0);
    dRdX = (R_meth_signed(X_hi, T, p) - R_meth_signed(X_lo, T, p)) / (X_hi - X_lo);
    sens = (C_Bo * dRdT) / (1.0 - C_Bo * dRdX); % Eq. (19)
end

function [has_mult, T_all, X_all] = simulate_cascade_branches(n_stages, T_cool_vec, p)
    paths = ones(3, n_stages);
    paths(2, 1) = 2; % [2-1-...]
    paths(3, 1) = 3; % [3-1-...]
    T_all = []; X_all = []; mult_pts = 0;
    
    for idx_T = 1:length(T_cool_vec)
        Tc = T_cool_vec(idx_T);
        valid_here = 0;
        for b = 1:3
            [valid, X_tot] = eval_single_path(Tc, paths(b, :), p);
            if valid
                T_all(end+1) = Tc; X_all(end+1) = X_tot; %#ok<AGROW>
                valid_here = valid_here + 1;
            end
        end
        if valid_here > 1, mult_pts = mult_pts + 1; end
    end
    has_mult = (mult_pts > 0);
end

function [valid, X_total] = eval_single_path(T_cool, state_seq, p)
    n_stages = length(state_seq);
    St_n = p.St / n_stages; C_n = p.C_1 / n_stages;
    X_cum = 0.0; T_prev = T_cool; valid = true;
    
    for i = 1:n_stages
        req = state_seq(i);
        dT_ad_i = (1.0 - X_cum) * p.dT_ad;
        T_op_i  = (T_prev + St_n * T_cool) / (1.0 + St_n);
        if dT_ad_i <= 1e-7
            if req > 1, valid = false; X_total = nan; return; end
            continue;
        end
        T_max = T_op_i + dT_ad_i / (1.0 + St_n);
        res_f = @(T) ((1+St_n)/dT_ad_i)*(T - T_op_i) - ...
            C_n * R_meth_signed(X_cum + ((1+St_n)/dT_ad_i)*(T - T_op_i), T, p);
        
        T_scan = linspace(T_op_i + 1e-7, T_max - 1e-6, 55);
        f_scan = arrayfun(res_f, T_scan);
        r_cnt  = 0; T_found = nan; X_found = nan;
        for k = 1:(length(T_scan)-1)
            if f_scan(k) * f_scan(k+1) <= 0
                r_cnt = r_cnt + 1;
                if r_cnt == req
                    T_found = fzero(res_f, [T_scan(k), T_scan(k+1)]);
                    X_found = max(0.0, ((1+St_n)/dT_ad_i)*(T_found - T_op_i));
                    break;
                end
            end
        end
        if isnan(T_found), valid = false; X_total = nan; return; end
        T_prev = T_found; X_cum = X_cum + X_found;
    end
    X_total = X_cum;
end

function R_meth = R_meth_signed(X_CO2_cum, T, p)
    X_CO2_cum = max(0.0, min(X_CO2_cum, 1.0 - 1e-14));
    w = -X_CO2_cum .* (p.nu ./ p.nu(1)) .* (p.M ./ p.M(1)) .* p.w_in(1) + p.w_in;
    moles_rel = w ./ p.M; x_alpha = moles_rel ./ sum(moles_rel);
    p_alpha = x_alpha .* p.p_in;
    p_CO2 = p_alpha(1); p_H2 = p_alpha(2); p_CH4 = p_alpha(3); p_H2O = p_alpha(4);
    
    inv_T = (1.0 / p.T_ref) - (1.0 / T);
    k     = p.k_ref     * exp((p.E_A    / p.R_gas) * inv_T);
    K_OH  = p.K_OH_ref  * exp((p.dH_OH  / p.R_gas) * inv_T);
    K_H2  = p.K_H2_ref  * exp((p.dH_H2  / p.R_gas) * inv_T);
    K_mix = p.K_mix_ref * exp((p.dH_mix / p.R_gas) * inv_T);
    K_eq  = 137.0 * (T ^ (-3.998)) * exp(158700.0 / (p.R_gas * T));
    
    DEN    = 1.0 + K_OH * (p_H2O / sqrt(p_H2)) + K_H2 * sqrt(p_H2) + K_mix * sqrt(p_CO2);
    r_meth = k * sqrt(p_CO2 * p_H2) * (1.0 - (p_CH4 * p_H2O^2) / (K_eq * p_CO2 * p_H2^4)) / (DEN^2);
    R_meth = (1.0 - p.eps) * p.eta_meth * (r_meth * p.rho_cat * 1000.0);
end