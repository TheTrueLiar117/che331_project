%% Heat Addition vs. Heat Rejection (van Heerden) Stability Proof for OP 2
%  Reference: Bremer & Sundmacher (2021), Front. Energy Res. 8:549298
%             van Heerden (1953), Ind. Eng. Chem. 45(6), 1242-1247
%
%  Plots:
%    1) Heat Addition Q_gen(T) = dT_ad * X_mass(T) vs.
%       Heat Rejection Q_rem(T) = (1 + St) * (T - T_op)  [van Heerden diagram]
%       with +delta T (ignition -> OP 3) and -delta T (extinction -> OP 1)
%       perturbation trajectories around OP 2.
%    2) Net Heat Accumulation dQ_net(T) = Q_gen(T) - Q_rem(T) (~ dT/dt)
%       showing positive slope at OP 2 (unstable) vs. negative at OP 1, OP 3.
%    3) Mathematical Slope Condition: dX/dT|_mass vs. (1 + St)/dT_ad.

clear; clc; close all;

%% 1. Parameters (Table S1 & Table S2, Koschany et al. 2016)
P.R      = 8.314;      % J/(mol K)
P.p      = 5;          % bar          inlet pressure
P.dTad   = 1074;       % K            adiabatic temperature rise
P.tau    = 0.2;        % s            residence time
P.eps    = 0.4;        % -            bed void fraction
P.eta    = 1;          % -            effectiveness factor
P.rhoCat = 2355.2;     % kg/m3        catalyst density
P.cin    = 40;         % mol/m3       inlet CO2 concentration

St0  = 3.4;            % -            Stanton number (reference)
Top0 = 450;            % K            operating temperature (reference)

% Kinetics, Table S2
P.Tref   = 550;                          % K
P.kref   = 3.46e-4;                      % mol/(bar g_cat s)
P.EA     = 77.5e3;                       % J/mol
P.KOHref = 0.5;   P.dHOH  =  22.4e3;     % bar^-0.5, J/mol
P.KH2ref = 0.44;  P.dHH2  =  -6.2e3;
P.Kmixref= 0.88;  P.dHmix = -10.0e3;

P.Cx = P.tau/P.eps * (1-P.eps)*P.eta*P.rhoCat*1000 / P.cin;

%% 2. Solve Mass-Balance Curve & Compute Heat Addition / Heat Rejection
Tg = linspace(420, 900, 3000);
Xm = zeros(size(Tg));
opts = optimset('TolX', 1e-12);
for i = 1:numel(Tg)
    f     = @(X) X - P.Cx * rate(Tg(i), X, P);
    Xm(i) = fzero(f, [0, 1 - 1e-6], opts);
end

dXm_dT = gradient(Xm, Tg);                 % Mass-based sensitivity dX/dT [1/K]
slope_eb = (1 + St0) / P.dTad;             % Energy-balance critical slope [1/K]

% Temperature-normalized heat rates [K] (Q_dot / (V_dot * rho * c_p))
Q_gen = P.dTad * Xm;                       % Heat addition (chemical generation)
Q_rem = (1 + St0) * (Tg - Top0);           % Heat rejection (convection + cooling)
Q_net = Q_gen - Q_rem;                     % Net accumulation (proportional to dT/dt)
dQnet_dT = gradient(Q_net, Tg);

%% 3. Locate Steady States (OP 1, OP 2, OP 3) & Print Proof Summary
idx = find(Q_net(1:end-1) .* Q_net(2:end) <= 0);
Tss = zeros(1, numel(idx));
Qss = zeros(1, numel(idx));
Xss = zeros(1, numel(idx));
dX_ss = zeros(1, numel(idx));
dQ_ss = zeros(1, numel(idx));

for k = 1:numel(idx)
    i0 = idx(k);
    Tss(k)   = Tg(i0) - Q_net(i0) * (Tg(i0+1) - Tg(i0)) / (Q_net(i0+1) - Q_net(i0));
    Qss(k)   = (1 + St0) * (Tss(k) - Top0);
    Xss(k)   = Qss(k) / P.dTad;
    dX_ss(k) = interp1(Tg, dXm_dT, Tss(k));
    dQ_ss(k) = interp1(Tg, dQnet_dT, Tss(k));
end

fprintf('=========================================================================\n');
fprintf('  OP 2 INSTABILITY VERIFICATION (Top = %g K, St = %g, dTad = %g K)\n', Top0, St0, P.dTad);
fprintf('  Critical Energy-Balance Slope: (1+St)/dTad = %.6f 1/K\n', slope_eb);
fprintf('-------------------------------------------------------------------------\n');
fprintf('%-6s %9s %9s %11s %12s %12s %-10s\n', ...
    'State', 'T_ss [K]', 'X_ss [-]', 'Q_ss [K]', 'dX/dT [1/K]', 'dQnet/dT', 'Stability');
fprintf('-------------------------------------------------------------------------\n');
for k = 1:numel(Tss)
    if dX_ss(k) < slope_eb
        st_str = 'STABLE';
    else
        st_str = 'UNSTABLE';
    end
    fprintf('OP %-3d %9.2f %9.4f %11.2f %12.6f %12.3f %-10s\n', ...
        k, Tss(k), Xss(k), Qss(k), dX_ss(k), dQ_ss(k), st_str);
end
fprintf('=========================================================================\n');

%% 4. Visualization (3-Panel Justification of OP 2 Instability)
figure('Name', 'OP 2 Instability: Heat Addition vs Heat Rejection', ...
       'Color', 'w', 'Position', [60 90 1400 470]);

% --- Panel 1: Heat Addition vs. Heat Rejection (van Heerden Diagram) ---
subplot(1, 3, 1); hold on; box on; grid on;

if numel(Tss) >= 3
    % Shade Q_rem > Q_gen between OP1 and OP2 (cooling / extinction zone)
    m12 = (Tg >= Tss(1)) & (Tg <= Tss(2));
    T12 = Tg(m12);
    fill([T12, fliplr(T12)], [Q_gen(m12), fliplr(Q_rem(m12))], ...
         [0.80 0.90 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.55, ...
         'DisplayName', 'Q_{rem} > Q_{gen} (dT/dt < 0)');

    % Shade Q_gen > Q_rem between OP2 and OP3 (heating / ignition zone)
    m23 = (Tg >= Tss(2)) & (Tg <= Tss(3));
    T23 = Tg(m23);
    fill([T23, fliplr(T23)], [Q_gen(m23), fliplr(Q_rem(m23))], ...
         [1.00 0.85 0.82], 'EdgeColor', 'none', 'FaceAlpha', 0.55, ...
         'DisplayName', 'Q_{gen} > Q_{rem} (dT/dt > 0)');
end

plot(Tg, Q_gen, 'r-',  'LineWidth', 2.4, 'DisplayName', 'Heat Addition: Q_{gen} = \DeltaT_{ad} X_{mass}(T)');
plot(Tg, Q_rem, 'b--', 'LineWidth', 2.0, 'DisplayName', 'Heat Rejection: Q_{rem} = (1+St)(T - T_{op})');

% Mark steady states and perturbation arrows around OP 2
for k = 1:numel(Tss)
    if dX_ss(k) < slope_eb
        plot(Tss(k), Qss(k), 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.15 0.65 0.25], ...
             'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    else
        plot(Tss(k), Qss(k), 's', 'MarkerSize', 10, 'MarkerFaceColor', 'w', ...
             'MarkerEdgeColor', [0.85 0.15 0.10], 'LineWidth', 2.0, 'HandleVisibility', 'off');
    end
    text(Tss(k) + 8, Qss(k) - 45, sprintf('OP %d', k), 'FontWeight', 'bold', 'FontSize', 10);
end

if numel(Tss) >= 3
    % Perturbation arrows at OP 2: -delta T (left) and +delta T (right)
    q_y2 = Qss(2);
    plot([Tss(2)-5, Tss(2)-48], [q_y2, q_y2], 'b-', 'LineWidth', 2.0, 'HandleVisibility', 'off');
    plot(Tss(2)-48, q_y2, '<', 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'HandleVisibility', 'off');
    text(Tss(2)-78, q_y2+55, '-\deltaT \rightarrow OP 1', 'Color', [0 0.25 0.7], 'FontWeight', 'bold', 'FontSize', 9);

    plot([Tss(2)+5, Tss(2)+65], [q_y2, q_y2], 'r-', 'LineWidth', 2.0, 'HandleVisibility', 'off');
    plot(Tss(2)+65, q_y2, '>', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'HandleVisibility', 'off');
    text(Tss(2)+15, q_y2+55, '+\deltaT \rightarrow OP 3', 'Color', [0.75 0.1 0.1], 'FontWeight', 'bold', 'FontSize', 9);
end

xlim([410 880]); ylim([-60 1150]);
xlabel('Reactor Temperature T [K]', 'FontWeight', 'bold');
ylabel('Normalized Heat Rate [K]', 'FontWeight', 'bold');
title('Heat Addition vs. Heat Rejection (van Heerden)', 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 8.2);

% --- Panel 2: Net Heat Accumulation (Q_gen - Q_rem) ---
subplot(1, 3, 2); hold on; box on; grid on;
plot([430 880], [0 0], 'k:', 'LineWidth', 1.3, 'HandleVisibility', 'off');
plot(Tg, Q_net, 'k-', 'LineWidth', 2.2, 'DisplayName', '\DeltaQ_{net} = Q_{gen} - Q_{rem} \propto dT/dt');

for k = 1:numel(Tss)
    if dX_ss(k) < slope_eb
        plot(Tss(k), 0, 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.15 0.65 0.25], ...
             'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    else
        plot(Tss(k), 0, 's', 'MarkerSize', 10, 'MarkerFaceColor', 'w', ...
             'MarkerEdgeColor', [0.85 0.15 0.10], 'LineWidth', 2.0, 'HandleVisibility', 'off');
    end
    text(Tss(k) + 8, -32, sprintf('OP %d', k), 'FontWeight', 'bold', 'FontSize', 10);
end

% Phase-line flow arrows along Delta Q_net = 0
if numel(Tss) >= 3
    T_mid12 = 0.5 * (Tss(1) + Tss(2));
    T_mid23 = 0.5 * (Tss(2) + Tss(3));
    plot(T_mid12, 0, '<', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', ...
         'DisplayName', 'Flow toward OP 1 (dT/dt < 0)');
    plot(T_mid23, 0, '>', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', ...
         'DisplayName', 'Flow toward OP 3 (dT/dt > 0)');
    text(Tss(2)-40, 120, 'd(\DeltaQ_{net})/dT > 0 \Rightarrow Unstable', ...
         'Color', [0.8 0.1 0.1], 'FontWeight', 'bold', 'FontSize', 9);
end

xlim([430 880]); ylim([-350 350]);
xlabel('Reactor Temperature T [K]', 'FontWeight', 'bold');
ylabel('Net Accumulation \DeltaQ_{net} [K]', 'FontWeight', 'bold');
title('Net Thermal Accumulation & Phase Flow', 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 8.5);

% --- Panel 3: Mathematical Slope Criterion ---
subplot(1, 3, 3); hold on; box on; grid on;
plot(Tg, dXm_dT * 1e3, 'k-', 'LineWidth', 2.2, 'DisplayName', 'Mass sensitivity dX/dT \times 10^3');
yline_val = slope_eb * 1e3;
plot([430 880], [yline_val yline_val], 'b--', 'LineWidth', 2.0, ...
     'DisplayName', sprintf('Cooling slope (1+St)/\\DeltaT_{ad} = %.2f \\times 10^{-3}', yline_val));

% Highlight unstable region where dX/dT > slope_eb
unstab = dXm_dT > slope_eb;
if any(unstab)
    Tu = Tg(unstab);
    dXu = dXm_dT(unstab) * 1e3;
    fill([Tu, fliplr(Tu)], [dXu, yline_val*ones(size(Tu))], ...
         [1.0 0.80 0.80], 'EdgeColor', 'none', 'FaceAlpha', 0.65, ...
         'DisplayName', 'Instability zone: dX/dT > (1+St)/\DeltaT_{ad}');
end

for k = 1:numel(Tss)
    if dX_ss(k) < slope_eb
        plot(Tss(k), dX_ss(k)*1e3, 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.15 0.65 0.25], ...
             'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    else
        plot(Tss(k), dX_ss(k)*1e3, 's', 'MarkerSize', 10, 'MarkerFaceColor', 'w', ...
             'MarkerEdgeColor', [0.85 0.15 0.10], 'LineWidth', 2.0, 'HandleVisibility', 'off');
    end
    text(Tss(k) + 8, dX_ss(k)*1e3 + 0.35, sprintf('OP %d', k), 'FontWeight', 'bold', 'FontSize', 10);
end

xlim([430 880]);
xlabel('Reactor Temperature T [K]', 'FontWeight', 'bold');
ylabel('Thermal Sensitivity [10^{-3} K^{-1}]', 'FontWeight', 'bold');
title('Mathematical Slope Condition (Eq. 7)', 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 8.2);

%% Local Rate Function (Koschany et al., 2016)
function r = rate(T, X, P)
    R    = P.R;
    nTot = 5 - 2*X;
    pCO2 = P.p*(1-X)./nTot;   pH2  = 4*P.p*(1-X)./nTot;
    pCH4 = P.p*X./nTot;       pH2O = 2*P.p*X./nTot;
    k    = P.kref   *exp( P.EA  /R*(1/P.Tref - 1./T));
    KOH  = P.KOHref *exp( P.dHOH /R*(1/P.Tref - 1./T));
    KH2  = P.KH2ref *exp( P.dHH2 /R*(1/P.Tref - 1./T));
    Kmix = P.Kmixref*exp( P.dHmix/R*(1/P.Tref - 1./T));
    Keq  = 137*T.^(-3.998).*exp(158.7e3./(R*T));
    drv  = 1 - pCH4.*pH2O.^2./(Keq.*pCO2.*pH2.^4);
    den  = (1 + KOH.*pH2O./sqrt(pH2) + KH2.*sqrt(pH2) + Kmix.*sqrt(pCO2)).^2;
    r    = k.*sqrt(pH2.*pCO2).*drv./den;
end
