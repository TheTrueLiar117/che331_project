%% Single non-isothermal CSTR - CO2 methanation  (version 2, parameters from Table S1/S2)
%  Bremer & Sundmacher (2021), Front. Energy Res. 8:549298
%  Reproduces Fig. 4 (variation of T_op) and Fig. 5 (variation of St).
%
%  Mass balance   (S6/S7, Eq. 5): X = tau/eps * R_meth(T,X) / c_CO2,in      (S-curve)
%  Energy balance (S20,  Eq. 6): X = (1+St)/dTad * (T - Top)                (line)
%  Stability      (Eq. 7)      : dX/dT|mass curve < (1+St)/dTad
%  R_meth = (1-eps)*eta*r_tilde,   r_tilde = r_meth*rho_cat*1000   (S3, S39)
%  Reaction rate: Koschany et al. (2016), Eqs. S34-S38, Table S2.

clear; clc; close all;

%% ---------------- Parameters: reference values of Table S1 ----------------
P.R      = 8.314;      % J/(mol K)
P.p      = 5;          % bar          inlet pressure
P.dTad   = 1074;       % K            adiabatic temperature rise
P.tau    = 0.2;        % s            residence time
P.eps    = 0.4;        % -            bed void fraction
P.eta    = 1;          % -            effectiveness factor
P.rhoCat = 2355.2;     % kg/m3        catalyst density
P.cin    = 40;         % mol/m3       inlet CO2 concentration
St0      = 3.4;        % -            Stanton number (reference)
Top0     = 450;        % K            operating temperature (reference)

% Kinetics, Table S2 (Koschany et al. 2016)
P.Tref   = 550;                          % K
P.kref   = 3.46e-4;                      % mol/(bar g_cat s)
P.EA     = 77.5e3;                       % J/mol
P.KOHref = 0.5;   P.dHOH  =  22.4e3;     % bar^-0.5, J/mol
P.KH2ref = 0.44;  P.dHH2  =  -6.2e3;
P.Kmixref= 0.88;  P.dHmix = -10.0e3;

% Prefactor:  X = Cx * r_meth   (r_meth in mol/(g_cat s))
P.Cx = P.tau/P.eps * (1-P.eps)*P.eta*P.rhoCat*1000 / P.cin;

%% ---------------- Mass-balance curve X(T)  (Eq. 5 / S10) ----------------
Tg  = linspace(400, 1000, 3000);        % reaction temperature range of Table S1
Xm  = zeros(size(Tg));
opts = optimset('TolX', 1e-12);
for i = 1:numel(Tg)
    f     = @(X) X - P.Cx*rate(Tg(i), X, P);
    Xm(i) = fzero(f, [0, 1-1e-6], opts);        % f(0)<0, f(~1)>0 -> unique root
end
dXm = gradient(Xm, Tg);                         % mass-based thermal sensitivity dX/dT
[dXmax, iStar] = max(dXm);
Tstar = Tg(iStar);                              % inflection point T*

%% ---------------- Uniqueness / multiplicity criterion (Eq. 9, 10) ----------------
slope0 = (1+St0)/P.dTad;
fprintf('Steepest mass-curve gradient : dX/dT = %.4f 1/K at T* = %.1f K\n', dXmax, Tstar);
fprintf('Energy-line slope            : (1+St)/dTad = %.4f 1/K  (St = %g)\n', slope0, St0);
if dXmax > slope0
    fprintf('-> multiplicity criterion fulfilled (Eq. 10)\n');
else
    fprintf('-> uniqueness criterion fulfilled (Eq. 9)\n');
end
fprintf('Uniqueness guaranteed for St >= %.1f (at any T_op)\n\n', dXmax*P.dTad - 1);

%% ================= FIGURE 4: variation of operating temperature =================
TopV = linspace(400, 700, 1200);
[parA, XA, TA, opA, nA] = sweep(TopV, 'Top', Tg, Xm, dXm, Tstar, Top0, St0, P);
m = find(nA >= 3);
fprintf('Fig.4: 3 steady states for T_op = %.1f ... %.1f K (St = %g)\n', TopV(m(1)), TopV(m(end)), St0);
TopEx = [TopV(max(m(1)-round(0.25*numel(m)),1)), Top0, TopV(min(m(end)+round(0.25*numel(m)),numel(TopV)))];

figure('Name','Fig 4','Position',[60 80 1150 470]);
subplot(1,2,1); hold on; box on; grid on;
plot(Tg, Xm, 'k', 'LineWidth', 2.2, 'DisplayName', 'mass balance');
cols = [0 0.45 0.74; 0.47 0.67 0.19; 0.85 0.33 0.10];
for j = 1:3
    plot(Tg, slope0*(Tg-TopEx(j)), '-', 'Color', cols(j,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('energy balance, T_{op} = %.0f K', TopEx(j)));
    [Ts, Xs, s] = steadyStates(Tg, Xm, dXm, TopEx(j), St0, P);
    plot(Ts(s),  Xs(s),  'o', 'MarkerFaceColor', cols(j,:), 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
    plot(Ts(~s), Xs(~s), 's', 'MarkerFaceColor', 'w',       'MarkerEdgeColor', cols(j,:), 'LineWidth', 1.5, 'HandleVisibility','off');
    if TopEx(j) == Top0
        op = classify(Ts, s, Tstar);
        for k = 1:numel(Ts)
            text(Ts(k)+8, Xs(k)-0.05, sprintf('OP %d', op(k)), 'FontWeight','bold');
        end
    end
end
xlabel('T / K'); ylabel('X_{CO_2}'); xlim([400 900]); ylim([0 1]);
title(sprintf('Mass and energy balance (St = %g)', St0)); legend('Location','northeast');

subplot(1,2,2); hold on; box on; grid on;
plotBranches(parA, XA, opA);
xlabel('T_{op} / K'); ylabel('X_{CO_2}'); ylim([0 1]); xlim([400 700]);
title('Equilibrium points vs. operating temperature'); legend('Location','southeast');

%% ================= FIGURE 5: variation of Stanton number =================
StV = logspace(log10(0.05), log10(1035), 1500);      % Table S1 range of St
[parB, XB, TB, opB, nB] = sweep(StV, 'St', Tg, Xm, dXm, Tstar, Top0, St0, P);
m = find(nB >= 3);
if isempty(m)
    fprintf('Fig.5: no multiplicity found for T_op = %g K\n', Top0);
    StEx = StV([300 700 1100]);
else
    fprintf('Fig.5: 3 steady states for St = %.2f ... %.2f (T_op = %g K)\n', StV(m(1)), StV(m(end)), Top0);
    StEx = [StV(max(m(1)-round(0.25*numel(m)),1)), St0, StV(min(m(end)+round(0.25*numel(m)),numel(StV)))];
end

figure('Name','Fig 5','Position',[100 60 1150 470]);
subplot(1,2,1); hold on; box on; grid on;
plot(Tg, Xm, 'k', 'LineWidth', 2.2, 'DisplayName', 'mass balance');
for j = 1:3
    sl = (1+StEx(j))/P.dTad;
    plot(Tg, sl*(Tg-Top0), '-', 'Color', cols(j,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('energy balance, St = %.3g', StEx(j)));
    [Ts, Xs, s] = steadyStates(Tg, Xm, dXm, Top0, StEx(j), P);
    plot(Ts(s),  Xs(s),  'o', 'MarkerFaceColor', cols(j,:), 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
    plot(Ts(~s), Xs(~s), 's', 'MarkerFaceColor', 'w',       'MarkerEdgeColor', cols(j,:), 'LineWidth', 1.5, 'HandleVisibility','off');
    if StEx(j) == St0
        op = classify(Ts, s, Tstar);
        for k = 1:numel(Ts)
            text(Ts(k)+8, Xs(k)-0.05, sprintf('OP %d', op(k)), 'FontWeight','bold');
        end
    end
end
xlabel('T / K'); ylabel('X_{CO_2}'); xlim([400 900]); ylim([0 1]);
title(sprintf('Mass and energy balance (T_{op} = %g K)', Top0)); legend('Location','northeast');

subplot(1,2,2); hold on; box on; grid on;
plotBranches(parB, XB, opB);
set(gca, 'XScale', 'log');
xlabel('St'); ylabel('X_{CO_2}'); ylim([0 1]); xlim([StV(1) StV(end)]);
title('Equilibrium points vs. Stanton number'); legend('Location','northeast');

%% ======================= local functions =======================
function r = rate(T, X, P)
% Koschany et al. (2016), Eqs. S34-S38, in mol/(g_cat s)
% basis 1 mol CO2 in (H2:CO2 = 4): CO2 (1-X), H2 4(1-X), CH4 X, H2O 2X, total 5-2X
    R    = P.R;
    nTot = 5 - 2*X;
    pCO2 = P.p*(1-X)./nTot;   pH2  = 4*P.p*(1-X)./nTot;
    pCH4 = P.p*X./nTot;       pH2O = 2*P.p*X./nTot;
    k    = P.kref   *exp( P.EA  /R*(1/P.Tref - 1./T));
    KOH  = P.KOHref *exp( P.dHOH /R*(1/P.Tref - 1./T));
    KH2  = P.KH2ref *exp( P.dHH2 /R*(1/P.Tref - 1./T));
    Kmix = P.Kmixref*exp( P.dHmix/R*(1/P.Tref - 1./T));
    Keq  = 137*T.^(-3.998).*exp(158.7e3./(R*T));          % bar^-2
    drv  = 1 - pCH4.*pH2O.^2./(Keq.*pCO2.*pH2.^4);        % distance from equilibrium
    den  = (1 + KOH.*pH2O./sqrt(pH2) + KH2.*sqrt(pH2) + Kmix.*sqrt(pCO2)).^2;
    r    = k.*sqrt(pH2.*pCO2).*drv./den;
end

function [Tss, Xss, stab] = steadyStates(Tg, Xm, dXm, Top, St, P)
% intersections of mass curve and energy line + stability criterion (Eq. 7)
    slope = (1+St)/P.dTad;
    g     = Xm - slope*(Tg - Top);
    idx   = find(g(1:end-1).*g(2:end) < 0);
    Tss   = Tg(idx) - g(idx).*(Tg(idx+1)-Tg(idx))./(g(idx+1)-g(idx));
    Xss   = slope*(Tss - Top);
    stab  = logical(interp1(Tg, dXm, Tss) < slope);
end

function op = classify(Tss, stab, Tstar)
% OP1 = stable extinguished, OP2 = unstable, OP3 = stable ignited
    op = 2*ones(size(Tss));
    op(stab & Tss <  Tstar) = 1;
    op(stab & Tss >= Tstar) = 3;
end

function [par, Xs, Ts, op, nSS] = sweep(vals, mode, Tg, Xm, dXm, Tstar, Top0, St0, P)
    par = []; Xs = []; Ts = []; op = []; nSS = zeros(size(vals));
    for i = 1:numel(vals)
        if strcmp(mode,'Top'), Top = vals(i); St = St0;
        else,                  Top = Top0;    St = vals(i); end
        [T, X, s] = steadyStates(Tg, Xm, dXm, Top, St, P);
        par = [par, vals(i)*ones(size(X))];  %#ok<AGROW>
        Xs  = [Xs, X];                        %#ok<AGROW>
        Ts  = [Ts, T];                        %#ok<AGROW>
        op  = [op, classify(T, s, Tstar)];    %#ok<AGROW>
        nSS(i) = numel(X);
    end
end

function plotBranches(par, X, op)
    plot(par(op==1), X(op==1), '.', 'Color', [0 0.45 0.74], 'MarkerSize', 9, 'DisplayName', 'OP 1 (stable)');
    plot(par(op==2), X(op==2), '.', 'Color', [0.85 0.33 0.10], 'MarkerSize', 9, 'DisplayName', 'OP 2 (unstable)');
    plot(par(op==3), X(op==3), '.', 'Color', [0.47 0.67 0.19], 'MarkerSize', 9, 'DisplayName', 'OP 3 (stable)');
end
