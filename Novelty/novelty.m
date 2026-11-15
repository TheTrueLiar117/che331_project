clear; clc; close all;

% quick script to look at unequal reactor volumes (Bremer & Sundmacher idea)
% distribs: uniform, increasing, decreasing, etc.

P.EoR   = 12000;  
P.k0    = 1e10;    
P.dTad  = 250;     
P.tau   = 100;     
P.St    = 1.0;     
P.n     = 4;       

Tcool_vec = linspace(300, 500, 201);

dist = struct('name',{}, 'w',{}, 'phi',{});

dist(1).name = 'Uniform  (1:1:1:1)';         dist(1).w = [1 1 1 1];
dist(2).name = 'Increasing  (1:2:3:4)';       dist(2).w = [1 2 3 4];
dist(3).name = 'Decreasing  (4:3:2:1)';       dist(3).w = [4 3 2 1];
dist(4).name = 'Middle-heavy  (1:4:4:1)';     dist(4).w = [1 4 4 1];
dist(5).name = 'End-heavy  (4:1:1:4)';        dist(5).w = [4 1 1 4];
dist(6).name = 'First-heavy  (7:1:1:1)';      dist(6).w = [7 1 1 1];

nD = numel(dist);
for d = 1:nD
    dist(d).phi = dist(d).w / sum(dist(d).w);
end

res = cell(nD,1);
for d = 1:nD
    phi = dist(d).phi;
    S   = struct('Tcool',[],'Xout',[],'Tmax',[],'Tout',[],...
        'Xall',[],'Tall',[],'stable',[]);
    for Tc = Tcool_vec
        sols = enumerate_steady_states(phi, Tc, P);
        k = size(sols,1);
        S.Tcool  = [S.Tcool;  Tc*ones(k,1)];
        S.Xout   = [S.Xout;   sols(:,1)];
        S.Tmax   = [S.Tmax;   sols(:,2)];
        S.Tout   = [S.Tout;   sols(:,3)];
        S.stable = [S.stable; sols(:,4)];
    end
    res{d} = S;
    fprintf('Done with dist %d: %s\n', d, dist(d).name);
end

fprintf('\n-----------------------------------------------------------------------\n');
fprintf('%-26s %7s %12s %12s %12s\n', ...
    'Distribution', 'Max#SS', 'Bistable[K]', 'PeakT(st)', 'BestX(st)');
fprintf('-----------------------------------------------------------------------\n');

Summary = zeros(nD,5);  
dTc = Tcool_vec(2) - Tcool_vec(1);

for d = 1:nD
    S = res{d};
    nSS = zeros(size(Tcool_vec));
    nSt = nSS;
    for j = 1:numel(Tcool_vec)
        m = abs(S.Tcool - Tcool_vec(j)) < 1e-6;
        nSS(j) = sum(m);
        nSt(j) = sum(m & S.stable);
    end
    bistableWidth = dTc * sum(nSt >= 2);
    multWidth     = dTc * sum(nSS >= 2);
    stab = S.stable == 1;
    if any(stab)
        peakT = max(S.Tmax(stab));
        bestX = max(S.Xout(stab));
    else
        peakT = NaN; bestX = NaN;
    end
    Summary(d,:) = [max(nSS), bistableWidth, multWidth, peakT, bestX];
    fprintf('%-26s %7d %12.1f %12.1f %12.4f\n', ...
        dist(d).name, max(nSS), bistableWidth, peakT, bestX);
end

figure('Name','Bifs','Position',[40 60 1400 900], 'Color','w');

nCols = 3; nRows = ceil(nD/nCols);
for d = 1:nD
    S  = res{d};
    st = S.stable == 1;
    ax = subplot(nRows, nCols, d); hold on; grid on; box on;

    plot(S.Tcool(~st), S.Xout(~st), 'o', ...
        'Color',[0.85 0.2 0.2], 'MarkerSize', 3, ...
        'MarkerFaceColor',[1 0.5 0.5]);
    plot(S.Tcool(st),  S.Xout(st),  '.', ...
        'Color',[0 0.35 0.75], 'MarkerSize', 10);

    xlabel('$T_{\mathrm{cool}} = T_{\mathrm{in}}$ [K]', 'Interpreter','latex','FontSize',11);
    ylabel('$X_n$ (outlet conversion)', 'Interpreter','latex','FontSize',11);
    title(sprintf('%s\n$\\phi$ = [%s]', dist(d).name, ...
        num2str(dist(d).phi,'%.2f  ')), 'Interpreter','latex','FontSize',12);
    ylim([0 1.05]); xlim([Tcool_vec(1) Tcool_vec(end)]);
    set(ax, 'FontSize', 10);
    if d == 1
        legend({'unstable','stable'}, 'Location','best','FontSize',9);
    end
end
sgtitle('Bifurcation Diagrams — Effect of Volume Distribution', 'FontSize',14, 'FontWeight','bold');

figure('Name','Peak Temp','Position',[80 80 900 500], 'Color','w');
hold on; grid on; box on;

markers = {'o','s','d','^','v','>'};
cols = lines(nD);
for d = 1:nD
    S  = res{d};
    st = S.stable == 1;
    plot(S.Tcool(st), S.Tmax(st), markers{d}, ...
        'Color', cols(d,:), 'MarkerSize', 4, ...
        'MarkerFaceColor', cols(d,:), ...
        'DisplayName', dist(d).name);
end
xlabel('$T_{\mathrm{cool}}$ [K]', 'Interpreter','latex','FontSize',12);
ylabel('Peak stage temperature [K]', 'FontSize',12);
title('Peak Temperature of Stable States', 'FontSize',13, 'FontWeight','bold');
legend('Location','northwest','FontSize',9);
xlim([Tcool_vec(1) Tcool_vec(end)]);

figure('Name','Outlet Temp','Position',[120 120 900 500], 'Color','w');
hold on; grid on; box on;
for d = 1:nD
    S  = res{d};
    st = S.stable == 1;
    plot(S.Tcool(st), S.Tout(st), markers{d}, ...
        'Color', cols(d,:), 'MarkerSize', 4, ...
        'MarkerFaceColor', cols(d,:), ...
        'DisplayName', dist(d).name);
end
xlabel('$T_{\mathrm{cool}}$ [K]', 'Interpreter','latex','FontSize',12);
ylabel('$T_{\mathrm{out}}$ [K]', 'FontSize',12);
title('Outlet Temperature of Stable States', 'FontSize',13, 'FontWeight','bold');
legend('Location','northwest','FontSize',9);
xlim([Tcool_vec(1) Tcool_vec(end)]);

figure('Name','Summary Bar','Position',[160 160 1200 400],'Color','w');

shortNames = cell(1,nD);
for d = 1:nD, shortNames{d} = dist(d).name; end

subplot(1,4,1);
bar(Summary(:,1), 'FaceColor', [0.3 0.5 0.8]);
set(gca,'XTickLabel',shortNames,'XTickLabelRotation',35,'FontSize',9);
ylabel('Count'); title('Max # Steady States','FontSize',11);
grid on; box on;

subplot(1,4,2);
bar(Summary(:,2), 'FaceColor', [0.85 0.4 0.2]);
set(gca,'XTickLabel',shortNames,'XTickLabelRotation',35,'FontSize',9);
ylabel('[K]'); title('Bistable Window Width','FontSize',11);
grid on; box on;

subplot(1,4,3);
bar(Summary(:,4), 'FaceColor', [0.8 0.2 0.2]);
set(gca,'XTickLabel',shortNames,'XTickLabelRotation',35,'FontSize',9);
ylabel('[K]'); title('Peak T (stable)','FontSize',11);
grid on; box on;

subplot(1,4,4);
bar(Summary(:,5), 'FaceColor', [0.2 0.65 0.3]);
set(gca,'XTickLabel',shortNames,'XTickLabelRotation',35,'FontSize',9);
ylabel('X'); title('Best Stable Conversion','FontSize',11);
grid on; box on;

sgtitle('Summary Comparison','FontSize',13,'FontWeight','bold');

figure('Name','Num SS vs Tcool','Position',[200 50 900 450], 'Color','w');
hold on; grid on; box on;
for d = 1:nD
    S = res{d};
    nSS_d = zeros(size(Tcool_vec));
    for j = 1:numel(Tcool_vec)
        nSS_d(j) = sum(abs(S.Tcool - Tcool_vec(j)) < 1e-6);
    end
    plot(Tcool_vec, nSS_d, '-', 'Color', cols(d,:), 'LineWidth', 1.5, ...
        'DisplayName', dist(d).name);
end
xlabel('$T_{\mathrm{cool}}$ [K]', 'Interpreter','latex','FontSize',12);
ylabel('Number of steady states', 'FontSize',12);
title('Multiplicity Landscape', 'FontSize',13, 'FontWeight','bold');
legend('Location','best','FontSize',9);
ylim([0 max(Summary(:,1))+1]);
xlim([Tcool_vec(1) Tcool_vec(end)]);

Tc_profile = 400;   % TODO: check other temps if needed

figure('Name','Stage Profiles','Position',[240 80 1300 500], 'Color','w');

for d = 1:nD
    phi = dist(d).phi;
    paths = enum_rec(1, 0, Tc_profile, [], [], phi, Tc_profile, P);

    subplot(2, nCols, d); hold on; grid on; box on;
    nPaths = size(paths,1);
    stageIdx = 0:P.n;  

    for r = 1:nPaths
        Xs = [0, paths(r, 1:P.n)];
        Ts = [Tc_profile, paths(r, P.n+1:2*P.n)];
        isStab = is_stable(paths(r,1:P.n), paths(r,P.n+1:2*P.n), phi, Tc_profile, P);
        if isStab
            plot(stageIdx, Ts, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
                'MarkerFaceColor','b');
        else
            plot(stageIdx, Ts, 'r--s', 'LineWidth', 1.0, 'MarkerSize', 4, ...
                'MarkerFaceColor',[1 0.5 0.5]);
        end
    end
    xlabel('Stage','FontSize',10);
    ylabel('T [K]','FontSize',10);
    title(sprintf('%s', dist(d).name), 'FontSize',10);
    xlim([-0.3 P.n+0.3]);
    if d == 1
        legend({'stable','unstable'}, 'Location','best','FontSize',8);
    end
end
sgtitle(sprintf('Stage Temperature Profiles at $T_{\\mathrm{cool}}$ = %d K', Tc_profile),...
    'Interpreter','latex','FontSize',13);

function sols = enumerate_steady_states(phi, Tc, P)
n = numel(phi);
paths = enum_rec(1, 0, Tc, [], [], phi, Tc, P);
m = size(paths,1);
sols = zeros(m,4);
for r = 1:m
    Xs = paths(r, 1:n);
    Ts = paths(r, n+1:2*n);
    sols(r,1) = Xs(end);        
    sols(r,2) = max(Ts);        
    sols(r,3) = Ts(end);        
    sols(r,4) = is_stable(Xs, Ts, phi, Tc, P);
end

if m > 1
    keep = true(m,1);
    for i = 2:m
        for j = 1:i-1
            if keep(j) && norm(sols(i,1:3) - sols(j,1:3)) < 1e-6
                keep(i) = false;
                break;
            end
        end
    end
    sols = sols(keep,:);
end
end

function out = enum_rec(i, Xin, Tin, Xhist, Thist, phi, Tc, P)
n = numel(phi);
if i > n
    out = [Xhist Thist];
    return
end
Troots = stage_roots(Xin, Tin, Tc, phi(i), P);
out = [];
for T = Troots(:)'
    Da = P.k0 * exp(-P.EoR/T) * phi(i) * P.tau;
    X  = (Xin + Da) / (1 + Da);
    out = [out; enum_rec(i+1, X, T, [Xhist X], [Thist T], phi, Tc, P)]; %#ok<AGROW>
end
end

function Troots = stage_roots(Xin, Tin, Tc, phi_i, P)
tau_i = phi_i * P.tau;
St_i  = P.St * phi_i;

g = @(T) T*(1 + St_i) - Tin - St_i*Tc - ...
    P.dTad * ( (Xin + P.k0*tau_i*exp(-P.EoR./T)) ./ ...
    (1   + P.k0*tau_i*exp(-P.EoR./T)) - Xin );

lo = min(Tin, Tc) - 5;
hi = max(Tin, Tc) + P.dTad + 50;
Tg = linspace(lo, hi, 1000);
gv = g(Tg);

idx = find(gv(1:end-1) .* gv(2:end) < 0);
Troots = zeros(1, numel(idx));
for k = 1:numel(idx)
    Troots(k) = fzero(g, [Tg(idx(k)) Tg(idx(k)+1)]);
end
end

function tf = is_stable(Xs, Ts, phi, Tc, P)
n  = numel(phi);
z0 = reshape([Xs; Ts], [], 1);   
f  = @(z) rhs(z, phi, Tc, P);

J = zeros(2*n);
for j = 1:2*n
    h    = 1e-7 * max(1, abs(z0(j)));
    e    = zeros(2*n,1); e(j) = h;
    J(:,j) = (f(z0+e) - f(z0-e)) / (2*h);
end
tf = max(real(eig(J))) < -1e-10;   
end

function dz = rhs(z, phi, Tc, P)
n  = numel(phi);
dz = zeros(2*n, 1);
Xprev = 0;
Tprev = Tc;       

for i = 1:n
    X     = z(2*i - 1);
    T     = z(2*i);
    tau_i = phi(i) * P.tau;
    r     = P.k0 * exp(-P.EoR / T) * (1 - X);

    dz(2*i - 1) = (Xprev - X)/tau_i + r;
    dz(2*i)     = (Tprev - T)/tau_i + P.dTad*r ...
        - (P.St / P.tau) * (T - Tc);

    Xprev = X;
    Tprev = T;
end
end
