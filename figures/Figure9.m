% =====================================================================
%  Figure9.m  —  Manuscript Figure 9
%  Cooperative model: correct polarity with a single stage (M=1) via
%  self-recruitment. Sweeps the self-recruitment rate a_self (= a_s) at
%  M=1, V||/Vg=0.2.
%
%  Repository: spindle1D  (reproducibility code for the spindle-polarity paper)
%
%  Panels
%    (a,b) steady-state fields + filament snapshot at a_self = 0 and 0.8
%    (c)   density rho(x) family over a_self
%    (d)   polarity p(x) family (x normalised by spindle length Ls)
%    (e)   axial polarity P vs a_self
%
%  Input : ../data/sim_M1_eps0p20_aself*.mat  (cooperative steady states with
%          sampled rods; produced by  PDEsolver.m -> filamentSteadyState.m,
%          run once per a_self in {0,0.2,0.4,0.6,0.7,0.8,0.9}).
%  Output: ../output/Figure9.pdf
%
%  Usage : run from anywhere (paths resolve relative to this file).
% =====================================================================

clear; clc
close all;

% --- repository paths (resolved relative to this file) ---
thisDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(thisDir, '..', 'data');
outDir  = fullfile(thisDir, '..', 'output');
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% ===== Load data =====
M_fixed    = 1;
eps_fixed  = 0.2;
a_self_list = [0 0.2 0.4 0.6 0.7 0.8 0.9];

eps_str = strrep(sprintf('%.2f', eps_fixed), '.', 'p');
N       = numel(a_self_list);
sims    = cell(1, N);
for k = 1:N
    if a_self_list(k) == 0
        % a_self = 0 is the noncooperative M=1 case (shipped as sim_M1_eps0p20.mat)
        fname = fullfile(dataDir, sprintf('sim_M%d_eps%s.mat', M_fixed, eps_str));
    else
        rs_str = strrep(sprintf('%.2f', a_self_list(k)), '.', 'p');
        fname  = fullfile(dataDir, sprintf('sim_M%d_eps%s_aself%s.mat', M_fixed, eps_str, rs_str));
    end
    if ~isfile(fname), error('Missing %s (run PDEsolver.m + filamentSteadyState.m for a_self=%.2f).', fname, a_self_list(k)); end
    sims{k} = load(fname);
end
a_selfvec = a_self_list;
a_self_a = 0.0;   % column (a) example — must match the '$a_s=...$' header string
a_self_b = 0.8;   % column (b) example
col1 = sims{find(abs(a_self_list - a_self_a) < 1e-9, 1)};   % outward / wrong polarity
col2 = sims{find(abs(a_self_list - a_self_b) < 1e-9, 1)};   % inward  / correct polarity
fprintf('Fig9 data: a_self = %s  (M=%d, eps=%.2f)\n', num2str(a_self_list), M_fixed, eps_fixed);

%% ===== Compute B(a_self) =====
cmap = interp1([0; 1], [0.85 0.15 0.15; 0.20 0.45 0.90], linspace(0, 1, N)');
Bvec = zeros(1, N);
for k = 1:N
    sim      = sims{k};
    xv       = sim.x(:);
    rho      = sim.phi(:,1) + sim.phi(:,2);
    p        = (sim.phi(:,1) - sim.phi(:,2)) ./ (rho+1e-12);
    mask     = rho > 0.2 * max(rho) | abs(xv) < 2;
    Bvec(k)  = mean(p(mask) .* sign(xv(mask)));
end

%% ===== Build combined figure =====
fs = 10;

% Dimensions: 17 cm wide x 12 cm tall  (7 cm top + 5 cm bottom)
h_fig2_cm = 7;  h_gap_cm = 0.2;  h_fig3_cm = 5;
h_tot_cm  = h_fig2_cm + h_gap_cm + h_fig3_cm;
h_top_trim_cm = 0.1;   % white space to clip from top of PDF
w_in       = 6.6929;
h_tot_in   = h_tot_cm / 2.54;
h_paper_in = (h_tot_cm - h_top_trim_cm) / 2.54;

fig = figure('Color', 'w');
set(fig, 'Units', 'inches', 'Position', [0 0 w_in h_tot_in]);
set(fig, 'PaperUnits', 'inches', 'PaperSize', [w_in h_paper_in], ...
    'PaperPosition', [0 0 w_in h_tot_in]);

% Normalised height fractions
s3     = h_fig3_cm / h_tot_cm;
s2     = h_fig2_cm / h_tot_cm;
y2_off = (h_fig3_cm + h_gap_cm) / h_tot_cm;

%% ---- Bottom section: panels (c)(d)(e) ----
xL1 = 0.075;  xL2 = 0.41;  xL3 = 0.745;
wp3 = 0.25;
yb3 = 0.28 * s3;
ht3 = 0.66 * s3;

ax1 = axes('Position', [xL1, yb3, wp3, ht3]);
ax2 = axes('Position', [xL2, yb3, wp3, ht3]);
ax3 = axes('Position', [xL3, yb3, wp3, ht3]);

% Panel (c) — density rho(x)
hold(ax1, 'on');
for k = 1:N
    sim   = sims{k};
    rho   = sim.phi(:,1) + sim.phi(:,2);
    rho_n = rho / max(rho);
    mask  = rho_n > 0.2 | abs(sim.x(:)) < 2;
    plot(ax1, sim.x(:), rho_n, 'Color', cmap(k,:), 'LineWidth', 0.5);
    plot(ax1, sim.x(mask), rho_n(mask), 'Color', cmap(k,:), 'LineWidth', 1.5);
end
plot(ax1, sims{1}.x(:), 0.2*ones(1,length(sims{1}.x)), 'LineStyle', '--', 'Color', 'k');
xline(ax1, 0, 'k:', 'LineWidth', 0.8);
hold(ax1, 'off');
xlim(ax1, [-10, 10]);
ylim(ax1, [0, 1.15]);
set(ax1, 'FontSize', fs, 'Box', 'on', 'LineWidth', 1, 'TickDir', 'out', ...
    'YTick', [0 0.5 1], 'YTickLabel', {'0', '0.5', '1'});
ylabelSplit(ax1, 'density', '$\rho(x)/\rho_\mathrm{max}$', fs, -0.26, 0.37);
xlabel(ax1, '$x\ (V_{\mathrm{g}}/\lambda)$', 'FontSize', fs, 'Interpreter', 'latex');

% Panel (d) — polarity p(x), x normalised by spindle length Ls
hold(ax2, 'on');
for k = 1:N
    sim  = sims{k};
    rho  = sim.phi(:,1) + sim.phi(:,2);
    p    = (sim.phi(:,1) - sim.phi(:,2)) ./ (rho+1e-12);
    mask = rho > 0.2*max(rho) | abs(sim.x(:)) < 2;
    Ls   = (sim.x(2)-sim.x(1))*(sum(mask)-1);
    plot(ax2, sim.x(mask)/Ls, p(mask), 'Color', cmap(k,:), 'LineWidth', 1.5);
end
yline(ax2, 0, 'k:', 'LineWidth', 0.8);
xline(ax2, 0, 'k:', 'LineWidth', 0.8);
hold(ax2, 'off');
xlim(ax2, [-0.5, 0.5]);
ylim(ax2, [-1.15, 1.15]);
set(ax2, 'FontSize', fs, 'Box', 'on', 'LineWidth', 1, 'TickDir', 'out', ...
    'YTick', [-1 0 1], 'YTickLabel', {'-1', '0', '1'});
ylabelSplit(ax2, 'polarity', '$p(x)$', fs, -0.24, 0.55);
xlabel(ax2, '$x\ (L_{\mathrm{s}})$', 'FontSize', fs, 'Interpreter', 'latex');

% Panel (e) — axial polarity B(a_self)
hold(ax3, 'on');
plot(ax3, a_selfvec, Bvec, 'k-', 'LineWidth', 1.2);
for k = 1:N
    plot(ax3, a_selfvec(k), Bvec(k), 'o', ...
        'MarkerSize', 5, 'LineWidth', 1.2, ...
        'MarkerFaceColor', cmap(k,:), 'MarkerEdgeColor', cmap(k,:)*0.7);
end
yline(ax3, 0, 'k:', 'LineWidth', 0.8);
hold(ax3, 'off');
xlim(ax3, [a_selfvec(1)-0.05, a_selfvec(end)+0.05]);
ylim(ax3, [-1.15, 1.15]);
set(ax3, 'FontSize', fs, 'Box', 'on', 'LineWidth', 1, 'TickDir', 'out', ...
    'YTick', [-1 0 1], 'YTickLabel', {'-1', '0', '1'}, ...
    'XTick', a_selfvec);
ylabelSplit(ax3, 'axial polarity', '$\mathcal{P}$', fs, -0.24, 0.71);
xlabel(ax3, '$a_{\mathrm{s}}/\lambda$', 'FontSize', fs, 'Interpreter', 'latex');

% Panel labels (c)(d)(e) — placed in the left margin of each panel
lbl_fs       = fs + 2;
lbl_dy       = 0.08 * s3;
lbl_w        = 0.03;
lbl_x_offset = 0.03;
xpos3        = [xL1, xL2, xL3];
for k = 1:3
    annotation(fig, 'textbox', [xpos3(k)-lbl_w-lbl_x_offset, yb3+ht3-lbl_dy, lbl_w, lbl_dy], ...
        'String', sprintf('(%s)', char('c'+k-1)), ...
        'FontSize', lbl_fs, 'FontWeight', 'bold', ...
        'EdgeColor', 'none', 'Interpreter', 'none', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FitBoxToText', 'off', 'Margin', 0);
end

%% ---- Top section: panels (a)(b) ----
xL = 0.09;  xR = 0.56;  wc = 0.41;
y1 = y2_off + 0.54 * s2;
h1 = 0.35   * s2;
y2 = y2_off + 0.15 * s2;
h2 = 0.35   * s2;

yh      = y2_off + 0.89 * s2;
hh      = 0.11   * s2;
dy_math = 0.005;            % LaTeX vertical offset below the text row (raise -> math down)
yh_math = yh - dy_math;
w_lbl   = 0.028;  w_txt = 0.3;  w_eq = 0.16;

% Column (a) header — descriptive text (normal) + a_s value in a latex box
annotation(fig, 'textbox', [xL, yh, w_lbl, hh], ...
    'String', '(a)', 'FontSize', lbl_fs, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FitBoxToText', 'off', 'Margin', 0);
annotation(fig, 'textbox', [xL+w_lbl, yh, w_txt, hh], ...
    'String', 'Plus-end outward polarity,', 'FontSize', fs, 'FontWeight', 'normal', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Interpreter', 'none', 'FitBoxToText', 'off');
annotation(fig, 'textbox', [xL+w_lbl+w_txt-0.053, yh_math, w_eq, hh], ...
    'String', '$a_{\mathrm{s}}/\lambda{=}0$', 'FontSize', fs, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Interpreter', 'latex', 'FitBoxToText', 'off');

% Column (b) header
annotation(fig, 'textbox', [xR, yh, w_lbl, hh], ...
    'String', '(b)', 'FontSize', lbl_fs, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FitBoxToText', 'off', 'Margin', 0);
annotation(fig, 'textbox', [xR+w_lbl, yh, w_txt, hh], ...
    'String', 'Plus-end inward polarity,', 'FontSize', fs, 'FontWeight', 'normal', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Interpreter', 'none', 'FitBoxToText', 'off');
annotation(fig, 'textbox', [xR+w_lbl+w_txt-0.066, yh_math, w_eq, hh], ...
    'String', ['$a_{\mathrm{s}}/\lambda{=}0.8$'], 'FontSize', fs, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Interpreter', 'latex', 'FitBoxToText', 'off');

% Row 1 — steady-state fields
plotFields(axes('Position', [xL, y1, wc, h1]), col1, true);
plotFields(axes('Position', [xR, y1, wc, h1]), col2, false);

% Row 2 — filament snapshots
plotFilamentsInterlace(axes('Position', [xL, y2, wc, h2]), col1, true);
plotFilamentsInterlace(axes('Position', [xR, y2, wc, h2]), col2, false);

%% ===== Export =====
drawnow;
out_file = fullfile(outDir, 'Figure9.pdf');
exportgraphics(fig, out_file, 'ContentType', 'vector');
fprintf('Saved  %s\n', out_file);


%% ====================================================================
%%  Helper functions
%% ====================================================================

function ylabelSplit(ax, word, mathstr, fs, xlab, y0)
% Split rotated y-label: descriptive WORD in the figure font (below) and
% MATHSTR in a LaTeX box (above), abutting at the junction y0.
% xlab (normalized x, left of the tick labels) and y0 are tunable knobs.
g   = 0.015;  % half-gap between word and math at the junction (raise for more gap)
dxm = 0.01;   % rightward shift (in figure) of the LaTeX box only (raise to move it right)
text(ax, 'Units', 'normalized', 'Position', [xlab, y0-g, 0], 'String', word, ...
    'Rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
    'FontSize', fs, 'Interpreter', 'none');
text(ax, 'Units', 'normalized', 'Position', [xlab+dxm, y0+g, 0], 'String', mathstr, ...
    'Rotation', 90, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', fs, 'Interpreter', 'latex');
end

function plotFields(ax, col, show_yl)
fs         = 10;
clr_blue   = [0.15 0.40 0.90];
clr_red    = [0.85 0.15 0.15];
clr_teal   = [0.00 0.62 0.62];
clr_purple = [0.50 0.00 0.75];
clr_teal_d = [0.00 0.35 0.35];
clr_teal_l = [0.55 0.88 0.88];
clr_purp_d = [0.28 0.00 0.42];
clr_purp_l = [0.78 0.55 0.90];

xv = col.x(:);
xr = xv';

s_n  = col.s(:)     / max(col.s);
c_n  = col.c(:);% / max(col.c(:));
phiM = max(col.phi(:));
etaM = max(col.eta(:));
p1   = col.phi(:,1) / phiM;
p2   = col.phi(:,2) / phiM;
e1   = col.eta(:,1) / etaM;
e2   = col.eta(:,2) / etaM;

ymax = 1.28;

hold(ax, 'on');
fill(ax, [xr, fliplr(xr)], [s_n', zeros(1,numel(s_n))], ...
    clr_blue, 'FaceAlpha', 0.40, 'EdgeColor', 'none');
plot(ax, xv, s_n,   'Color', clr_blue,   'LineWidth', 1);
plot(ax, xv, c_n,   'Color', clr_red,    'LineWidth', 1);
plot(ax, xv, e1,    'Color', clr_teal,   'LineStyle', '--', 'LineWidth', 1.5);
plot(ax, xv, e2,    'Color', clr_purple, 'LineStyle', '--', 'LineWidth', 1.5);
% Faint solid lines for each intermediate chain stage phi^(m), m=1..M
if isfield(col, 'chain') && ~isempty(col.chain)
    M_ch = size(col.chain, 1);
    for m = 1:M_ch
        ref  = max(col.chain(m,:,:), [], 'all');  ref = max(ref, 1e-12);
        ch_p = squeeze(col.chain(m,:,1))' / ref;
        ch_m = squeeze(col.chain(m,:,2))' / ref;
        plot(ax, xv, ch_p, 'Color', [clr_teal,   0.22], 'LineStyle', '-', 'LineWidth', 1);
        plot(ax, xv, ch_m, 'Color', [clr_purple, 0.22], 'LineStyle', '-', 'LineWidth', 1);
    end
end
plot(ax, xv, p1,    'Color', clr_teal,   'LineStyle', '-',  'LineWidth', 2);
plot(ax, xv, p2,    'Color', clr_purple, 'LineStyle', '-',  'LineWidth', 2);
hw     = 1;
ya     = 1.14;
dy_arr = ymax * 0.014;
[~, ip] = max(e1);  xpeak_p = xv(ip);
[~, im] = max(e2);  xpeak_m = xv(im);
patch(ax, 'Faces', [1 2 3 4], ...
    'Vertices', [xpeak_p-hw, ya-dy_arr; xpeak_p+0*hw, ya-dy_arr; xpeak_p+0*hw, ya+dy_arr; xpeak_p-hw, ya+dy_arr], ...
    'FaceVertexCData', [clr_teal_d; clr_teal_l; clr_teal_l; clr_teal_d], ...
    'FaceColor', 'interp', 'EdgeColor', 'none');
plot(ax, xpeak_p+0*hw, ya, '->', 'MarkerSize', 7, ...
    'MarkerFaceColor', clr_teal_l, 'MarkerEdgeColor', 'none');
patch(ax, 'Faces', [1 2 3 4], ...
    'Vertices', [xpeak_m-0*hw, ya-dy_arr; xpeak_m+hw, ya-dy_arr; xpeak_m+hw, ya+dy_arr; xpeak_m-0*hw, ya+dy_arr], ...
    'FaceVertexCData', [clr_purp_l; clr_purp_d; clr_purp_d; clr_purp_l], ...
    'FaceColor', 'interp', 'EdgeColor', 'none');
plot(ax, xpeak_m-0*hw, ya, '<-', 'MarkerSize', 7, ...
    'MarkerFaceColor', clr_purp_l, 'MarkerEdgeColor', 'none');
xline(ax, 0, 'k:', 'LineWidth', 0.8);
hold(ax, 'off');

xlim(ax, [-10, 10]);
ylim(ax, [0, ymax]);
set(ax, 'XTickLabel', {}, 'FontSize', fs, 'Box', 'on', ...
    'LineWidth', 1, 'TickDir', 'out', ...
    'YTick', [0 0.5 1.0], 'YTickLabel', {'0', '0.5', '1.0'});
if show_yl
    ylabel(ax, 'norm. amplitude', 'FontSize', fs);
else
    set(ax, 'YTickLabel', {});
end
end

% -----------------------------------------------------------------
function gradRods(ax, rx, ry, row_minus, row_plus, clr_dark, clr_light)
Nf = size(rx, 2);
dy = 0.015;
xm = rx(row_minus, :)';
xp = rx(row_plus,  :)';
yy = ry(1, :)';

all_v = zeros(4*Nf, 2);
all_v(1:4:end, :) = [xm, yy - dy];
all_v(2:4:end, :) = [xp, yy - dy];
all_v(3:4:end, :) = [xp, yy + dy];
all_v(4:4:end, :) = [xm, yy + dy];

all_c = zeros(4*Nf, 3);
all_c(1:4:end, :) = repmat(clr_dark,  Nf, 1);
all_c(2:4:end, :) = repmat(clr_light, Nf, 1);
all_c(3:4:end, :) = repmat(clr_light, Nf, 1);
all_c(4:4:end, :) = repmat(clr_dark,  Nf, 1);

faces = repmat([1 2 3 4], Nf, 1) + repmat((0:Nf-1)' * 4, 1, 4);

patch(ax, 'Faces', faces, 'Vertices', all_v, ...
    'FaceVertexCData', all_c, 'FaceColor', 'interp', 'EdgeColor', 'none');
end

% -----------------------------------------------------------------
function plotFilamentsInterlace(ax, col, show_yl)
fs = 10;

clr_teal_d = [0.00 0.35 0.35];
clr_teal_l = [0.55 0.88 0.88];
clr_purp_d = [0.28 0.00 0.42];
clr_purp_l = [0.78 0.55 0.90];

Nf_p = size(col.rod_x_p, 2);
Nf_m = size(col.rod_x_m, 2);
N    = max(Nf_p, Nf_m);

hold(ax, 'on');
for k = 1:N
    if k <= Nf_p
        gradRods(ax, col.rod_x_p(:,k), col.rod_y_p(:,k), 1, 2, clr_teal_d, clr_teal_l);
        plot(ax, col.rod_x_p(2,k), col.rod_y_p(2,k), '>', ...
            'MarkerSize', 3, 'MarkerFaceColor', clr_teal_l, ...
            'MarkerEdgeColor', 'none', 'LineStyle', 'none');
    end
    if k <= Nf_m
        gradRods(ax, col.rod_x_m(:,k), col.rod_y_m(:,k), 2, 1, clr_purp_d, clr_purp_l);
        plot(ax, col.rod_x_m(1,k), col.rod_y_m(1,k), '<', ...
            'MarkerSize', 3, 'MarkerFaceColor', clr_purp_l, ...
            'MarkerEdgeColor', 'none', 'LineStyle', 'none');
    end
end

s_n   = col.s(:)' / max(col.s);
xv_s  = col.x(:)';
clr_s = [0.15 0.40 0.90];
img_s = repmat(reshape(clr_s, 1, 1, 3), 2, numel(s_n), 1);
ih    = image(ax, 'XData', [xv_s(1), xv_s(end)], 'YData', [-1, 1], 'CData', img_s);
set(ih, 'AlphaData', repmat(s_n * 0.65, 2, 1));

hold(ax, 'off');

set(ax, 'Color', [0.94 0.94 0.94]);
xlim(ax, [-10, 10])
ylim(ax, [-1, 1]);
xlabel(ax, '$x\ (V_{\mathrm{g}}/\lambda)$', 'FontSize', fs, 'Interpreter', 'latex');
set(ax, 'FontSize', fs, 'Box', 'on', 'LineWidth', 1, ...
    'TickDir', 'out', 'YTick', []);
if ~show_yl
    set(ax, 'YTickLabel', {});
end
end