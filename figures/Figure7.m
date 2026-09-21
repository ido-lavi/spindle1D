% =====================================================================
%  Figure7.m  —  Manuscript Figure 7
%  State diagram of spindle polarity: axial polarity P over the (V||/Vg, M) plane.
%
%  Repository: spindle1D  (reproducibility code for the spindle-polarity paper)
%
%    x-axis : velocity ratio  V||/Vg
%    y-axis : sequential assembly steps M
%    colour : axial polarity P = mean_{spindle} p(x) sgn(x),  P in [-1, 1]
%
%  Input : ../data/sim_M*_eps*.mat  (one steady state per (M, V||/Vg);
%          produced by  PDEsolver.m -> fieldSteadyState.m or filamentSteadyState.m)
%  Output: ../output/Figure7.pdf
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

%% ===== Load all available snapshots =====
files = dir(fullfile(dataDir, 'sim_M*_eps*.mat'));
files = files(~contains({files.name}, 'aself'));   % noncooperative snapshots only (exclude a_self>0 cooperative runs)
if isempty(files)
    error('No sim_M*_eps*.mat found in %s', dataDir);
end

nb = 270;    % filament count prefactor  (n_b in paper)
a0 = 1000;    % a_0: lattice-site production rate (= ell/ell_0)

Ms    = [];
epss  = [];
PLs    = [];
deads = [];

for k = 1:numel(files)
    d  = load(fullfile(dataDir, files(k).name), 'x', 'phi', 'M', 'eps' ,'chain');
    xv = d.x(:);
    rho = d.phi(:,1) + d.phi(:,2);
    p   = (d.phi(:,1) - d.phi(:,2)) ./ (rho + 1e-12);
    mask = rho > 0.2 * max(rho) | abs(xv) < 2;
    B    = mean(p(mask) .* sign(xv(mask)));
    % Dead: no spindle if normalised peak density is below threshold
    Ms(end+1)    = d.M;                               %#ok<AGROW>
    epss(end+1)  = d.eps;                             %#ok<AGROW>
    PLs(end+1)    = B;                                 %#ok<AGROW>
    deads(end+1) = max(rho) / a0 / nb < 0.01;
end

%% ===== Build grid =====
eps_vals = sort(unique(round(epss, 4)));
M_vals   = sort(unique(Ms));
Ne = numel(eps_vals);
Nm = numel(M_vals);

PLgrid     = NaN(Nm, Ne);
dead_grid = false(Nm, Ne);
for k = 1:numel(PLs)
    ie = find(abs(eps_vals - epss(k)) < 1e-6);
    im = find(M_vals == Ms(k));
    PLgrid(im, ie)     = PLs(k);
    dead_grid(im, ie) = deads(k);
end

%% ===== Build figure =====
fs = 10;   % match other figures

% Layout in absolute inches — width fixed at 85 mm (PNAS single col)
fig_w_in = 85/25.4;
ml = 0.36;   % left  margin  (room for ylabel)
mr = 0.09;   % right margin
mb = 0.52;   % bottom margin (xlabel + tick labels)
mt = 0.60;   % top margin    (colorbar + label)

ax_w_in  = fig_w_in - ml - mr;
ax_h_in  = ax_w_in * (7/9);
fig_h_in = mb + ax_h_in + mt;

fig6 = figure('Color', 'w');
set(fig6, 'Units', 'inches', 'Position', [0 0 fig_w_in fig_h_in]);
set(fig6, 'PaperUnits', 'inches', 'PaperSize', [fig_w_in fig_h_in], ...
    'PaperPosition', [0 0 fig_w_in fig_h_in]);

% Normalised axes position derived from absolute margins
ax_l = ml / fig_w_in;
ax_b = mb / fig_h_in;
ax_w = ax_w_in / fig_w_in;
ax_h = ax_h_in / fig_h_in;
ax = axes('Position', [ax_l, ax_b, ax_w, ax_h]);

%% ----- Heatmap -----
set(ax, 'YDir', 'normal');

% Diverging colormap: blue (B=-1) -> white (B=0) -> red (B=+1)
n   = 256;
blu = [0.20 0.45 0.90];
red = [0.85 0.15 0.15];
wht = [1.00 1.00 1.00];
t1  = linspace(0, 1, n/2)';
t2  = linspace(0, 1, n/2)';
cmap = [blu + (wht - blu) .* t1; wht + (red - wht) .* t2];
colormap(ax, cmap);
clim(ax, [-1 1]);

% NaN cells as mid-gray
set(ax, 'Color', [0.96 0.96 0.96]);

%% ----- Filled circles at each data point -----
[eps_grid, M_grid] = meshgrid(eps_vals, M_vals);
valid = ~isnan(PLgrid);
dead  = dead_grid & valid;
alive = valid & ~dead;
hold(ax, 'on');
scatter(ax, eps_grid(alive), M_grid(alive), 100, PLgrid(alive), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
% Dead spindles: dagger marker, no polarity color
text(eps_grid(dead) - 0.001, M_grid(dead) - 0.04, repmat({char(8224)}, numel(eps_grid(dead)), 1), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Color', 'k', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none', 'Parent', ax);
text(eps_grid(dead) - 0.001, M_grid(dead) - 0.04, repmat({char(8224)}, numel(eps_grid(dead)), 1), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Color', blu, 'FontSize', 16, 'Interpreter', 'none', 'Parent', ax);

%% ----- Theoretical border M = 1/eps - 2 -----
eps_th = linspace(max(eps_vals(1), 0.01), eps_vals(end), 300);
M_th   = 1 ./ eps_th - 2;
plot(ax, eps_th, M_th, 'k--', 'LineWidth', 1.2);

hold(ax, 'off');

%% ----- Region labels -----
text(ax, 0.185, 0.53, 'plus-end outward', ...
    'Color', [0.85 0.15 0.15], 'FontSize', fs, 'FontAngle', 'italic', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
text(ax, 0.43, 3.53, 'plus-end inward', ...
    'Color', [0.20 0.45 0.90], 'FontSize', fs, 'FontAngle', 'italic', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

%% ----- Axes styling -----
ylim(ax, [0 6]);
set(ax, 'FontSize', fs, 'Box', 'on', 'LineWidth', 1, 'TickDir', 'out', ...
    'XTick', eps_vals, 'YTick', 0:6, ...
    'XTickLabel', arrayfun(@(v) sprintf('%.1f', v), eps_vals, 'UniformOutput', false));
xlabel(ax, '');
ylabelSplit(ax, 'sequential stages', '$M$', fs, -0.11, 0.68);

%% ----- Colorbar: manual, above axes, 60% width, centered -----
cb_h_in   = 0.10;
cb_gap_in = 0.16;
cb_h  = cb_h_in  / fig_h_in;
cb_gap = cb_gap_in / fig_h_in;
cb_w  = 0.6 * ax_w;
cb_x  = ax_l + (ax_w - cb_w) / 2;
cb_y  = ax_b + ax_h + cb_gap;
cb = colorbar(ax, 'Location', 'manual', 'Position', [cb_x, cb_y, cb_w, cb_h], 'Orientation', 'horizontal');
set(ax, 'Position', [ax_l, ax_b, ax_w, ax_h]);   % restore if MATLAB resized
cb.FontSize   = fs - 2;
cb.Ticks      = [-1 -0.5 0 0.5 1];
cb.TickLabels = {'-1', '-0.5', '0', '0.5', '1'};
cb.Label.String = '';   % suppress built-in label; use annotation instead
cb_lbl_y = cb_y + cb_h + 0.052;
cb_xj    = cb_x + 0.70 * cb_w;   % word/math junction (tune to center the label)
annotation(fig6, 'textbox', [cb_x - 0.2*cb_w, cb_lbl_y, cb_xj - (cb_x - 0.2*cb_w), 0.05], ...
    'String', 'axial polarity ', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'EdgeColor', 'none', 'FontSize', fs, 'Interpreter', 'none', ...
    'FitBoxToText', 'off', 'Margin', 0);
annotation(fig6, 'textbox', [cb_xj + 0.006, cb_lbl_y - 0.004, 0.3*cb_w, 0.05], ...
    'String', '$\mathcal{P}$', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
    'EdgeColor', 'none', 'FontSize', fs, 'Interpreter', 'latex', ...
    'FitBoxToText', 'off', 'Margin', 0);

%% ----- xlabel: text + math side-by-side (one visual line) -----
y_ann  = ax_b * 0.04;
h_ann  = 0.07;
txt_w  = ax_w * 0.70;
gap    = 0.008;
x_off  = 0.025;   % shift xlabel right
annotation(fig6, 'textbox', [ax_l + x_off, y_ann, txt_w, h_ann], ...
    'String', 'sliding/growth velocity ratio ', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontSize', fs, 'Interpreter', 'none', ...
    'FitBoxToText', 'off', 'Margin', 0);
annotation(fig6, 'textbox', [ax_l + x_off + txt_w + gap, y_ann - 0.006, ax_w - txt_w, h_ann], ...
    'String', '$V_{\parallel}/V_{\mathrm{g}}$', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontSize', fs, 'Interpreter', 'latex', ...
    'FitBoxToText', 'off', 'Margin', 0);

%% ===== Export =====
drawnow;
out_file = fullfile(outDir, 'Figure7.pdf');
exportgraphics(fig6, out_file, 'ContentType', 'vector');
fprintf('Saved  %s\n', out_file);

%% ====================================================================
%%  Helper functions
%% ====================================================================

function ylabelSplit(ax, word, mathstr, fs, xlab, y0)
% Split rotated y-label: descriptive WORD in the figure font (below) and
% MATHSTR in a LaTeX box (above), abutting at the junction y0.
% xlab (normalized x, left of the tick labels) and y0 are tunable knobs.
g   = 0.01;  % half-gap between word and math at the junction (raise for more gap)
dxm = 0.01;   % rightward shift (in figure) of the LaTeX box only (raise to move it right)
text(ax, 'Units', 'normalized', 'Position', [xlab, y0-g, 0], 'String', word, ...
    'Rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
    'FontSize', fs, 'Interpreter', 'none');
text(ax, 'Units', 'normalized', 'Position', [xlab+dxm, y0+g, 0], 'String', mathstr, ...
    'Rotation', 90, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', fs, 'Interpreter', 'latex');
end