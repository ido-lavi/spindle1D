% =====================================================================
%  filamentSteadyState.m  —  reconstruct steady state psi, sample rods, and save
%
%  Repository: spindle1D  (reproducibility code for the spindle-polarity paper)
%
%  Run in the SAME workspace right after PDEsolver.m reaches steady state.
%  Reconstructs the full steady-state number density psi_pm(x,l) from the
%  nucleation field, samples individual filaments (rods), shows diagnostic
%  figures (also written to ../output), and saves a snapshot WITH rods to
%  ../data/sim_M<M>_eps<eps>.mat ('_aself<a_self>' appended when a_self > 0).
%
%  Saved: x, L, eps, M, phi, eta, s, c, pol  (+ chain if M > 0)
%         + rod_x_p, rod_y_p, rod_x_m, rod_y_m.
%
%  Steady-state formula (dimensionless, Vg = lambda = 1):
%    psi_+(x,l) = n_+(x + eps*l) exp(-l)   [col 1, right-pointing]
%    psi_-(x,l) = n_-(x - eps*l) exp(-l)   [col 2, left-pointing]
%
%  This is the field+filament saver (Figures 3, 6, 9). For the fields-only
%  snapshot (Figure 7 grid, no rods) use fieldSteadyState.m instead.
%
%  Requires the workspace left by PDEsolver.m: nuc, hat_nuc, xi, eta, phi, x, k,
%  del, eps, N, L, denom_pm, a0.
%
%  Usage:  PDEsolver.m  ->  filamentSteadyState.m
% =====================================================================
close all;

% --- repository paths (resolved relative to this file) ---
thisDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(thisDir, '..', 'data');
outDir  = fullfile(thisDir, '..', 'output');
if ~exist(dataDir, 'dir'); mkdir(dataDir); end
if ~exist(outDir,  'dir'); mkdir(outDir);  end

%% Length grid
lmax = 6;
Nl   = 400;
l    = linspace(0, lmax, Nl);
dl   = l(2) - l(1);

%% Spectral reconstruction of psi(x, l): N × Nl × 2
Psi = zeros(N, Nl, 2);
for i = 1:Nl
    decay = exp(-l(i));
    Psi(:,i,1) = real(ifft(hat_nuc(:,1) .* exp( 1i*k'*eps*l(i)))) * decay;
    Psi(:,i,2) = real(ifft(hat_nuc(:,2) .* exp(-1i*k'*eps*l(i)))) * decay;
end

%% Check consistency: recover xi (minus ends), eta (plus ends), phi (mass) from Psi
% xi_+(x) = int_0^inf psi_+(x,l) dl   [integrate over l, no shift]
% xi_-(x) = int_0^inf psi_-(x,l) dl   [integrate over l, no shift]
xi_rec = squeeze(sum(Psi, 2)) * dl;     % N×2 x-polarity

% eta_+(x) = int_0^inf psi_+(x-l, l) dl   [plus end at x => shift by -l]
% eta_-(x) = int_0^inf psi_-(x+l, l) dl   [plus end at x => shift by +l]
hat_Psi1 = fft(Psi(:,:,1));     % N×Nl k-l space
hat_Psi2 = fft(Psi(:,:,2));     % N×Nl k-l space
phase_neg = exp(-1i * k' * l);  % N×Nl  shift by -l (col 1)
phase_pos = exp( 1i * k' * l);  % N×Nl  shift by +l (col 2)
hat_eta_rec = [sum(hat_Psi1 .* phase_neg, 2), ...
               sum(hat_Psi2 .* phase_pos, 2)] * dl;    % N×2 k-polarity
eta_rec = real(ifft(hat_eta_rec));  % N×2 x-polarity

% phi^(0)_+(x) = a0 * int_0^inf dl int_0^l ds psi_+(x-s, l)   [Eq. 21]
% phi^(0)_-(x) = a0 * int_0^inf dl int_0^l ds psi_-(x+s, l)
% Inner integral over s done analytically in Fourier:
%   col 1: int_0^l e^{-iks} ds = (1 - e^{-ikl}) / (ik),  limit l at k=0
%   col 2: int_0^l e^{+iks} ds = (e^{+ikl} - 1) / (ik),  limit l at k=0
kk  = k';           % N×1
ll  = l;            % 1×Nl
nz  = kk ~= 0;      % non-zero wavenumbers

fac1 = zeros(N, Nl);
fac2 = zeros(N, Nl);
fac1(nz,:)  = (1 - exp(-1i*kk(nz)*ll)) ./ (1i*kk(nz));
fac1(~nz,:) = repmat(ll, sum(~nz), 1);
fac2(nz,:)  = (exp(1i*kk(nz)*ll) - 1) ./ (1i*kk(nz));
fac2(~nz,:) = repmat(ll, sum(~nz), 1);

hat_phi_rec = [sum(hat_Psi1 .* fac1, 2), sum(hat_Psi2 .* fac2, 2)] * (dl * a0);
phi_rec = real(ifft(hat_phi_rec));

% Comparison plots
fs = 16;
figure;

subplot(3,1,1)
plot(x, xi(:,1),     'k',   'LineWidth', 2); hold on
plot(x, xi_rec(:,1), '--r', 'LineWidth', 1.5);
plot(x, xi(:,2),     'k--', 'LineWidth', 2);
plot(x, xi_rec(:,2), ':r',  'LineWidth', 1.5);
lgd = legend('\xi_+ sim','\xi_+ rec','\xi_- sim','\xi_- rec');
lgd.FontSize = fs;
title('Minus-end density \xi_\pm'); xlabel('x'); ylabel('AU')
set(gca, 'FontSize', fs)

subplot(3,1,2)
plot(x, eta(:,1),     'k',   'LineWidth', 2); hold on
plot(x, eta_rec(:,1), '--r', 'LineWidth', 1.5);
plot(x, eta(:,2),     'k--', 'LineWidth', 2);
plot(x, eta_rec(:,2), ':r',  'LineWidth', 1.5);
lgd = legend('\eta_+ sim','\eta_+ rec','\eta_- sim','\eta_- rec');
lgd.FontSize = fs;
title('Plus-end density \eta_\pm'); xlabel('x'); ylabel('AU')
set(gca, 'FontSize', fs)

subplot(3,1,3)
plot(x, phi(:,1),     'k',   'LineWidth', 2); hold on
plot(x, phi_rec(:,1), '--r', 'LineWidth', 1.5);
plot(x, phi(:,2),     'k--', 'LineWidth', 2);
plot(x, phi_rec(:,2), ':r',  'LineWidth', 1.5);
lgd = legend('\phi_+ sim','\phi_+ rec','\phi_- sim','\phi_- rec');
lgd.FontSize = fs;
title('MT lattice \phi_\pm'); xlabel('x'); ylabel('AU')
set(gca, 'FontSize', fs)

%% Sample and plot filaments as rods
%
% Sampling from the 2-D density Psi(x, l):
%   Psi(i,j,pm) is the number density of filaments per unit x per unit l at
%   grid point (x_i, l_j).  Multiplied by the cell area dx*dl it gives the
%   expected filament count in that cell.  We treat the N×Nl grid values as
%   an unnormalized probability mass function over discrete (x,l) pairs,
%   clamp any numerical negatives to zero, and normalise to sum 1.
%   randsample then draws Nf independent indices from this distribution
%   (with replacement), and ind2sub converts each flat index back to a grid
%   position (ix, il).  The result is a set of (x_minus, l) pairs whose
%   joint histogram converges to Psi as Nf -> inf.

dx      = x(2) - x(1);
Nf_p    = ceil(sum(nuc(:,1)) * dx);   % steady-state right-pointing count (R_plus / lambda, lambda=1)
Nf_m    = ceil(sum(nuc(:,2)) * dx);   % steady-state left-pointing count
y_range = [-1, 1];      % vertical range for random rod placement

% max(...,0): element-wise clamp of negatives (FFT noise); does NOT find a max.
% w(:): flatten N×Nl matrix to N*Nl vector; divide by sum to get a PMF.
w_p = max(Psi(:,:,1), 0);  w_p = w_p(:) / sum(w_p(:));
w_m = max(Psi(:,:,2), 0);  w_m = w_m(:) / sum(w_m(:));

% Draw Nf flat indices from 1:N*Nl with replacement; cell i is chosen with
% probability w(i), so dense regions of Psi are sampled more often.
idx_p = randsample(N*Nl, Nf_p, true, w_p);
idx_m = randsample(N*Nl, Nf_m, true, w_m);

% Convert each flat index back to (x-grid, l-grid) subscripts.
[ix_p, il_p] = ind2sub([N, Nl], idx_p);
[ix_m, il_m] = ind2sub([N, Nl], idx_m);

% right-pointing: minus end at x(ix), plus end at x(ix) + l(il)
% left-pointing:  minus end at x(ix), plus end at x(ix) - l(il)
xm_p = x(ix_p)';   lp = l(il_p)';
xm_m = x(ix_m)';   lm = l(il_m)';

yp = y_range(1) + diff(y_range) * rand(Nf_p, 1);
ym = y_range(1) + diff(y_range) * rand(Nf_m, 1);

rod_x_p = [xm_p,      xm_p + lp, NaN(Nf_p,1)]';   % 3 × Nf_p
rod_y_p = [yp,        yp,        NaN(Nf_p,1)]';
rod_x_m = [xm_m - lm, xm_m,     NaN(Nf_m,1)]';   % 3 × Nf_m
rod_y_m = [ym,        ym,        NaN(Nf_m,1)]';

%% Static 4-panel figure matching filament_anim layout
fig_snap = figure;
fig_snap.Position = [50 30 1400 1000];
set(fig_snap, 'color', 'k');

% Panel 1 (top): nucleation rate
ax_nuc = axes('Position', [0.05, 0.71, 0.90, 0.18]);
set(ax_nuc, 'color', 'k', 'XColor', 'w', 'YColor', 'w', ...
    'YAxisLocation', 'origin', 'FontSize', 16);
hold(ax_nuc, 'on');
plot(ax_nuc, x, nuc(:,1), 'Color', 'c', 'LineWidth', 2.5);
plot(ax_nuc, x, nuc(:,2), 'Color', 'm', 'LineWidth', 2.5);
xlim(ax_nuc, [-L/2, L/2]);
ylim(ax_nuc, [0, 1.1*max(nuc(:))]);
title(ax_nuc, sprintf('Steady state nucleation rate  (M=%d)', M), 'Color', 'w', 'FontSize', 20);
xlabel(ax_nuc, 'x  (V_g/\lambda)', 'Color', 'w', 'FontSize', 16);
pbaspect(ax_nuc, [1, diff(y_range)/L, 1]);

% Panel 2: minus ends only
ax_xm = axes('Position', [0.05, 0.42, 0.90, 0.18]);
set(ax_xm, 'color', 'k', 'XColor', 'w', 'YColor', 'none', 'FontSize', 16);
hold(ax_xm, 'on');
plot(ax_xm, xm_p, yp, 'o', 'Color', 'c', 'MarkerSize', 4, 'MarkerFaceColor', 'c', 'LineStyle', 'none');
plot(ax_xm, xm_m, ym, 'o', 'Color', 'm', 'MarkerSize', 4, 'MarkerFaceColor', 'm', 'LineStyle', 'none');
xlim(ax_xm, [-L/2, L/2]); ylim(ax_xm, y_range + [-0.05 0.05]);
pbaspect(ax_xm, [1, diff(y_range)/L, 1]);
title(ax_xm, sprintf('Steady-state sample:  N = %d', Nf_p + Nf_m), 'Color', 'w', 'FontSize', 20);
set(ax_xm, 'XTickLabel', []);

% Panel 3: plus ends only
ax_xp = axes('Position', [0.05, 0.24, 0.90, 0.18]);
set(ax_xp, 'color', 'k', 'XColor', 'w', 'YColor', 'none', 'FontSize', 16);
hold(ax_xp, 'on');
plot(ax_xp, xm_p + lp, yp, '>', 'Color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'c', 'LineStyle', 'none');
plot(ax_xp, xm_m - lm, ym, '<', 'Color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'm', 'LineStyle', 'none');
xlim(ax_xp, [-L/2, L/2]); ylim(ax_xp, y_range + [-0.05 0.05]);
pbaspect(ax_xp, [1, diff(y_range)/L, 1]);
set(ax_xp, 'XTickLabel', []);

% Panel 4 (bottom): full filaments + both ends
ax_all = axes('Position', [0.05, 0.06, 0.90, 0.18]);
set(ax_all, 'color', 'k', 'XColor', 'w', 'YColor', 'none', 'FontSize', 16);
hold(ax_all, 'on');
plot(ax_all, rod_x_p(:), rod_y_p(:), 'Color', 'c', 'LineWidth', 1.8);
plot(ax_all, xm_p + lp, yp, '>', 'Color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'c', 'LineStyle', 'none');
plot(ax_all, xm_p,      yp, 'o', 'Color', 'c', 'MarkerSize', 4, 'MarkerFaceColor', 'c', 'LineStyle', 'none');
plot(ax_all, rod_x_m(:), rod_y_m(:), 'Color', 'm', 'LineWidth', 1.8);
plot(ax_all, xm_m - lm, ym, '<', 'Color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'm', 'LineStyle', 'none');
plot(ax_all, xm_m,      ym, 'o', 'Color', 'm', 'MarkerSize', 4, 'MarkerFaceColor', 'm', 'LineStyle', 'none');
xlim(ax_all, [-L/2, L/2]); ylim(ax_all, y_range + [-0.05 0.05]);
pbaspect(ax_all, [1, diff(y_range)/L, 1]);
set(ax_all, 'XTickLabel', [], 'YTick', [], 'YColor', 'none');

text(ax_xm,  -9.8, 0.88, 'minus-ends', 'Color', 'w', 'FontSize', 20, 'VerticalAlignment', 'top');
text(ax_xp,  -9.8, 0.88, 'plus-ends',  'Color', 'w', 'FontSize', 20, 'VerticalAlignment', 'top');
text(ax_all, -9.8, 0.88, 'filaments',  'Color', 'w', 'FontSize', 20, 'VerticalAlignment', 'top');
linkaxes([ax_nuc, ax_xm, ax_xp, ax_all], 'x');
drawnow;
exportgraphics(fig_snap, fullfile(outDir, sprintf('psi_M%d_eps%.3f.png', M, eps)), 'BackgroundColor', 'k', 'Resolution', 150);

%% Compute histograms and coarse exact psi for both polarities
Nf_hist = 30000;
Nbx     = 80;
Nbl     = 60;

idx_hist_p = randsample(N*Nl, Nf_hist, true, w_p);
[ix_h_p, il_h_p] = ind2sub([N, Nl], idx_hist_p);
x_hist_p = x(ix_h_p)';   l_hist_p = l(il_h_p)';

idx_hist_m = randsample(N*Nl, Nf_hist, true, w_m);
[ix_h_m, il_h_m] = ind2sub([N, Nl], idx_hist_m);
x_hist_m = x(ix_h_m)';   l_hist_m = l(il_h_m)';

x_edges = linspace(-L/2, L/2, Nbx+1);
l_edges = linspace(0, lmax, Nbl+1);
x_ctrs  = 0.5*(x_edges(1:end-1) + x_edges(2:end));
l_ctrs  = 0.5*(l_edges(1:end-1) + l_edges(2:end));
dx_bin  = x_edges(2) - x_edges(1);
dl_bin  = l_edges(2) - l_edges(1);

H_p = histcounts2(x_hist_p, l_hist_p, x_edges, l_edges);
H_m = histcounts2(x_hist_m, l_hist_m, x_edges, l_edges);
H_density_p = H_p * (sum(Psi(:,:,1),'all')*dx*dl) / (Nf_hist * dx_bin * dl_bin);
H_density_m = H_m * (sum(Psi(:,:,2),'all')*dx*dl) / (Nf_hist * dx_bin * dl_bin);

bx = min(Nbx, max(1, ceil((x(:) + L/2) / dx_bin)));
bl = min(Nbl, max(1, ceil( l(:)        / dl_bin)));
[II, JJ] = ndgrid(1:N, 1:Nl);
bx_all = bx(II(:));   bl_all = bl(JJ(:));
vals_p = Psi(:,:,1);
vals_m = Psi(:,:,2);
Psi_coarse_p = accumarray([bx_all, bl_all], vals_p(:), [Nbx, Nbl], @mean, 0);
Psi_coarse_m = accumarray([bx_all, bl_all], vals_m(:), [Nbx, Nbl], @mean, 0);
cmax_p = max(Psi_coarse_p(:));
cmax_m = max(Psi_coarse_m(:));

n_cm = 256;
cm_cyan    = [zeros(n_cm,1), linspace(0,1,n_cm)', linspace(0,1,n_cm)'];
cm_magenta = [linspace(0,1,n_cm)', zeros(n_cm,1), linspace(0,1,n_cm)'];

%% Figure 3: nuc (top) | exact psi_-, psi_+ (middle) | histograms (bottom)
fig_psi = figure;
fig_psi.Position = [50 50 1400 1200];
set(fig_psi, 'color', 'k');

% Top: nucleation rate (full width)
ax_nuc3 = axes('Position', [0.05, 0.71, 0.90, 0.18]);
set(ax_nuc3, 'color', 'k', 'XColor', 'w', 'YColor', 'w', ...
    'YAxisLocation', 'origin', 'FontSize', 16);
hold(ax_nuc3, 'on');
plot(ax_nuc3, x, nuc(:,1), 'Color', 'c', 'LineWidth', 2.5);
plot(ax_nuc3, x, nuc(:,2), 'Color', 'm', 'LineWidth', 2.5);
xlim(ax_nuc3, [-L/2, L/2]);  ylim(ax_nuc3, [0, 1.1*max(nuc(:))]);
title(ax_nuc3, sprintf('Steady state nucleation rate  (M=%d)', M), 'Color', 'w', 'FontSize', 20);
xlabel(ax_nuc3, 'x  (V_g/\lambda)', 'Color', 'w', 'FontSize', 16);
pbaspect(ax_nuc3, [1, diff(y_range)/L, 1]);

% Middle-left: exact psi_- (magenta colormap)
ax_em = axes('Position', [0.07, 0.35, 0.41, 0.30]);
imagesc(ax_em, x_ctrs, l_ctrs, Psi_coarse_m');
axis(ax_em, 'xy');
colormap(ax_em, cm_magenta);
cb = colorbar(ax_em);  cb.Color = 'w';
xlabel(ax_em, 'minus-end position x (V_g/\lambda)', 'Color', 'w', 'FontSize', 16);  ylabel(ax_em, 'length l  (V_g/\lambda)', 'Color', 'w', 'FontSize', 16);
title(ax_em, '\psi_-(x,l)  [exact]', 'Color', 'w', 'FontSize', 16);
pbaspect(ax_em, [1, lmax/L, 1]);
set(ax_em, 'color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 16, 'Box', 'on', ...
    'CLim', [0, cmax_m], 'TickDir', 'out');

% Middle-right: exact psi_+ (cyan colormap)
ax_ep = axes('Position', [0.54, 0.35, 0.41, 0.30]);
imagesc(ax_ep, x_ctrs, l_ctrs, Psi_coarse_p');
axis(ax_ep, 'xy');
colormap(ax_ep, cm_cyan);
cb = colorbar(ax_ep);  cb.Color = 'w';
xlabel(ax_ep, 'minus-end position x (V_g/\lambda)', 'Color', 'w', 'FontSize', 16);  ylabel(ax_ep, 'length l  (V_g/\lambda)', 'Color', 'w', 'FontSize', 16);
title(ax_ep, '\psi_+(x,l)  [exact]', 'Color', 'w', 'FontSize', 16);
pbaspect(ax_ep, [1, lmax/L, 1]);
set(ax_ep, 'color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 16, 'Box', 'on', ...
    'CLim', [0, cmax_p], 'TickDir', 'out');

% Characteristic arrow: (0,0) -> (-eps*5, 5);  pbaspect makes x and l have equal pixel/unit
hold(ax_ep, 'on');
quiver(ax_ep, 0, 0, -eps*5, 5, 0, 'Color', 'c', 'LineWidth', 1, 'MaxHeadSize', 0.3);
quiver(ax_ep, 5, 0, -eps*5, 5, 0, 'Color', 'c', 'LineWidth', 1, 'MaxHeadSize', 0.3);
quiver(ax_ep, -5, 0, -eps*5, 5, 0, 'Color', 'c', 'LineWidth',1, 'MaxHeadSize', 0.3);
ang_char = atan2d(1, -eps);   % screen angle in degrees
% "characteristic" tilted, to the left of arrow (~l=2.2)
text(ax_ep, -eps*2.5 - 5.5, 2.5, 'characteristics', ...
    'Color', 'c', 'FontSize', 16, 'Rotation', ang_char, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Characteristic arrow: (0,0) -> (eps*5, 5);  pbaspect makes x and l have equal pixel/unit
hold(ax_em, 'on');
quiver(ax_em, 0, 0, eps*5, 5, 0, 'Color', 'm', 'LineWidth', 1, 'MaxHeadSize', 0.3);
quiver(ax_em, 5, 0, eps*5, 5, 0, 'Color', 'm', 'LineWidth', 1, 'MaxHeadSize', 0.3);
quiver(ax_em, -5, 0, eps*5, 5, 0, 'Color', 'm', 'LineWidth',1, 'MaxHeadSize', 0.3);
ang_char = atan2d(1, eps);   % screen angle in degrees
% "characteristic" tilted, to the left of arrow (~l=2.2)
text(ax_em, eps*2.5 - 5.5, 2.5, 'characteristics', ...
    'Color', 'm', 'FontSize', 16, 'Rotation', ang_char, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Bottom-left: histogram psi_-
ax_hm = axes('Position', [0.07, 0.03, 0.41, 0.30]);
imagesc(ax_hm, x_ctrs, l_ctrs, H_density_m');
axis(ax_hm, 'xy');
colormap(ax_hm, cm_magenta);
cb = colorbar(ax_hm);  cb.Color = 'w';
xlabel(ax_hm, 'minus-end position x  (V_g/\lambda)', 'Color', 'w', 'FontSize', 16);  ylabel(ax_hm, 'length l  (V_g/\lambda)', 'Color', 'w', 'FontSize', 16);
title(ax_hm, sprintf('\\psi_-(x,l)  [sampled]', Nf_hist), 'Color', 'w', 'FontSize', 16);
pbaspect(ax_hm, [1, lmax/L, 1]);
set(ax_hm, 'color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 16, 'Box', 'on', ...
    'CLim', [0, cmax_m], 'TickDir', 'out');

% Bottom-right: histogram psi_+
ax_hp = axes('Position', [0.54, 0.03, 0.41, 0.30]);
imagesc(ax_hp, x_ctrs, l_ctrs, H_density_p');
axis(ax_hp, 'xy');
colormap(ax_hp, cm_cyan);
cb = colorbar(ax_hp);  cb.Color = 'w';
xlabel(ax_hp, 'minus-end position x  (V_g/\lambda)', 'Color', 'w', 'FontSize', 16);  ylabel(ax_hp, 'length l  (V_g/\lambda)', 'Color', 'w', 'FontSize', 16);
title(ax_hp, sprintf('\\psi_+(x,l)  [sampled]', Nf_hist), 'Color', 'w', 'FontSize', 16);
pbaspect(ax_hp, [1, lmax/L, 1]);
set(ax_hp, 'color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 16, 'Box', 'on', ...
    'CLim', [0, cmax_p], 'TickDir', 'out');

drawnow;
exportgraphics(fig_psi, fullfile(outDir, sprintf('sample_M%d_eps%.3f.png', M, eps)), 'BackgroundColor', 'k', 'Resolution', 150);


%% ── Save steady-state snapshot WITH sampled rods ──────────────────────────
eps_str = strrep(sprintf('%.2f', eps), '.', 'p');
if exist('a_self', 'var') && a_self ~= 0
    aself_str = strrep(sprintf('%.2f', a_self), '.', 'p');
    filename  = fullfile(dataDir, sprintf('sim_M%d_eps%s_aself%s.mat', M, eps_str, aself_str));
else
    filename  = fullfile(dataDir, sprintf('sim_M%d_eps%s.mat', M, eps_str));
end

if exist('chain', 'var') && M > 0
    save(filename, 'x', 'L', 'eps', 'M', 'phi', 'eta', 's', 'c', 'pol', 'chain', ...
        'rod_x_p', 'rod_y_p', 'rod_x_m', 'rod_y_m', 'params');
else
    save(filename, 'x', 'L', 'eps', 'M', 'phi', 'eta', 's', 'c', 'pol', ...
        'rod_x_p', 'rod_y_p', 'rod_x_m', 'rod_y_m', 'params');
end

fprintf('Saved  %s\n', filename);