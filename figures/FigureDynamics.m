% =====================================================================
%  FigureDynamics.m  —  Manuscript Figures 2 (M=0) and 5 (M=5)
%  Transient relaxation to steady state: time snapshots x 2 columns.
%    Left  column : continuum density fields
%    Right column : stochastic filament snapshot at matching time (Gillespie)
%
%  Repository: spindle1D  (reproducibility code for the spindle-polarity paper)
%
%  Set matBase to select the run:
%    'M0_eps0.200'  -> Figure 2   (single-step)
%    'M5_eps0.200'  -> Figure 5   (multistep)
%
%  Input : ../data/fields_<matBase>.mat and nuc_t_<matBase>.mat
%          (time-series. The M0 pair (Figure 2) is shipped in data/; the larger
%           M5 pair (Figure 5) is not — regenerate it with PDEsolver.m.)
%  Output: ../output/FigureDynamics_<matBase>.pdf
%
%  Usage : set matBase and fi_list below, then run from anywhere.
% =====================================================================

clear; clc; close all;

% --- repository paths (resolved relative to this file) ---
thisDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(thisDir, '..', 'data');
outDir  = fullfile(thisDir, '..', 'output');
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% ── Settings ─────────────────────────────────────────────────────────────
matBase = 'M0_eps0.200';   % 'M0_eps0.200' -> Fig 2 ; 'M5_eps0.200' -> Fig 5

% Frame indices for the five time snapshots (1-based; t = (fi-1)*dt_save).
fi_list = [1, 51, 101, 1001, 2001];   % at dt_save=0.1: t = 0, 5, 10, 100, 200 — adjust to your run

% Ran-source signal width for the display overlay (matches PDEsolver.m ls)
ls_src = 0.5;

%% ── Load ─────────────────────────────────────────────────────────────────
fld = load(fullfile(dataDir, sprintf('fields_%s.mat', matBase)));
nuc = load(fullfile(dataDir, sprintf('nuc_t_%s.mat',  matBase)));

x       = fld.x(:);
L       = fld.L;
Np      = fld.N;
eps_val = fld.eps;
dt_sv   = fld.dt_save;
n_fr    = size(fld.phi_t, 3);
M_val   = fld.M;
has_chain = isfield(fld, 'chain_t') && ~isempty(fld.chain_t) && M_val > 0;

fi_list = min(fi_list, n_fr);
t_list  = (fi_list - 1) * dt_sv;

s    = exp(-(1/ls_src) * sqrt(x.^2));
dx_g = x(2) - x(1);

% c normalised to 1 at steady state
c_ref   = max(fld.c_t(:,   fi_list(end)));  c_ref   = max(c_ref,   1e-9);
phi_ref = max(fld.phi_t(:,:,fi_list(end)), [], 'all');  phi_ref = max(phi_ref, 1e-9);
eta_ref = max(fld.eta_t(:,:,fi_list(end)), [], 'all');  eta_ref = max(eta_ref, 1e-9);

% Per-stage chain reference: max of each stage at final frame
chain_refs = [];
if has_chain
    chain_refs = zeros(M_val, 1);
    for m = 1:M_val
        chain_refs(m) = max(fld.chain_t(m, :, :, fi_list(end)), [], 'all');
        chain_refs(m) = max(chain_refs(m), 1e-9);
    end
end

%% ── IC filaments from ψ(x,l) at frame 1 ────────────────────────────────
% Infer an effective birth field from the plus-end density eta at frame 1
% using the steady-state spectral relation below, then reconstruct ψ_±(x,l)
% and sample an initial filament population. This uses eta, not nuc_t(:,:,1):
% the default nonzero initial eta supplies filaments even when the actual
% initial nucleation rate is zero. In the formulas below, n_± denotes this
% inferred effective birth field.
%
% Formula (Vg = λ = 1):
%   ψ_+(x,l) = n_+(x + ε·l) · exp(−l)
%   ψ_−(x,l) = n_−(x − ε·l) · exp(−l)
% implemented spectrally via Fourier shift theorem.

rng(42);
y_range = [-1, 1];

lmax   = 6;
Nl_ic  = 400;
l_ic   = linspace(0, lmax, Nl_ic);
dl_ic  = l_ic(2) - l_ic(1);

k_ic   = (2*pi/L) * [0:(Np/2-1), (-Np/2):(-1)]';   % N×1

% Recover effective birth rate from η at frame 1 via steady-state spectral relation:
%   (1−ε)∂_x η_± + η_± = n_±  →  n̂_+(k) = η̂_+(k)·(1 + ik(1−ε))
%                                  n̂_−(k) = η̂_−(k)·(1 − ik(1−ε))
hat_eta_ic  = fft(fld.eta_t(:, :, 1));   % N×2
hat_nuc_ic  = zeros(Np, 2);
hat_nuc_ic(:,1) = hat_eta_ic(:,1) .* (1 + 1i*k_ic*(1-eps_val));
hat_nuc_ic(:,2) = hat_eta_ic(:,2) .* (1 - 1i*k_ic*(1-eps_val));
nuc_ic_real = real(ifft(hat_nuc_ic));    % N×2  for expected-count calculation

Psi_ic = zeros(Np, Nl_ic, 2);
for il = 1:Nl_ic
    decay = exp(-l_ic(il));
    Psi_ic(:,il,1) = real(ifft(hat_nuc_ic(:,1) .* exp( 1i*k_ic*eps_val*l_ic(il)))) * decay;
    Psi_ic(:,il,2) = real(ifft(hat_nuc_ic(:,2) .* exp(-1i*k_ic*eps_val*l_ic(il)))) * decay;
end

% Expected alive counts: ∫ψ dx dl = R_± · ⟨lifetime⟩ = R_± · 1
Nf_p_ic = ceil(sum(max(nuc_ic_real(:,1), 0)) * dx_g);
Nf_m_ic = ceil(sum(max(nuc_ic_real(:,2), 0)) * dx_g);

w_p_ic = max(Psi_ic(:,:,1), 0);  w_p_ic = w_p_ic(:) / (sum(w_p_ic(:)) + 1e-12);
w_m_ic = max(Psi_ic(:,:,2), 0);  w_m_ic = w_m_ic(:) / (sum(w_m_ic(:)) + 1e-12);

% Sample (xm, l) pairs from ψ  — guard against all-zero weights
if Nf_p_ic > 0 && any(w_p_ic > 0)
    idx_p    = randsample(Np*Nl_ic, Nf_p_ic, true, w_p_ic);
    [ix_p, il_p] = ind2sub([Np, Nl_ic], idx_p);
    xm_p_ic  = x(ix_p)';   lp_ic = l_ic(il_p)';
else
    xm_p_ic = [];  lp_ic = [];  Nf_p_ic = 0;
end

if Nf_m_ic > 0 && any(w_m_ic > 0)
    idx_m    = randsample(Np*Nl_ic, Nf_m_ic, true, w_m_ic);
    [ix_m, il_m] = ind2sub([Np, Nl_ic], idx_m);
    xm_m_ic  = x(ix_m)';   lm_ic = l_ic(il_m)';
else
    xm_m_ic = [];  lm_ic = [];  Nf_m_ic = 0;
end

% Convert to running-filament parameterisation used by the Gillespie loop:
%   age = l  →  t_birth = -l
%   x_birth_right = xm + ε·l   (so xm = x_birth - ε·age at render time)
%   x_birth_left  = xm - ε·l
%   remaining lifetime drawn from Exp(1)  (memoryless property)
t_births = [-lp_ic(:)',             -lm_ic(:)'];
t_deaths = [exprnd(1, [1,Nf_p_ic]), exprnd(1, [1,Nf_m_ic])];
x_birth  = [xm_p_ic(:)' + eps_val*lp_ic(:)', ...
            xm_m_ic(:)' - eps_val*lm_ic(:)'];
y_birth  = [y_range(1) + diff(y_range)*rand(1,Nf_p_ic), ...
            y_range(1) + diff(y_range)*rand(1,Nf_m_ic)];
pol      = [ones(1,Nf_p_ic), 2*ones(1,Nf_m_ic)];

fprintf('IC: %d right + %d left filaments sampled from psi at frame 1.\n', Nf_p_ic, Nf_m_ic);

snap   = cell(numel(fi_list), 1);
fi_ptr = 1;

for fi = 1 : max(fi_list)
    t_now   = (fi - 1) * dt_sv;
    nuc_now = nuc.nuc_t(:, :, fi);

    % Snapshot at each requested frame
    if fi_ptr <= numel(fi_list) && fi == fi_list(fi_ptr)
        alive = find(t_births <= t_now & t_deaths > t_now);
        snap{fi_ptr} = makeSnap(alive, t_now, eps_val, ...
                                x_birth, y_birth, pol, t_births);
        fi_ptr = fi_ptr + 1;
    end
    if fi == max(fi_list), break; end

    % Sample births over [t_now, t_now + dt_sv)
    R_p   = sum(nuc_now(:,1)) * dx_g;
    R_m   = sum(nuc_now(:,2)) * dx_g;
    R_tot = R_p + R_m;
    if R_tot > 1e-12
        w_p = max(nuc_now(:,1), 0);  w_p = w_p / (sum(w_p) + 1e-12);
        w_m = max(nuc_now(:,2), 0);  w_m = w_m / (sum(w_m) + 1e-12);
        t_rem = dt_sv;  t_ev = t_now;
        while true
            delt = exprnd(1 / R_tot);
            if delt > t_rem, break; end
            t_rem = t_rem - delt;
            t_ev  = t_ev  + delt;
            pm = 1 + (rand >= R_p / R_tot);
            if pm == 1
                xb = x(randsample(Np, 1, true, w_p));
            else
                xb = x(randsample(Np, 1, true, w_m));
            end
            t_births(end+1) = t_ev;                              %#ok<AGROW>
            t_deaths(end+1) = t_ev + exprnd(1);                  %#ok<AGROW>
            x_birth(end+1)  = xb;                                %#ok<AGROW>
            y_birth(end+1)  = y_range(1) + diff(y_range)*rand;   %#ok<AGROW>
            pol(end+1)      = pm;                                 %#ok<AGROW>
        end
    end

    keep     = t_deaths > t_now;
    t_births = t_births(keep);  t_deaths = t_deaths(keep);
    x_birth  = x_birth(keep);   y_birth  = y_birth(keep);
    pol      = pol(keep);
end
fprintf('Gillespie done (%d snapshots taken).\n', numel(fi_list));

%% ── Figure layout ────────────────────────────────────────────────────────
fig_w_in = 170/25.4;   % 6.693 in  (PNAS two-column)

ml_in = 0.55;  mr_in = 0.10;  mb_in = 0.38;  mt_in = 0.32;
gc_in = 0.18;  gr_in = 0.10;

n_rows   = numel(fi_list);
pw_in    = (fig_w_in - ml_in - mr_in - gc_in) / 2;
ph_in    = 0.35 * pw_in;   % aspect ratio height/width = 0.35
fig_h_in = mb_in + mt_in + n_rows * ph_in + (n_rows-1) * gr_in;

ml = ml_in / fig_w_in;   mb = mb_in / fig_h_in;
gc = gc_in / fig_w_in;   gr = gr_in / fig_h_in;
pw = pw_in / fig_w_in;   ph = ph_in / fig_h_in;
xl = [ml, ml + pw + gc];   % left edge of each column

fig_d = figure('Color', 'w');
set(fig_d, 'Units', 'inches', 'Position', [0 0 fig_w_in fig_h_in]);
set(fig_d, 'PaperUnits', 'inches', 'PaperSize', [fig_w_in fig_h_in], ...
    'PaperPosition', [0 0 fig_w_in fig_h_in]);

fs     = 10;

% Column headers
hdr_y = mb + n_rows*(ph + gr) - gr + 0.004;
hdr_h = mt_in/fig_h_in * 0.80;
 annotation(fig_d, 'textbox', [xl(1), hdr_y, pw, hdr_h], ...
     'String', 'Continuum fields', 'FontSize', fs+1, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'EdgeColor', 'none', 'FitBoxToText', 'off', 'Interpreter', 'none');
 annotation(fig_d, 'textbox', [xl(2), hdr_y, pw, hdr_h], ...
     'String', 'Stochastic filaments', 'FontSize', fs+1, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'EdgeColor', 'none', 'FitBoxToText', 'off', 'Interpreter', 'none');

%% ── Draw panels ─────────────────────────────────────────────────────────
for ri = 1:n_rows
    row_y  = mb + (n_rows - ri) * (ph + gr);   % ri=1 → top, ri=4 → bottom
    fi_r   = fi_list(ri);
    t_r    = t_list(ri);
    is_bot = (ri == n_rows);

    phi_r = fld.phi_t(:, :, fi_r);
    eta_r = fld.eta_t(:, :, fi_r);
    c_r   = fld.c_t(:,      fi_r);
    chain_r = [];
    if has_chain
        chain_r = fld.chain_t(:, :, :, fi_r);   % M x N x 2
    end

    % Left panel: fields
    ax_L = axes('Position', [xl(1), row_y, pw, ph]);  %#ok<LAXES>
    plotFieldsSnap(ax_L, x, L, s, phi_r, eta_r, c_r, chain_r, chain_refs, phi_ref, eta_ref, c_ref, is_bot, fs);

    % Right panel: stochastic filaments
    ax_R = axes('Position', [xl(2), row_y, pw, ph]);  %#ok<LAXES>
    plotFilamSnap(ax_R, snap{ri}, x, L, s, is_bot, fs);

    % Time stamp  (inside left panel, top-right in data coords)
    text(ax_L, 10 - 0.4, 1.2, sprintf('$t = %.0f\\,\\lambda^{-1}$', t_r), ...
        'Interpreter', 'latex', 'FontSize', fs, 'Color', 'k', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
end

%% ── Shared y-axis label (left / continuum-fields column) ────────────────
col_h  = n_rows*ph + (n_rows-1)*gr;      % normalised height of left column
yl_ctr = mb + col_h/2;                   % vertical centre of the column
xl_pos = 0.13 / fig_w_in;                % ~0.13 in from the figure's left edge
ax_yl  = axes('Position', [0 0 1 1], 'Visible', 'off', 'HitTest', 'off');
set(ax_yl, 'XLim', [0 1], 'YLim', [0 1]);
text(ax_yl, xl_pos, yl_ctr, 'normalized amplitudes', ...
    'Rotation', 90, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', fs, 'Interpreter', 'none');

%% ── Export ──────────────────────────────────────────────────────────────
drawnow;
out_file = fullfile(outDir, sprintf('FigureDynamics_%s.pdf', matBase));
exportgraphics(fig_d, out_file, 'ContentType', 'vector');
fprintf('Saved  %s\n', out_file);


%% ════════════════════════════════════════════════════════════════════════
%%  Local functions
%% ════════════════════════════════════════════════════════════════════════

function sn = makeSnap(alive, t_now, eps_val, x_birth, y_birth, pol, t_births)
% Build rod endpoint arrays for alive filaments at t_now.
sn.rod_x_p = zeros(3, 0);  sn.rod_y_p = zeros(3, 0);
sn.rod_x_m = zeros(3, 0);  sn.rod_y_m = zeros(3, 0);
if isempty(alive), return; end

ages = t_now - t_births(alive);
is_r = pol(alive) == 1;

if any(is_r)
    ag  = ages(is_r);
    xb  = x_birth(alive(is_r));
    yb  = y_birth(alive(is_r));
    xm  = xb - eps_val * ag;
    xp  = xm + ag;
    Nr  = sum(is_r);
    sn.rod_x_p = [xm; xp; NaN(1,Nr)];
    sn.rod_y_p = [yb; yb; NaN(1,Nr)];
end

if any(~is_r)
    ag  = ages(~is_r);
    xb  = x_birth(alive(~is_r));
    yb  = y_birth(alive(~is_r));
    xm  = xb + eps_val * ag;
    xp  = xm - ag;
    Nl  = sum(~is_r);
    sn.rod_x_m = [xp; xm; NaN(1,Nl)];   % row 1 = plus end (for '<' arrow)
    sn.rod_y_m = [yb; yb; NaN(1,Nl)];
end
end


% ─────────────────────────────────────────────────────────────────────────
function plotFieldsSnap(ax, x, L, s, phi, eta, c, chain_r, chain_refs, phi_ref, eta_ref, c_ref, show_xl, fs)
clr_blue   = [0.15 0.40 0.90];
clr_red    = [0.85 0.15 0.15];
clr_teal   = [0.00 0.62 0.62];
clr_purple = [0.50 0.00 0.75];
clr_teal_d = [0.00 0.35 0.35];
clr_teal_l = [0.55 0.88 0.88];
clr_purp_d = [0.28 0.00 0.42];
clr_purp_l = [0.78 0.55 0.90];

xv = x(:);  xr = xv';

s_n  = s(:) / max(s(:));
c_n  = c(:);% / c_ref;
p1 = phi(:,1) / phi_ref;   p2 = phi(:,2) / phi_ref;
e1 = eta(:,1) / eta_ref;   e2 = eta(:,2) / eta_ref;

ymax = 1.28;

hold(ax, 'on');
fill(ax, [xr, fliplr(xr)], [s_n', zeros(1,numel(s_n))], ...
    clr_blue, 'FaceAlpha', 0.40, 'EdgeColor', 'none');
plot(ax, xv, s_n, 'Color', clr_blue,   'LineWidth', 1.0);
plot(ax, xv, c_n, 'Color', clr_red,    'LineWidth', 1.0);
plot(ax, xv, e1,  'Color', clr_teal,   'LineStyle', '--', 'LineWidth', 1.5);
plot(ax, xv, e2,  'Color', clr_purple, 'LineStyle', '--', 'LineWidth', 1.5);
% Faint solid lines for each chain stage
M_ch = size(chain_r, 1);
for m = 1:M_ch
    ch_p = squeeze(chain_r(m,:,1))' / chain_refs(m);
    ch_m = squeeze(chain_r(m,:,2))' / chain_refs(m);
    plot(ax, xv, ch_p, 'Color', [clr_teal,   0.22], 'LineStyle', '-', 'LineWidth', 1);
    plot(ax, xv, ch_m, 'Color', [clr_purple, 0.22], 'LineStyle', '-', 'LineWidth', 1);
end
plot(ax, xv, p1,  'Color', clr_teal,   'LineStyle', '-',  'LineWidth', 2.0);
plot(ax, xv, p2,  'Color', clr_purple, 'LineStyle', '-',  'LineWidth', 2.0);

% Direction arrows: shown only when spindle is spatially structured
if std(p1) > 0.08
    [~, ip] = max(e1);  xpk_p = xv(ip);
    [~, im] = max(e2);  xpk_m = xv(im);
    hw = 1;  ya = 1.14;  dy = ymax * 0.013;
    patch(ax, 'Faces', [1 2 3 4], ...
        'Vertices', [xpk_p-hw, ya-dy; xpk_p+0*hw, ya-dy; xpk_p+0*hw, ya+dy; xpk_p-hw, ya+dy], ...
        'FaceVertexCData', [clr_teal_d; clr_teal_l; clr_teal_l; clr_teal_d], ...
        'FaceColor', 'interp', 'EdgeColor', 'none');
    plot(ax, xpk_p+0*hw, ya, '->', 'MarkerSize', 7, ...
        'MarkerFaceColor', clr_teal_l, 'MarkerEdgeColor', 'none');
    patch(ax, 'Faces', [1 2 3 4], ...
        'Vertices', [xpk_m-0*hw, ya-dy; xpk_m+hw, ya-dy; xpk_m+hw, ya+dy; xpk_m-0*hw, ya+dy], ...
        'FaceVertexCData', [clr_purp_l; clr_purp_d; clr_purp_d; clr_purp_l], ...
        'FaceColor', 'interp', 'EdgeColor', 'none');
    plot(ax, xpk_m-0*hw, ya, '<-', 'MarkerSize', 7, ...
        'MarkerFaceColor', clr_purp_l, 'MarkerEdgeColor', 'none');
end

xline(ax, 0, 'k:', 'LineWidth', 0.8);
hold(ax, 'off');

xlim(ax, [-10, 10]);%[-L/2, L/2]);
ylim(ax, [0, ymax]);
set(ax, 'FontSize', fs, 'Box', 'on', 'LineWidth', 1, 'TickDir', 'out', ...
    'YTick', [0 0.5 1.0], 'YTickLabel', {'0','0.5','1'});
% (shared y-label added once in main script, spanning the whole column)
if show_xl
    xlabel(ax, '$x\ (V_{\mathrm{g}}/\lambda)$', 'FontSize', fs, 'Interpreter', 'latex');
else
    set(ax, 'XTickLabel', {});
end
end


% ─────────────────────────────────────────────────────────────────────────
function plotFilamSnap(ax, sn, x, L, s, show_xl, fs)
clr_teal_d = [0.00 0.35 0.35];
clr_teal_l = [0.55 0.88 0.88];
clr_purp_d = [0.28 0.00 0.42];
clr_purp_l = [0.78 0.55 0.90];

Nf_p = size(sn.rod_x_p, 2);
Nf_m = size(sn.rod_x_m, 2);

hold(ax, 'on');

% Interleave teal and purple rods so neither population is always on top
N_draw = max(Nf_p, Nf_m);
for k = 1:N_draw
    if k <= Nf_p
        gradRods(ax, sn.rod_x_p(:,k), sn.rod_y_p(:,k), 1, 2, clr_teal_d, clr_teal_l);
        plot(ax, sn.rod_x_p(2,k), sn.rod_y_p(1,k), '>', ...
            'MarkerSize', 3, 'MarkerFaceColor', clr_teal_l, ...
            'MarkerEdgeColor', 'none', 'LineStyle', 'none');
    end
    if k <= Nf_m
        gradRods(ax, sn.rod_x_m(:,k), sn.rod_y_m(:,k), 2, 1, clr_purp_d, clr_purp_l);
        plot(ax, sn.rod_x_m(1,k), sn.rod_y_m(1,k), '<', ...
            'MarkerSize', 3, 'MarkerFaceColor', clr_purp_l, ...
            'MarkerEdgeColor', 'none', 'LineStyle', 'none');
    end
end

% Ran-source blue-haze overlay
s_n   = s(:)' / max(s(:));
xv_s  = x(:)';
clr_s = [0.15 0.40 0.90];
img_s = repmat(reshape(clr_s, 1, 1, 3), 2, numel(s_n), 1);
ih    = image(ax, 'XData', [xv_s(1), xv_s(end)], 'YData', [-1, 1], 'CData', img_s);
set(ih, 'AlphaData', repmat(s_n * 0.65, 2, 1));

hold(ax, 'off');

set(ax, 'Color', [0.94 0.94 0.94]);
xlim(ax, [-10, 10]);%[-L/2, L/2]);
ylim(ax, [-1, 1]);
set(ax, 'FontSize', fs, 'Box', 'on', 'LineWidth', 1, 'TickDir', 'out', 'YTick', []);
if show_xl
    xlabel(ax, '$x\ (V_{\mathrm{g}}/\lambda)$', 'FontSize', fs, 'Interpreter', 'latex');
else
    set(ax, 'XTickLabel', {});
end
end


% ─────────────────────────────────────────────────────────────────────────
function gradRods(ax, rx, ry, row_minus, row_plus, clr_dark, clr_light)
% Render rods as gradient-filled rectangles (dark minus end → light plus end).
% rx, ry : 3×Nf  (row 3 is NaN sentinel; only rows 1–2 are used for geometry)
Nf = size(rx, 2);
if Nf == 0, return; end
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
