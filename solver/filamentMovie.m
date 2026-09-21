% =====================================================================
%  filamentMovie.m  —  stochastic filament movie (time-dependent)
%
%  Repository: spindle1D  (reproducibility code for the spindle-polarity paper)
%
%  Renders the Supplementary Movies: an inhomogeneous-Poisson filament
%  simulation driven frame-by-frame by the resolved nucleation-rate history
%  nuc_t(x,t) from the continuum solver. Births are sampled from the
%  instantaneous nuc_t(x,t); exponential (Poisson) death times are drawn at
%  birth. While alive, each minus end moves at speed -eps and each plus end
%  at 1-eps.
%
%  Input : ../data/nuc_t_<...>.mat   (produced by PDEsolver.m; set matName below)
%  Output: ../output/filamentMovie_M<M>_eps<eps>.gif
%
%  Usage : set matName, then run from anywhere.
% =====================================================================

clear; clc;
close all;

% --- repository paths (resolved relative to this file) ---
thisDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(thisDir, '..', 'data');
outDir  = fullfile(thisDir, '..', 'output');
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% Load nuc_t data  —  set filename here
matName = 'nuc_t_M0_eps0.200.mat';

fprintf('Loading %s\n', matName);
load(fullfile(dataDir, matName), 'nuc_t', 'x', 'eps', 'L', 'N', 'M', 'dt_save');

%% Parameters
dt_frame  = dt_save;
y_range   = [-1, 1];
gifFile   = fullfile(outDir, sprintf('filamentMovie_M%d_eps%.3f.gif', M, eps));

dx       = x(2) - x(1);
n_frames = size(nuc_t, 3);
nmax_nuc = max(nuc_t(:)) * 1.1;

%% Initialise running filament arrays (empty: simulation starts from scratch)
t_births = [];
t_deaths = [];
x_birth  = [];
y_birth  = [];
polarity = [];   % 1 = right-pointing (cyan), 2 = left-pointing (magenta)

%% Build figure
fig_anim = figure;
fig_anim.Position = [50 100 1400 1000];
set(fig_anim, 'color', 'k');

ax_nuc = axes('Position', [0.05, 0.71, 0.90, 0.18]);
set(ax_nuc, 'color', 'k', 'XColor', 'w', 'YColor', 'w', ...
    'YAxisLocation', 'origin', 'FontSize', 16);

ax_xm = axes('Position', [0.05, 0.42, 0.90, 0.18]);
set(ax_xm, 'color', 'k', 'XColor', 'w', 'YColor', 'none');

ax_xp = axes('Position', [0.05, 0.24, 0.90, 0.18]);
set(ax_xp, 'color', 'k', 'XColor', 'w', 'YColor', 'none');

ax_all = axes('Position', [0.05, 0.06, 0.90, 0.18]);
set(ax_all, 'color', 'k', 'XColor', 'w', 'YColor', 'none');

linkaxes([ax_nuc, ax_xp, ax_xm, ax_all], 'x');

%% Animation loop
for fi = 1 : n_frames
    t_now   = (fi - 1) * dt_frame;
    nuc_now = nuc_t(:, :, fi);

    %--- 1. Render state at t_now ---
    alive = find(t_births <= t_now & t_deaths > t_now);

    cla(ax_nuc); hold(ax_nuc, 'on'); set(ax_nuc, 'color', 'k');
    plot(ax_nuc, x, nuc_now(:,1), 'Color', 'c', 'LineWidth', 2.5);
    plot(ax_nuc, x, nuc_now(:,2), 'Color', 'm', 'LineWidth', 2.5);
    xlim(ax_nuc, [-10, 10]);%[-L/2, L/2]);
    ylim(ax_nuc, [0, nmax_nuc]);
    title(ax_nuc, sprintf('nucleation rate   t = %.2f   (M=%d, %c=%g)', t_now, M, char(949), eps), ...
        'Color', 'w', 'FontSize', 20);
    xlabel(ax_nuc, 'x  (V_g/\lambda)', 'Color', 'w', 'FontSize', 16);
    pbaspect(ax_nuc, [1, 2*diff(y_range)/L, 1]);

    cla(ax_xp);  hold(ax_xp,  'on');  set(ax_xp,  'color', 'k');
    cla(ax_xm);  hold(ax_xm,  'on');  set(ax_xm,  'color', 'k');
    cla(ax_all); hold(ax_all, 'on');  set(ax_all, 'color', 'k');

    xm_r = []; xp_r = []; yy_r = []; Nr = 0;
    xm_l = []; xp_l = []; yy_l = []; Nl = 0;

    if ~isempty(alive)
        ages = t_now - t_births(alive);
        is_r = polarity(alive) == 1;

        if any(is_r) % right-pointing
            ag   = ages(is_r);
            xm_r = x_birth(alive(is_r)) - eps * ag;
            xp_r = xm_r + ag;
            yy_r = y_birth(alive(is_r));
            Nr   = sum(is_r);
            plot(ax_xm, xm_r, yy_r, 'o', 'Color', 'c', 'MarkerSize', 4, 'MarkerFaceColor', 'c', 'LineStyle', 'none');
            plot(ax_xp, xp_r, yy_r, '>', 'Color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'c', 'LineStyle', 'none');
        end

        if any(~is_r) % left-pointing
            ag   = ages(~is_r);
            xm_l = x_birth(alive(~is_r)) + eps * ag;
            xp_l = xm_l - ag;
            yy_l = y_birth(alive(~is_r));
            Nl   = sum(~is_r);
            plot(ax_xm, xm_l, yy_l, 'o', 'Color', 'm', 'MarkerSize', 4, 'MarkerFaceColor', 'm', 'LineStyle', 'none');
            plot(ax_xp, xp_l, yy_l, '<', 'Color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'm', 'LineStyle', 'none');
        end
    end

    % Interleaved draw on ax_all so neither population dominates
    N_draw = max(Nr, Nl);
    for kk = 1:N_draw
        if kk <= Nr
            plot(ax_all, [xm_r(kk) xp_r(kk)], [yy_r(kk) yy_r(kk)], 'Color', 'c', 'LineWidth', 1.8);
            plot(ax_all, xp_r(kk), yy_r(kk), '>', 'Color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'c', 'LineStyle', 'none');
            plot(ax_all, xm_r(kk), yy_r(kk), 'o', 'Color', 'c', 'MarkerSize', 4, 'MarkerFaceColor', 'c', 'LineStyle', 'none');
        end
        if kk <= Nl
            plot(ax_all, [xp_l(kk) xm_l(kk)], [yy_l(kk) yy_l(kk)], 'Color', 'm', 'LineWidth', 1.8);
            plot(ax_all, xp_l(kk), yy_l(kk), '<', 'Color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'm', 'LineStyle', 'none');
            plot(ax_all, xm_l(kk), yy_l(kk), 'o', 'Color', 'm', 'MarkerSize', 4, 'MarkerFaceColor', 'm', 'LineStyle', 'none');
        end
    end

    for ax = {ax_xm, ax_xp, ax_all}
        xlim(ax{1}, [-10, 10]);%[-L/2, L/2]);
        ylim(ax{1}, y_range + [-0.05 0.05]);
        pbaspect(ax{1}, [1, 2*diff(y_range)/L, 1]);
        set(ax{1}, 'FontSize', 16, 'YTick', [], 'YColor', 'none');
    end
    set(ax_xp,  'XTickLabel', []);
    set(ax_xm,  'XTickLabel', []);
    set(ax_all, 'XTickLabel', []);
    title(ax_xm, sprintf('N_{alive} = %d', numel(alive)), 'Color', 'w', 'FontSize', 20);
    text(ax_xp,  -9.8, 0.88, 'plus-ends',  'Color', 'w', 'FontSize', 20, 'VerticalAlignment', 'top');
    text(ax_xm,  -9.8, 0.88, 'minus-ends', 'Color', 'w', 'FontSize', 20, 'VerticalAlignment', 'top');
    text(ax_all, -9.8, 0.88, 'filaments',  'Color', 'w', 'FontSize', 20, 'VerticalAlignment', 'top');
    drawnow;

    frame_data = getframe(fig_anim);
    im = frame2im(frame_data);
    [imind, cm_gif] = rgb2ind(im, 256);
    if fi == 1
        imwrite(imind, cm_gif, gifFile, 'gif', 'DelayTime', dt_frame/3, 'Loopcount', inf);
    else
        imwrite(imind, cm_gif, gifFile, 'gif', 'DelayTime', dt_frame/3, 'WriteMode', 'append');
    end

    %--- 2. Sample births for the next interval [t_now, t_now + dt_frame) ---
    if fi < n_frames
        R_plus  = sum(nuc_now(:,1)) * dx;
        R_minus = sum(nuc_now(:,2)) * dx;
        R_nuc   = R_plus + R_minus;

        w_p = max(nuc_now(:,1), 0);  w_p = w_p / sum(w_p);
        w_m = max(nuc_now(:,2), 0);  w_m = w_m / sum(w_m);

        t_rem = dt_frame;
        t_ev  = t_now;
        while true
            tau   = exprnd(1 / R_nuc);
            if tau > t_rem, break; end
            t_rem = t_rem - tau;
            t_ev  = t_ev  + tau;

            pm = 1 + (rand >= R_plus / R_nuc);   % 1=right, 2=left
            if pm == 1
                xb = x(randsample(N, 1, true, w_p));
            else
                xb = x(randsample(N, 1, true, w_m));
            end

            t_births(end+1) = t_ev;                               
            t_deaths(end+1) = t_ev + exprnd(1);                   
            x_birth(end+1)  = xb;                                 
            y_birth(end+1)  = y_range(1) + diff(y_range) * rand;  
            polarity(end+1) = pm;                                 
        end
    end

    %--- 3. Cull filaments that died before t_now (will never be alive again) ---
    keep     = t_deaths > t_now;
    t_births = t_births(keep);
    t_deaths = t_deaths(keep);
    x_birth  = x_birth(keep);
    y_birth  = y_birth(keep);
    polarity = polarity(keep);
end
fprintf('Saved %d frames to %s\n', n_frames, gifFile);