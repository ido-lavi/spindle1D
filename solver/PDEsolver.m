% =====================================================================
%  PDEsolver.m  —  Continuum solver for the 1D spindle branching-nucleation model
%
%  Repository: spindle1D  (reproducibility code for the spindle-polarity paper)
%
%  Pseudo-spectral (integrating-factor AB2) integrator for the reduced-dimension
%  kinetic system: MT plus-end density, tubulin lattice density, M sequential
%  branching-complex stages, and the Ran-activated nucleator c.
%
%  Set the parameters in the SETUP block below (M, eps = V||/Vg, a_self, ...),
%  then run. Writes to ../data:
%    fields_<tag>.mat  — field time-series        (input to FigureDynamics.m)
%    nuc_t_<tag>.mat   — nucleation-rate series   (input to filamentMovie.m)
%  where <tag> = 'M<M>_eps<eps>'  (with '_aself<a_self>' appended when a_self > 0).
%  Diagnostic GIFs are written to ../output only when makeGif = true (off by default).
%
%  Steady-state snapshots for the figure scripts are then saved by running
%  fieldSteadyState.m  or  filamentSteadyState.m  in the same workspace.
%
%  Usage : run from anywhere (paths resolve relative to this file).
% =====================================================================

%% SETUP
clear; clc
close all;

% --- repository paths (resolved relative to this file) ---
thisDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(thisDir, '..', 'data');
outDir  = fullfile(thisDir, '..', 'output');
if ~exist(dataDir, 'dir'); mkdir(dataDir); end
if ~exist(outDir,  'dir'); mkdir(outDir);  end

% Physical dimensionless params
L     = 40;      % system size (units of vg/lambda)
eps   = 0.2;     % V_parallel/V_g
alpha = 20;      % alpha: binding affinity of c to the mature (stage-M) complex
nb    = 270;     % n_b: nucleation prefactor (units of lambda)
ls    = 0.5;     % l_s: RanGTP signal decay length (units of vg/lambda)
M     = 0;       % chain steps beyond phi (0: nucleation from phi)
a_m   = 0.1;     % a_m: per-stage association rate (chain step m -> m+1, units of lambda)
a0    = 1000;    % a_0: lattice-site production rate by plus ends (units of lambda)
a_self = 0.0;    % a_self: cooperative self-recruitment rate (units of lambda)
lambda_c = 0.2;  % lambda_c: nucleator degradation rate (units of lambda)
D_c      = 7.5;  % D: nucleator diffusion coefficient (units of ell^2*lambda = V_g^2/lambda)
s_0 = 0.01;      % s_0: RanGTP activation half-saturation
a_c = 0.2;         % a_c: RanGTP activation plateau (units of lambda)

% Rescaling alpha to fix loop gain: the entered alpha is a reference value, rescaled
% by 1/(nb*a0*a_m^M), with nb controlling the overall filament-count scale. The rescaled value
% (not the entered number) is the affinity actually used by the solver.
alpha = alpha/nb/a0/a_m^M;

% Numerical params
N      = 1024;
dt     = 0.02;
Tfinal = 200;
Nt     = round(Tfinal / dt);

% Diagnostics: write the relaxation GIFs to ../output (slow rendering; purely visual,
% not read by any figure script). Set true to regenerate the animations.
makeGif = false;

% Provenance: record all run parameters (saved into every snapshot below).
% alpha here is the rescaled value actually used by the solver (see rescale above).
params = struct('M', M, 'eps', eps, 'a_self', a_self, ...
                'a_m', a_m, 'a0', a0, 'nb', nb, 'ls', ls, 'D_c', D_c, ...
                'lambda_c', lambda_c, 's_0', s_0, 'a_c', a_c, 'alpha', alpha, ...
                'L', L, 'N', N, 'dt', dt, 'Tfinal', Tfinal);

% Discretization
x  = (-(N/2):(N/2-1)) * L/N;              % 1×N periodic grid
dx = x(2) - x(1);
k  = (2*pi/L) * [0:(N/2-1), (-N/2):(-1)];
del = 1i*k'; del(N/2+1) = 0;              % N×1, Nyquist zeroed
delsq = -(k.^2)'; 
DA  = [ones(1+floor(N/3),1); zeros(N-1-2*floor(N/3),1); ones(floor(N/3),1)]; % N×1

% Enslaving denominators 
denom_pm = [1 - eps*del,     1 + eps*del];
% Integrating factors (N×2: col 1=plus, col 2=minus)
IFeta    = [exp(((eps-1)*del - 1)*dt),  exp((-(eps-1)*del - 1)*dt)];
IF0      = [exp(( eps*del    - 1)*dt),  exp((-eps*del    - 1)*dt)];
IFm      = [exp(( eps*del + a_self - 1)*dt),  exp((-eps*del + a_self - 1)*dt)];
% Integrating factor c
IFc      = exp((D_c*delsq -lambda_c)*dt);

% Nucleation signal
s     = exp(-(1/ls) * sqrt(x.^2))';           % N×1
sact  = a_c* s./(s_0+s);
hat_s = fft(sact); 

% Initial conditions (N×2, col 1 = plus, col 2 = minus)
xi  = 5*ones(N,2);
eta = 5*ones(N,2);
phi = zeros(N,2);
% Initial condition c
c    = zeros(N,1);

% Fourier ICs
hat_xi = fft(xi);
hat_eta = fft(eta);
hat_phi = fft(phi);
% Fourier IC c
hat_c   = fft(c);

% Chain: M lattice-bound steps beyond phi  (M×N×2, [] when M=0)
if M > 0
    hat_chain = zeros(M, N, 2);
else
    hat_chain = [];
end

%% Enslave ICs
hat_phi = a0 * hat_eta ./ denom_pm;
for m = 1:M
    if m == 1
        prev = hat_phi;
    else
        prev = squeeze(hat_chain(m-1,:,:));
    end
    hat_chain(m,:,:) = a_m * reshape(prev ./ denom_pm, 1, N, 2);
end

% Time-stepping current state
hat_Xi   = hat_xi;
hat_Eta  = hat_eta;
hat_Phi  = hat_phi;
hat_Chain = hat_chain;   % M×N×2 (or [])
% Time-stepping current state c
hat_C    = hat_c;

% AB2 history (previous step)
Xi_old   = hat_xi;
Eta_old  = hat_eta;
Phi_old  = hat_phi;
Chain_old = hat_chain;   % M×N×2 (or [])
nuc_old   = zeros(N,2);
% AB2 c transport
cTransport_old = zeros(N,1);

% Centroid init (for plot arrows)
cm_val = 0;

% Colors and plot scale
clr.purple = [0.5, 0, 0.8];
clr.teal   = [0, 0.6, 0.6];
clr.mag    = [0.90, 0, 0.90];
clr.blue   = [81,80,195]/255;
clr.orange = [218,163,0]/255;
clr.olive  = [0.467, 0.576, 0.235];
clr.gray = [.7 .7 .7];
maxphi = 20;

% Chain step colors
if M > 0
    palette    = [clr.orange; clr.blue; clr.mag; clr.teal; clr.purple];      % 5
    chain_clrs = palette(mod((0:M-1)', 5) + 1, :);    % M×5
else
    chain_clrs = [];
end

%% Initial plot (real-space fields at ICs, before loop)
rho = phi(:,1) + phi(:,2);            % zeros initially
pol = (phi(:,1) - phi(:,2)) ./ rho;  % NaN initially (matches original)
chain_r = zeros(M, N, 2);               % zeros initially ([] when M=0 is fine too)
if M == 0, chain_r = []; end

if makeGif
    f = figure; f.Position = [100 10 1080 720];
    set(gcf, 'color', 'white');
    doPlot(x, L, s, c, xi, eta, phi, chain_r, pol, rho, cm_val, maxphi, clr, chain_clrs, eps, nb, a0, a_m);

    gifBase  = fullfile(outDir, sprintf('PDEsolver_M%d_eps%.3f', M, eps));
    gifFirst = [gifBase '_first.gif'];
    gifFile  = [gifBase '.gif'];
    gifLast  = [gifBase '_last.gif'];
    frame = getframe(1);
    im = frame2im(frame);
    [imind, cmap] = rgb2ind(im, 256);
    imwrite(imind, cmap, gifFirst, 'gif', 'DelayTime', 0.04, 'Loopcount', 0);
    imwrite(imind, cmap, gifFile,  'gif', 'DelayTime', 0.04, 'Loopcount', 0);
end

%% nuc_t accumulator (for filamentMovie.m)
dt_save         = 0.1;
n_save_interval = round(dt_save / dt);
n_save_total    = floor(Nt / n_save_interval) + 1;
nuc_t           = zeros(N, 2, n_save_total);
nuc_t_idx       = 0;

%% Field snapshot accumulators (for FigureDynamics.m)
phi_t = zeros(N, 2, n_save_total);
eta_t = zeros(N, 2, n_save_total);
c_t   = zeros(N,    n_save_total);
if M > 0
    chain_t = zeros(M, N, 2, n_save_total);
end

%% TIME STEPPING
for n = 0:Nt
    t = n * dt;

    %% Real-space fields
    xi  = real(ifft(hat_xi));
    eta = real(ifft(hat_eta));
    phi = real(ifft(hat_phi));
    chain = zeros(M, N, 2);
    for m = 1:M
        chain(m,:,:) = reshape(real(ifft(squeeze(hat_chain(m,:,:)))), 1, N, 2);
    end
    % Real-space c (antialiased since c eq is nonlinear)
    c   = real(ifft(DA.*hat_c));

    %% Nucleation driver: last chain step, or phi if M=0. Anti-Alias since nucleation is nonlinear.
    if M == 0
        phi_M = real(ifft(DA.*hat_phi));
    else
        phi_M = real(ifft(DA.*squeeze(hat_chain(M,:,:))));   % N×2
    end

    %% Nucleation rate
    denom   = 1 + alpha * (phi_M(:,1) + phi_M(:,2));
    cu      = c ./ denom;                   % unbound c (N×1)
    cb      = alpha * phi_M .* cu;          % bound c on + and - MTs (N×2)
    nuc     = nb * cb;                      % nucleation rate
    nuc     = real(ifft(DA.*fft(nuc)));
    % Fourier
    hat_nuc = fft(nuc);
    if n == 0
        nuc_old = hat_nuc;
    end

    %% Save nuc + field snapshots
    if mod(n, n_save_interval) == 0
        nuc_t_idx = nuc_t_idx + 1;
        nuc_t(:,:,nuc_t_idx) = nuc;
        phi_t(:,:,nuc_t_idx) = phi;
        eta_t(:,:,nuc_t_idx) = eta;
        c_t(:,  nuc_t_idx)   = c;
        if M > 0
            chain_t(:,:,:,nuc_t_idx) = chain;
        end
    end

    %% IF-AB2 updates
    % xi (minus-end density, driven by nuc, uses IF0)
    hat_xi = IF0 .* (hat_Xi + (dt/2) * (3*hat_nuc  - IF0 .* nuc_old));
    % eta (plus-end density, driven by nuc, uses IFeta)
    hat_eta = IFeta .* (hat_Eta + (dt/2) * (3*hat_nuc  - IFeta .* nuc_old));
    % phi (MT lattice, driven by eta, uses IF0)
    hat_phi = IF0 .* (hat_Phi + (dt/2) * a0 * (3*hat_Eta - IF0 .* Eta_old));
    % chain steps (lattice-bound, each driven by previous, all use IFm)
    for m = 1:M
        if m == 1
            drv     = hat_Phi;
            drv_old = Phi_old;
        else
            drv     = squeeze(hat_Chain(m-1,:,:));
            drv_old = squeeze(Chain_old(m-1,:,:));
        end
        curr = squeeze(hat_Chain(m,:,:));
        hat_chain(m,:,:) = reshape(IFm .* (curr + (dt/2) * a_m * (3*drv - IFm.*drv_old)), 1, N, 2);
    end
    % c (Ran-activated species, total concentration: unbound + bound_+-)
    hat_advc     = eps * del .* fft(cb(:,1) - cb(:,2));    % advection rhs (Fourier)
    hat_diffc    = D_c * delsq.*fft(cu - c);               % diffusion rhs (Fourier)
    hat_cTransport = hat_advc + hat_diffc;          % current transport
    if n == 0
        cTransport_old = hat_cTransport;
    end
    hat_c = IFc .* (hat_C + (dt/2) * (3*hat_cTransport - IFc .* cTransport_old) + (dt/2) * (3 - IFc) .* hat_s);

    %% Update AB2 history
    Eta_old  = hat_Eta;
    Phi_old  = hat_Phi;
    Chain_old = hat_Chain;
    nuc_old   = hat_nuc;
    % AB2 c transport
    cTransport_old = hat_cTransport;

    %% Advance current state
    hat_Xi   = hat_xi;
    hat_Eta  = hat_eta;
    hat_Phi  = hat_phi;
    hat_Chain = hat_chain;
    % Advance current c
    hat_C    = hat_c;

    %% Density and polarity
    rho = phi(:,1) + phi(:,2);
    pol = (phi(:,1) - phi(:,2)) ./ rho;

    %% Centroid of the plus-end lattice (positions the plot arrows; every 2 steps)
    if makeGif && mod(n,2) == 0
        cm_val = sum(x' .* phi(:,1)) / sum(phi(:,1));
    end

    %% Plot every 8 steps (diagnostic GIF frames)
    if makeGif && mod(n,8) == 0 && n > 0
        doPlot(x, L, s, c, xi, eta, phi, chain, pol, rho, cm_val, maxphi, clr, chain_clrs, eps, nb, a0, a_m);
        frame = getframe(1);
        im = frame2im(frame);
        [imind, cmap] = rgb2ind(im, 256);
        imwrite(imind, cmap, gifFile, 'gif', 'DelayTime', 0.05, 'WriteMode', 'append');
    end
end

%% Save nuc_t for filamentMovie.m
% Tag: plain 'M<M>_eps<eps>' for non-cooperative runs; append '_aself<a_self>'
% only when a_self > 0, so the main figures load the plain-named files.
if a_self == 0
    tag = sprintf('M%d_eps%.3f', M, eps);
else
    tag = sprintf('M%d_eps%.3f_aself%.2f', M, eps, a_self);
end
nuc_t      = nuc_t(:,:,1:nuc_t_idx);
nuc_t_file = fullfile(dataDir, sprintf('nuc_t_%s.mat', tag));
save(nuc_t_file, 'nuc_t', 'x', 'eps', 'L', 'N', 'M', 'dt_save', 'params');
fprintf('Saved nuc_t (%d frames, dt=%.3f) to %s\n', nuc_t_idx, dt_save, nuc_t_file);

%% Save field snapshots for FigureDynamics.m
phi_t       = phi_t(:,:,1:nuc_t_idx);
eta_t       = eta_t(:,:,1:nuc_t_idx);
c_t         = c_t(:,    1:nuc_t_idx);
fields_file = fullfile(dataDir, sprintf('fields_%s.mat', tag));
if M > 0
    chain_t = chain_t(:,:,:,1:nuc_t_idx);
    save(fields_file, 'phi_t', 'eta_t', 'c_t', 'chain_t', 'x', 'eps', 'L', 'N', 'M', 'dt_save', 'params');
else
    save(fields_file, 'phi_t', 'eta_t', 'c_t', 'x', 'eps', 'L', 'N', 'M', 'dt_save', 'params');
end
fprintf('Saved fields (%d frames, dt=%.3f) to %s\n', nuc_t_idx, dt_save, fields_file);

%% -- Plot last (diagnostic GIF) ---
if makeGif
    doPlot(x, L, s, c, xi, eta, phi, chain, pol, rho, cm_val, maxphi, clr, chain_clrs, eps, nb, a0, a_m);
    frame = getframe(1);
    im = frame2im(frame);
    [imind, cmap] = rgb2ind(im, 256);
    imwrite(imind, cmap, gifLast, 'gif', 'DelayTime', 0.05, 'Loopcount', 0);
end

%% ── Local function ───────────────────────────────────────────────────────────
function doPlot(x, L, s, c, xi, eta, phi, chain, pol, rho, cm_val, maxphi, clr, chain_clrs, eps, nb, a0, a_m)
    M = size(chain, 1);   % 0 when chain is []
    plot(x, s, 'Color', 'b', 'LineWidth', 1.6);
    hold on
    fill([x' fliplr(x')], [s zeros(size(s))], ...
     'b','FaceAlpha',0.3,'EdgeColor','none')
    ylim manual; ylim([-1.01 1.3]) %[-1, 1.4] is pol is plotted
    xlim manual; xlim([-L/2 L/2])
    xlabel('x (Vg/lambda)'); ylabel('AU');
    title(sprintf('M=%d model  (Vk/Vg=1/5)', M));
    ax = gca; ax.FontSize = 20;
    pbaspect([2.0 1 1]); %[1.7 1 1] i pol is plotted
    % c (Ran-activated, red)
    plot(x, c, 'Color', 'r', 'LineWidth', 1.6);
    % eta (MT plus ends, gray)
    plot(x, eta(:,1) * (10/nb),       'Color', clr.gray, 'LineWidth', 2.6);
    plot(x, eta(:,2) * (10/nb), '--', 'Color', clr.gray, 'LineWidth', 2.6);
    % phi (MT lattice, olive)
    plot(x, phi(:,1) * (10/nb) / a0,       'Color', clr.olive, 'LineWidth', 2.6);
    plot(x, phi(:,2) * (10/nb) / a0, '--', 'Color', clr.olive, 'LineWidth', 2.6);
    % chain steps (orange → magenta gradient)
    for m = 1:M
        col = chain_clrs(m,:);
        plot(x, squeeze(chain(m,:,1))' * (10/nb) / a0/a_m^m,       'Color', col, 'LineWidth', 2.6);
        plot(x, squeeze(chain(m,:,2))' * (10/nb) / a0/a_m^m, '--', 'Color', col, 'LineWidth', 2.6);
    end
    % polarization (black)
    plot(x, pol, 'Color', 'k', 'LineWidth', 2.6);
    % double-headed arrows at ±cm_val
    hw = 0.4; ya = 1.066;
    plot( cm_val+hw, ya, '->',                              'MarkerSize', 20, 'MarkerFaceColor', clr.olive, 'MarkerEdgeColor', 'none');
    plot([ cm_val-hw,  cm_val+hw], [ya, ya], '-', 'Color', clr.olive, 'LineWidth', 4);
    plot(-cm_val-hw, ya, '<-',                              'MarkerSize', 20, 'MarkerFaceColor', clr.olive, 'MarkerEdgeColor', 'none');
    plot([-cm_val-hw, -cm_val+hw], [ya, ya], '-', 'Color', clr.olive, 'LineWidth', 4);
    hold off
    drawnow;
end
