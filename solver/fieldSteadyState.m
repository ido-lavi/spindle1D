% =====================================================================
%  fieldSteadyState.m  —  save a steady-state snapshot (fields only, no rods)
%
%  Repository: spindle1D  (reproducibility code for the spindle-polarity paper)
%
%  Run in the SAME workspace right after PDEsolver.m has reached steady state.
%  Writes ../data/sim_M<M>_eps<eps>.mat  (with '_aself<a_self>' when a_self > 0).
%  Saved: x, L, eps, M, phi, eta, s, c, pol  (+ chain if M > 0).
%
%  Use this light saver for the phase-diagram grid (Figure7.m). For the
%  field+filament panels (Figures 3, 6, 9) use filamentSteadyState.m instead, which
%  also stores the sampled rods.
%
%  Usage:  PDEsolver.m  ->  fieldSteadyState.m
% =====================================================================

% --- data path (resolved relative to this file) ---
thisDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(thisDir, '..', 'data');
if ~exist(dataDir, 'dir'); mkdir(dataDir); end

eps_str = strrep(sprintf('%.2f', eps), '.', 'p');
if exist('a_self', 'var') && a_self ~= 0
    aself_str = strrep(sprintf('%.2f', a_self), '.', 'p');
    filename  = fullfile(dataDir, sprintf('sim_M%d_eps%s_aself%s.mat', M, eps_str, aself_str));
else
    filename  = fullfile(dataDir, sprintf('sim_M%d_eps%s.mat', M, eps_str));
end

if exist('chain', 'var') && M > 0
    save(filename, 'x', 'L', 'eps', 'M', 'phi', 'eta', 's', 'c', 'pol', 'chain', 'params');
else
    save(filename, 'x', 'L', 'eps', 'M', 'phi', 'eta', 's', 'c', 'pol', 'params');
end

fprintf('Saved  %s\n', filename);