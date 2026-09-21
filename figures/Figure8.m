% =====================================================================
%  Figure8.m  —  Manuscript Figure 8
%  Branching statistics in the multistep model. All panels ANALYTIC
%  (dimensionless, lambda = Vg = 1); no data files needed.
%
%  Repository: spindle1D  (reproducibility code for the spindle-polarity paper)
%
%  Panels
%    (a) Age PDF of branching events: Gamma(M+2,1) = tau^(M+1) e^{-tau}/Gamma(M+2);
%        dashed reference = bare MT survival factor e^{-tau}. Mean = M+2.
%    (b) Age CDF: regularized lower incomplete gamma, gammainc(tau, M+2).
%    (c) Mature-complex density vs normalized arclength u = s/l (u=0 minus end,
%        u=1 plus end): (M+1)(1-u)^M. Mean from minus end = 1/(M+2).
%    (d) Branch-position PDF at distance d from the plus end, marginalized over
%        the steady-state length distribution e^{-l}: d^M e^{-d}/Gamma(M+1). Mean = M+1.
%    Colored markers in (a),(c),(d) indicate each distribution's mean;
%    markers in (b) indicate the median age.
%
%  Input : none (analytic).
%  Output: ../output/Figure8.pdf
%
%  Usage : run from anywhere (paths resolve relative to this file).
% =====================================================================

clear; clc; close all;

% --- repository paths (resolved relative to this file) ---
thisDir = fileparts(mfilename('fullpath'));
outDir  = fullfile(thisDir, '..', 'output');
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% ---- style ----
fs    = 10;
Mlist = 0:4;
nM    = numel(Mlist);
clr   = interp1([0;1], [0.85 0.15 0.15; 0.20 0.45 0.90], linspace(0,1,nM)');  % red -> blue (M-sweep, as in Figure3/5)
clr_ref = [0 0 0];    % dashed survival reference (black, panels a & b)

%% ---- figure geometry (~14.5 cm wide, 2x2 panels; canvas scaled to 0.85 so panels shrink while fonts stay fixed) ----
w_in = 6.6929;  h_in = 5.2;
fig = figure('Color','w');
set(fig,'Units','inches','Position',[1 1 w_in h_in]);
set(fig,'PaperUnits','inches','PaperSize',[w_in h_in],'PaperPosition',[0 0 w_in h_in]);
drawnow;   % settle the figure at its target size before any text-width measurements

xL1 = 0.095; xL2 = 0.595; wp = 0.375;
yTop = 0.585; yBot = 0.110; hp = 0.355;
posA = [xL1 yTop wp hp];
posB = [xL2 yTop wp hp];
posC = [xL1 yBot wp hp];
posD = [xL2 yBot wp hp];

%% ---- grids ----
tau = linspace(0, 12, 700);
u   = linspace(0,  1, 700);
dd  = linspace(0, 12, 700);

%% ================= Panel (a): age PDF =================
axA = axes('Position', posA); hold(axA,'on');
xline(axA, 1, '|', 'Color',[0.3 0.3 0.3], 'LineWidth',0.2);
hRef = plot(axA, tau, exp(-tau), '--', 'Color', clr_ref, 'LineWidth', 1.8);
hC   = gobjects(1,numel(Mlist));
for i = 1:numel(Mlist)
    M = Mlist(i);
    pdf = tau.^(M+1).*exp(-tau)/gamma(M+2);
    hC(i) = plot(axA, tau, pdf, 'Color', clr(i,:), 'LineWidth', 1.8);
    mu  = M+2;                                  % mean of Gamma(M+2,1)
    plot(axA, mu, 0.82, 'o', 'MarkerSize',5, 'MarkerFaceColor',clr(i,:), ...   % mean marker, leveled at ~80% ymax
         'MarkerEdgeColor',clr(i,:)*0.6, 'LineWidth',0.8);
end
xlim(axA,[0 12]); ylim(axA,[0 1.02]);
text(axA, 4, 0.91, '$\langle\tau\rangle= (M+2)/\lambda$', 'Interpreter','latex', ...
     'FontSize', fs-1, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
% ---- manual legend: M column (sample line + label per M) ----
xM1 = 8.5; xM2 = 9.2; xMt = 9.4;      % M-column sample + label x
yTopL = 0.94; dyL = 0.105;            % top row (M=0) and row spacing
for i = 1:nM
    yy = yTopL - (i-1)*dyL;
    plot(axA, [xM1 xM2], [yy yy], '-', 'Color', clr(i,:), 'LineWidth', 1.8);
    text(axA, xMt, yy, sprintf('$M=%d$', Mlist(i)), 'Interpreter','latex', ...
         'FontSize', fs-1, 'VerticalAlignment','middle', 'HorizontalAlignment','left');
end
hold(axA,'off');
% x-label placed in the x-labels block at the end (measured on the finished figure)
ylabel(axA,'branching PDF','FontSize',fs);
styleAxis(axA, fs);

%% ================= Panel (c): arclength density =================
axB = axes('Position', posB); hold(axB,'on');
for i = 1:numel(Mlist)
    M  = Mlist(i);
    fb = (M+1)*(1-u).^M;
    plot(axB, u, fb, 'Color', clr(i,:), 'LineWidth', 1.8);
    mu  = 1/(M+2);                            % mean position from minus end
    plot(axB, mu, 4.24, 'o', 'MarkerSize',5, 'MarkerFaceColor',clr(i,:), ...   % mean marker, leveled at ~80% ymax
         'MarkerEdgeColor',clr(i,:)*0.6, 'LineWidth',0.8);
end
text(axB, 0.35, 4.65, '$\langle s/l \rangle =1/(M+2)$', 'Interpreter','latex', ...
     'FontSize', fs-1, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
hold(axB,'off');
xlim(axB,[0 1]); ylim(axB,[0 5.3]);
ylabel(axB,'norm. complex density','FontSize',fs);
styleAxis(axB, fs);

%% ================= Panel (b): age CDF =================
axC = axes('Position', posC); hold(axC,'on');
xline(axC, 1, '|', 'Color',[0.3 0.3 0.3], 'LineWidth',0.2);
plot(axC, tau, 1-exp(-tau), '--', 'Color', clr_ref, 'LineWidth', 1.8);   % CDF of the exponential reference: 1 - e^{-tau}
for i = 1:numel(Mlist)
    M   = Mlist(i);
    cdf = gammainc(tau, M+2);                 % regularized lower incomplete gamma
    plot(axC, tau, cdf, 'Color', clr(i,:), 'LineWidth', 1.8);
    tmed = gammaincinv(0.5, M+2);             % median age: CDF = 0.5
    plot(axC, tmed, 0.5, 'o', 'MarkerSize',5, 'MarkerFaceColor',clr(i,:), ...   % median marker, on the curve
         'MarkerEdgeColor',clr(i,:)*0.6, 'LineWidth',0.8);
end
text(axC, 6.1, 0.505, '$\mathrm{med}(\tau)$', 'Interpreter','latex', ...   % inline, right of the median markers (y=0.5)
     'FontSize', fs-1, 'HorizontalAlignment','left', 'VerticalAlignment','middle');
hold(axC,'off');
xlim(axC,[0 12]); ylim(axC,[0 1.02]);
ylabel(axC,'branching CDF','FontSize',fs);
styleAxis(axC, fs);

%% ================= Panel (d): plus-end distance PDF =================
axD = axes('Position', posD); hold(axD,'on');
xline(axD, 1, '|', 'Color',[0.3 0.3 0.3], 'LineWidth',0.2);
for i = 1:numel(Mlist)
    M  = Mlist(i);
    fd = dd.^M.*exp(-dd)/gamma(M+1);
    plot(axD, dd, fd, 'Color', clr(i,:), 'LineWidth', 1.8);
    mu  = M+1;                                % mean distance from plus end
    plot(axD, mu, 0.82, 'o', 'MarkerSize',5, 'MarkerFaceColor',clr(i,:), ...   % mean marker, leveled at ~80% ymax
         'MarkerEdgeColor',clr(i,:)*0.6, 'LineWidth',0.8);
end
text(axD, 3, 0.91, '$\langle d \rangle= (M+1)V_\mathrm{g}/\lambda$', 'Interpreter','latex', ...
     'FontSize', fs-1, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
hold(axD,'off');
xlim(axD,[0 12]); ylim(axD,[0 1.02]);
ylabel(axD,'branching probability','FontSize',fs);
styleAxis(axD, fs);
set(axD,'XDir','reverse');   % plus end on the right, minus-end direction to the left

%% ---- bold panel labels (a)-(d) ----
lbl = {'(a)','(c)','(b)','(d)'};   % column-major: top-right ->(c), bottom-left ->(b)
pos = {posA,posB,posC,posD};
for j = 1:4
    annotation(fig,'textbox',[pos{j}(1)-0.075, pos{j}(2)+pos{j}(4)-0.010, 0.05, 0.04], ...
        'String',lbl{j},'FontSize',fs+2,'FontWeight','bold', ...
        'EdgeColor','none','HorizontalAlignment','left','VerticalAlignment','top', ...
        'Margin',0,'FitBoxToText','off');
end

%% ---- x-labels (placed last: widths measured on the fully rendered figure) ----
drawnow;
% warm-up: the first text-Extent query in a figure can be unreliable; prime it so
% panel (a)'s measurement is as trustworthy as the rest.
axw = axes('Parent',fig,'Position',[0 0 1 1],'Units','normalized', ...
           'Visible','off','XLim',[0 1],'YLim',[0 1]);
tw  = text(axw, 0, 0, 'normalized arclength', 'Interpreter','none', 'FontSize',fs);
drawnow;  tmp = tw.Extent;  %#ok<NASGU>
delete(tw); delete(axw);
xlabSplit(fig, posA, 'parent age',             '$\tau\ (1/\lambda)$',      fs);
xlabSplit(fig, posB, 'normalized arclength',   '$s/L$',                 fs);
xlabSplit(fig, posC, 'parent age',             '$\tau\ (1/\lambda)$',      fs);
xlabSplit(fig, posD, 'distance from plus end', '$d\ (V_\mathrm{g}/\lambda)$', fs);

%% ---- export ----
drawnow;
out = fullfile(outDir, 'Figure8.pdf');
exportgraphics(fig, out, 'ContentType','vector');
fprintf('Saved %s\n', out);

%% ---- local helpers ----
function styleAxis(ax, fs)
    set(ax,'FontSize',fs,'Box','on','LineWidth',1,'TickDir','out');
end

function xlabSplit(fig, pos, txt, mathstr, fs)
    % Two-part xlabel (a la Figure6.m): descriptive txt in the default (non-LaTeX)
    % font, math/units in LaTeX. Rendered widths are measured so the text+math pair
    % is centered on the panel; the LaTeX part is nudged down to match the baseline.
    gap    = 0.008;                % seam gap between text and math (figure fraction)
    dyMath = 0.0031;                % LaTeX part sits slightly lower than the text
    yb     = pos(2) - 0.067;       % vertical position (more spacing below the graph)
    cx     = pos(1) + pos(3)/2;    % panel-center x

    % invisible full-figure axes so that data coordinates == figure fraction
    axm = axes('Parent',fig,'Position',[0 0 1 1],'Units','normalized', ...
               'Visible','off','XLim',[0 1],'YLim',[0 1]);
    tT = text(axm, 0, 0, txt,     'Interpreter','none',  'FontSize',fs, ...
              'HorizontalAlignment','left','VerticalAlignment','middle');
    tM = text(axm, 0, 0, mathstr, 'Interpreter','latex', 'FontSize',fs, ...
              'HorizontalAlignment','left','VerticalAlignment','middle');
    drawnow;
    wT = tT.Extent(3);  wM = tM.Extent(3);
    xleft = cx - (wT + gap + wM)/2;              % center the whole label on the panel
    tT.Position = [xleft,            yb,          0];
    tM.Position = [xleft + wT + gap, yb - dyMath, 0];
end
