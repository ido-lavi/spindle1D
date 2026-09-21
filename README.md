# spindle1D

MATLAB code accompanying **Proliferating against the current: how multistep branching nucleation can sort microtubule polarity in the spindle**. The repository includes a continuum
solver, stochastic filament visualizations, and scripts and data for reproducing
the manuscript’s numerical and analytical figures.

The model is a closed hierarchy of coupled one-dimensional transport equations.
Microtubule plus-end density sources lattice density, which feeds a sequence of
branching-complex assembly stages. The final stage couples to a Ran-activated
nucleator field to produce new plus ends, closing the feedback loop. Growth and
sliding transport the populations at different velocities, while microtubule
turnover removes lattice sites and their associated complexes. A separate
reaction–diffusion–advection equation governs the nucleator field.

## Requirements

MATLAB **R2022a or newer** (the phase-diagram script uses `clim`; the figures use
`exportgraphics`).

The **Statistics and Machine Learning Toolbox** is required by the stochastic-filament
scripts — `solver/filamentSteadyState.m`, `solver/filamentMovie.m`, and
`figures/FigureDynamics.m` (they use `randsample` / `exprnd`). The continuum solver
`PDEsolver.m` and the field-only / analytic figures do not need it.

## Layout

```
spindle1D/
├── solver/          continuum solver + stochastic-filament steady states / movies
│   ├── PDEsolver.m             pseudo-spectral PDE solver (set params, run)
│   ├── fieldSteadyState.m      save a steady-state snapshot: continuum fields only
│   ├── filamentSteadyState.m   reconstruct psi(x,l), sample filaments, save snapshot
│   └── filamentMovie.m         stochastic-filament movie from the saved nuc_t(x,t)
├── figures/         one script per manuscript figure (self-contained)
│   ├── FigureDynamics.m        Figures 2 & 5
│   ├── Figure3.m               Figure 3
│   ├── Figure6.m               Figure 6
│   ├── Figure7.m               Figure 7
│   ├── Figure8.m               Figure 8
│   └── Figure9.m               Figure 9
├── data/            steady-state snapshots and time-series inputs
└── output/          generated PDFs / movies (git-ignored)
```

Every script resolves `data/` and `output/` relative to its own location, so it
can be run from any working directory.

## Pipeline

All quantities are dimensionless: time in `1/λ`, space in `ℓ = V_g/λ`
(`V_g = λ = 1`). The two control parameters are `eps = V∥/V_g` and the number of
sequential assembly steps `M`. The binding affinity `alpha` set in `PDEsolver.m` is a
reference value that the solver rescales by `1/(n_b·a_0·a_m^M)`, so that varying
`n_b` changes the overall filament-count without affecting the solution shape.

`PDEsolver.m` writes the time-series directly. To make a steady-state snapshot,
run **one** of the two savers afterwards in the same workspace — they are
alternatives, both writing `sim_M<M>_eps<eps>.mat`:

```
PDEsolver.m ──► data/fields_*.mat , data/nuc_t_*.mat   (time-series → FigureDynamics, filamentMovie)
      │
      │   then, in the same workspace, run ONE saver:
      ├──► fieldSteadyState.m     ──► data/sim_*.mat    (continuum fields only)
      └──► filamentSteadyState.m  ──► data/sim_*.mat    (fields + sampled filaments; psi reconstructed internally, not saved)
```

- `fieldSteadyState.m` writes a light snapshot — enough for the Figure 7
  phase-diagram grid, which reads only `phi`/`chain`.
- `filamentSteadyState.m` is the superset: it reconstructs `psi(x,l)`, samples the
  discrete filaments (`rod_*`), and saves them alongside the fields — needed by the
  field+filament panels of Figures 3, 6, 9.
- Snapshot filenames are `sim_M<M>_eps<eps>.mat`, with `_aself<a_self>` appended
  only when the cooperative self-recruitment rate `a_self > 0`.

## Reproduce figures

The repository ships all steady-state snapshots (`sim_*.mat`, ~8.5 MB —
noncooperative plus the cooperative Figure 9 sweep) and the Figure 2 time-series
(`fields_M0_eps0.200.mat` and `nuc_t_M0_eps0.200.mat`, ~82 MB together) in `data/`.
Only the larger time-series (Figure 5 and the movies, up to ~450 MB each) are
regenerated with the solver; those are git-ignored.

**Figures 3, 6, 7, 9.** Run any of them:

```matlab
run figures/Figure3.m     % Fig 3
run figures/Figure6.m     % Fig 6
run figures/Figure7.m     % Fig 7
run figures/Figure9.m     % Fig 9
```

**Figure 8** is fully analytic — it needs no data files:

```matlab
run figures/Figure8.m     % Fig 8
```

**Figure 2 (M=0 relaxation):**

```matlab
run figures/FigureDynamics.m     % matBase = 'M0_eps0.200' (default) -> Fig 2
```

**Figure 5 (M=5).** Generate its time-series with the solver, then plot:

```matlab
cd solver
% edit PDEsolver.m: set M = 5, eps = 0.2, then run
PDEsolver                       % writes data/fields_M5_eps0.200.mat and nuc_t_...
cd ../figures
% set matBase = 'M5_eps0.200' in FigureDynamics.m, then run
FigureDynamics
```

For a movie, run `solver/filamentMovie.m` after setting `matName` to the
corresponding `nuc_t_*.mat`.

**Figure 9 (cooperative model)** runs directly from the shipped cooperative
snapshots (listed above). To **regenerate** them, sweep `a_self` over
`{0.2, 0.4, 0.6, 0.7, 0.8, 0.9}` (with `M=1`, `eps=0.2`); the `a_self = 0` point
reuses the noncooperative `sim_M1_eps0p20.mat`:

```matlab
cd solver
% edit PDEsolver.m: M = 1, eps = 0.2, a_self = <value>, then
PDEsolver; filamentSteadyState;      % writes data/sim_M1_eps0p20_aself<...>.mat
cd ../figures
Figure9
```

(Figures 1 and 4 in the manuscript are schematics, not produced here.)

## Saved variables

Each snapshot / time-series `.mat` stores dimensionless fields on the periodic grid
`x` (1×`N`). Orientation index **1 denotes right-pointing MTs** (plus ends toward
increasing `x`); index **2 denotes left-pointing MTs** (plus ends toward decreasing
`x`). Orientation is the second dimension of `phi` and `eta`, and the third
dimension of `chain`.

| Variable | Size | Meaning |
|---|---|---|
| `x` | 1×N | periodic spatial grid (`ℓ = V_g/λ`) |
| `phi` | N×2 | microtubule lattice density |
| `eta` | N×2 | plus-end density |
| `chain` | M×N×2 | the `M` sequential branching-complex stages (M>0) |
| `c` | N×1 | Ran-activated nucleator |
| `s` | N×1 | RanGTP signal |
| `pol` | N×1 | polarity `(phi_+ − phi_−)/(phi_+ + phi_−)` |
| `rod_x_p`,`rod_y_p`,`rod_x_m`,`rod_y_m` | 3×Nf | sampled filament endpoints (plus / minus populations) |
| `phi_t`,`eta_t` | N×2×T | field time-series (`fields_*.mat`) |
| `c_t` | N×T | nucleator time-series |
| `chain_t` | M×N×2×T | chain-stage time-series (M>0) |
| `nuc_t` | N×2×T | nucleation-rate series (`nuc_t_*.mat`) |
| `dt_save` | scalar | time between saved frames |
| `params` | struct | all run parameters (see below) |

`psi(x,l)` (the full length-resolved number density) is reconstructed internally by
`filamentSteadyState.m` / `FigureDynamics.m` but is **not** saved. Newly generated
snapshots carry a `params` struct (`M, eps, a_self, a_m, a0, nb, ls, D_c, lambda_c,
s_0, a_c, alpha, L, N, dt, Tfinal`) recording the generating run; snapshots shipped
before this field was added do not have it.

## Reproducibility notes

- **"Steady state"** here means the state at the final integration time
  (`Tfinal = 200`). The solver does not test for convergence, so after changing
  parameters, check that the fields have settled (e.g. inspect the diagnostic GIF)
  before trusting a snapshot.
- **Continuum fields are deterministic and reproducible.** The **sampled filament
  arrangements** (`rod_*`, and the movies) are stochastic illustrations that vary from
  run to run — `solver/filamentSteadyState.m` and `solver/filamentMovie.m` are not
  seeded. `figures/FigureDynamics.m` seeds its sampler (`rng(42)`), so its filament
  panels are reproducible.

## Citation

If you use this code, please cite:

> I. Lavi, S. Fürthauer, D. Needleman, and M. Shelley,
> "*Proliferating against the current: how multistep branching nucleation can sort microtubule polarity in the spindle*."
