# Birth–Death Gaussian Graphical Models

Code for *Birth–Death Gaussian Graphical Models for Structurally Evolving
Systems* (Bortolato, Hansen, Rossell). The model is a Gaussian graphical model
whose set of variables changes over time. The precision matrix on the active
set is low rank plus diagonal, and a spike-and-slab prior on rowwise loading
changes, informed by birth and death counts, drives its evolution. Estimation
uses a generalized EM algorithm.

## Layout

```text
.
├── Project.toml                         Julia dependencies
├── src/
│   ├── BirthDeathGraphicalModels.jl     estimator, simulator, diagnostics
│   └── ControlledPrecisionExperiments.jl  scenarios defined directly by precision matrices
├── experiments/
│   ├── birthdeath/run_experiment2.jl    Experiment II (paper, Section 7.4)
│   └── constant_precision/run_constant_precision_p5.jl
│                                        check that no changes are invented when the truth is constant
├── competitors/
│   └── FactorPrecision.R                precision-factor (PF) competitor sampler
├── results/
│   └── experiment2/                     CSV output of run_experiment2.jl
└── legacy/                              earlier MCMC samplers, kept for reference
```

The scripts for Experiment I (piecewise-constant precision), Experiment III
(endogenous exit), the competitor fits and the prior plots are not in the
repository yet.

## Setup

Install Julia 1.10 or later, then from the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`src/BirthDeathGraphicalModels.jl` loads Plots for its plotting helpers. On a
machine without a display, set `GKSwstype=100` before running.

## Running the experiments

All commands are run from the repository root.

**Experiment II.** One CSV row per simulated trajectory:

```bash
julia --project=. experiments/birthdeath/run_experiment2.jl \
    p=10 seeds=1:20 tau_active=0.75 minimum_active=6 \
    out=results/experiment2/p10_tau075_min6.csv
```

| Argument | Meaning | Default |
|---|---|---|
| `p` | number of variables (10 or 50) | 10 |
| `seeds` | `a:b` range or comma-separated list | `1:5` |
| `tau_active` | fitted active innovation scale τ_a | 0.75 |
| `minimum_active` | lower bound on the number of active variables p_t | k + 1 = 6 |
| `out` | output CSV path | none, prints only |

Every other setting is the paper's, and is fixed in the script.

**Constant precision, p = 5:**

```bash
julia --project=. experiments/constant_precision/run_constant_precision_p5.jl
```

It writes summaries and per-time metrics to
`experiments/constant_precision/constant_precision_p5_output/`.

## Evaluation metrics

`active_precision_errors` reports, at each time point and on the active set:

- `rv`: RV coefficient of the full precision matrices;
- `offdiagonal_rv`: RV on the off-diagonal entries only;
- `partial_correlation_rv`: RV on the off-diagonal partial correlations;
- relative Frobenius error, mean absolute error and RMSE.

The full-matrix RV is dominated by the diagonal. It can stay high when the
off-diagonal structure is poorly recovered, so report the off-diagonal or
partial-correlation RV alongside it.

## Implementation notes

- **Binary-splitting initialization.** `masked_segment_bic` is a deliberate
  pseudo-BIC, not a standard BIC; see its docstring. In practice it accepts
  every split, so the recursion is stopped by `min_segment_length` and
  `changepoint_max_depth`.
- **Minimum active set.** `generate_birth_death_times(...; minimum_active)`
  postpones a death that would leave fewer than `minimum_active` variables.
  Use `k + 1` to respect assumption (S4), p_t > k. The default of 1 only rules
  out empty time points. `EM_time_varying` warns when the data contain time
  points with p_t ≤ k or p_t = 0.
- **Births in Experiment II.** Each variable has a single active spell, so the
  number of births per trajectory is fixed: p minus the number initially
  active. Only the birth times are random.
- **`legacy/`** holds earlier MCMC samplers: `mcmc_sampler.jl`,
  `hmm_sampler_backup.jl` and `toy_example_sampler.R`. They are not used by the
  current estimator.
