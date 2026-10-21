# Changelog

## [Unreleased] (1.1)

### Fixed
- `GMGTS()` called the undefined function `ODEIQM` and passed positional
  arguments (`t`, `observed`, `init`) on to the Name-Value parser; it now
  works with both the struct and the positional calling syntax.
- `Generator` called `ODEIQM` when given a model file name, and rejected
  `error_const=0` / `error_std=0` although `0` is the default.
- `Estimator` crashed when numeric data was passed without positional
  arguments; `out.dfitted` contained fitted states instead of gradients;
  `loglik()`/`loglik_lognormal()` used GMGTS cell estimates when `"GTS"` was
  selected; derivative plots were titled with the wrong state; contours of
  log-normal random effect distributions used an incorrect density; LB/UB
  given as column vectors produced wrong constraints.
- `FirstStage` (partially observed systems): the final uncertainty update
  returned after the first cell, so the remaining cells kept the
  uncertainty estimates of the previous iteration. **May change results.**
- `FirstStage` (fully observed systems): crashed when the first stage
  stopped after a single iteration.
- `Optimization.QPGLS` always overwrote its `prior` argument, so the `Prior`
  option was ignored in the gradient-matching step (no script in this
  repository uses `Prior`).
- `System`: changing `system.fixed.values` after construction only affected
  ODE integration, not the symbolic right-hand side used for gradient
  matching. Fixed parameter values are now changed with
  `system.set_fixed_values(names, values)`, which rebuilds all derived
  functions; `fixed` is read-only. This affected
  `inference_experimental.m`, where FP-specific dilution rates were set
  this way. **May change results.**
- `System` now accepts model file names on the MATLAB path and model files
  in other folders (the MEX name was derived from the full path).
- `tryinv` silently returned nothing for unexpected errors; warnings were
  globally re-enabled after the second stage.
- Experiment scripts: missing `hedist` function, missing file extension of
  the tdTomato data, results folder created only after the first save,
  repressilator script loaded an untracked `.mat` file, inconsistent
  (unnormalized) reference distance in the maturation accuracy study.
- Example scripts: `Method=`/`False` typos, comments not matching settings.
- Removed `BSpline.spline_ddiff` (called a nonexistent method) and leftover
  debug output.

### Changed
- Repository layout: `src/`, `models/`, `examples/`, `experiments/`,
  `data/`, `docs/`; new `setup_GMGTS.m` for path setup and dependency
  checks. Experiment scripts write to `results/` and can be run from any
  folder. `maturation_accuracy2.mat` is now `maturation_accuracy.mat`.
- Documentation: rewritten README (requirements, installation,
  reproduction guide, corrected defaults), LICENSE file, data README,
  valid CITATION.cff with the article as preferred citation.

## [1.0] - 2025-02-17

Initial release accompanying the manuscript (Zenodo DOI
10.5281/zenodo.14884458).
