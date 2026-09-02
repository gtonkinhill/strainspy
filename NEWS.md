# strainspy news

## strainspy 0.99.0

Version set to 0.99.0 in preparation for Bioconductor submission; Bioconductor
requires a `0.99.z` version for new packages and increments it to 1.0.0 at the
next release.

### Bug fixes

- `hadjust()` no longer fails when `taxonomy` is given as a file path. The
  resolved taxonomy table is now passed to the internal ordering step, rather
  than the raw argument.
- `glmZiBFit()` and `glmObFit()` now honour `MAP_prior = NULL`, the documented
  unpenalised maximum-likelihood fit, instead of failing on a `NULL` prior.
- `resolve_priors()` now signals an error for an unrecognised `MAP_prior`.
  Previously any unrecognised value silently fell through to an unpenalised
  fit, so a typo, factor or list quietly changed the model that was fitted.
- Failed per-feature fits are now identified by a single shared
  `failed_fit_index()` helper across all four fitting functions, which handles
  both backend failure conventions. The previous per-function checks each
  covered only one convention.

### New validity checks

- `strainspy_fit` gained a validity method requiring every per-feature slot to
  be aligned with `row_data`, so a length mismatch can no longer silently
  attribute one strain's statistics to another.
- `strainspy_priors` gained a validity method checking `method` and the columns
  `glmmTMB` requires in `priors_df`. A malformed prior previously produced a fit
  object containing zero features with no error surfaced.
- `read_sylph()`, `read_metaphlan()` and `read_sourmash()` now stop when sample
  or feature names are duplicated, which could previously happen when
  `clean_names = TRUE` reduced two input paths to the same basename.

### Dependencies

- `readr` is no longer required; metadata is read and written with `utils`.
- `gamlss` and `gamlss.dist` moved to `Suggests`, as they are only needed for
  the non-default `glmZiBFit(method = "gamlss")` backend.
- `graphics`, `grDevices`, `tools` and `utils` are now declared in `Imports`,
  and functions previously registered with `globalVariables()` are declared
  with `@importFrom`.

### Other changes

- The unused `assay` slot has been removed from `strainspy_fit`. It was never
  populated, so `fit@assay` always returned an empty matrix. Per-feature
  metadata is available from `row_data`.
- Progress output is now emitted with `message()` rather than `cat()`, so it can
  be silenced with `suppressMessages()`.
- Example data in `inst/extdata` has been reduced to keep the source package
  within the Bioconductor size limit. The sylph profile and GTDB taxonomy are
  unchanged; the sourmash example is now a 20-sample subset.

## strainspy 0.1.0

- Initial public release.
- Added expanded test coverage across readers, modeling workflows, outputs, plotting, and utilities.
- Added Bioconductor readiness updates: metadata enhancements, vignette scaffolding, and documentation improvements.
