# Changelog

## leunbachR 0.1.0

- First CRAN submission.
- Made `mirai` an optional Suggest; parallel processing is opt-in via
  `parallel = TRUE` and an explicit `n_cores` value. The package now has
  no hard dependencies beyond base R.
- Moved example datasets to `inst/extdata/` and updated the vignette to
  load them via
  [`system.file()`](https://rdrr.io/r/base/system.file.html).
- Replaced non-ASCII characters in documentation and console output.
- Added unit tests under `tests/testthat/` and runnable examples for all
  exported functions.

## leunbachR 0.0.3.1 (2026-01-27)

- Removed theta values output from indirect equating print functions.

## leunbachR 0.0.3 (2026-01-23)

- Fix for LRT df calculations.
- Fix for gamma test SE, p-value, reverse one-sided test.
- Fix for indirect equating boundary case.
- Fix for text labels of tests in indirect equating output.

## leunbachR 0.0.2 (2026-01-20)

- Theta scores added to output.
- Vignette describes how to export bootstrap equating table to CSV file.

## leunbachR 0.0.1 (2026-01-15)

- Initial release
