# cran-comments.md

## Submission

This is the first submission of `leunbachR` to CRAN.

`leunbachR` implements the Leunbach test equating method (Adroher et al.
2019, <doi:10.1186/s12874-019-0768-y>), following the DIGRAM software
written by Svend Kreiner. The package supports both direct and indirect
equating with parametric bootstrap standard errors.

## Test environments

* local macOS, R 4.5.3
* `devtools::check_win_devel()`
* R-hub: ubuntu-latest, windows-latest, macos-latest

## R CMD check results

0 errors | 0 warnings | 1 notes

Note on CRAN submission.

## Reverse dependencies

This is a new release, so there are no reverse dependencies to check.

## Notes

The package has no hard dependencies beyond base R. The `mirai` package is
listed in Suggests and used only when the user explicitly requests parallel
processing in the bootstrap functions.
