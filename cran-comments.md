# cran-comments.md

## Submission

Thanks for the review. In response to the comment about writing to the user's home filespace:

The two write.csv() calls in vignettes/leunbachR.Rmd (which generate inst/doc/leunbachR.R) now write to file.path(tempdir(), ...) instead of the working directory.
No other code in the package writes outside tempdir().

## Test environments

* local macOS, R 4.5.3
* Win devel and release
* R-hub: ubuntu-latest, windows-latest, macos-latest

## R CMD check results

0 errors | 0 warnings | 1 notes

Note on CRAN submission and spelling in DESCRIPTION, the latter regarding names.

## Reverse dependencies

This is a new release, so there are no reverse dependencies to check.

## Notes

The package has no hard dependencies beyond base R. The `mirai` package is
listed in Suggests and used only when the user explicitly requests parallel
processing in the bootstrap functions.
