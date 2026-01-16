## Leunbach test equating using R

This package aims to implement the Leunbach test equating procedures in the software [DIGRAM](https://biostat.ku.dk/DIGRAM/), created by Svend Kreiner and Karl Bang Christensen. Both direct and indirect equating is implemented.

The equating procedure is described in detail in the [Additional file 1](https://link.springer.com/article/10.1186/s12874-019-0768-y#additional-information) related to this article:
Adroher, N. D., Kreiner, S., Young, C., Mills, R., & Tennant, A. (2019). Test equating sleep scales: Applying the Leunbach’s model. *BMC Medical Research Methodology, 19*(1), 141. <https://doi.org/10.1186/s12874-019-0768-y>

## Installation

This package is not yet available on [CRAN](https://cran.r-project.org/) and needs to be installed directly from this code repository. You can do this in multiple ways, I suggest either using `pak` or `remotes` (or `devtools`).

If you want to use parallel processing, which speeds up the bootstrap procedure considerably, you should also install the package `mirai`.

First, install the [`pak`](https://pak.r-lib.org/) package, if you don't already have it installed:
```r
install.packages('pak')
```

Then install the package.
```r
pak::pkg_install("pgmj/leunbachR")
```

If you have the `remotes` or `devtools` package installed, you can instead use:
```r
remotes::install_github("pgmj/leunbachR")
```
or
```r
devtools::install_github("pgmj/leunbachR")
```

There is an [intro article](https://pgmj.github.io/leunbachR/articles/intro.html) under the heading Articles at the top of this page, showing the basic functionality.

## Credits 

- This work was funded by the [Swedish Defence University](https://www.fhs.se/en/swedish-defence-university.html).
- Svend Kreiner kindly shared the source code for the DIGRAM test equating procedure and provided helpful information and guidance.
- Claude Opus 4.5 produced most of the code, based on the DIGRAM pascal code.

## Author

[Magnus Johansson](https://ki.se/en/people/magnus-johansson-3) is a licensed psychologist with a PhD in behavior analysis. He works as a research specialist at [Karolinska Institutet](https://ki.se/en/cns/research/centre-for-psychiatry-research), Department of Clinical Neuroscience, Center for Psychiatry Research.

- ORCID: [0000-0003-1669-592X](https://orcid.org/0000-0003-1669-592X)
- Bluesky: [@pgmj.bsky.social](https://bsky.app/profile/pgmj.bsky.social) 
- Mastodon: [@pgmj@scicomm.xyz](https://scicomm.xyz/@pgmj)