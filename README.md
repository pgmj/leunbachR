# leunbachR: Leunbach test equating using R

<!-- badges: start -->
[![CRAN Version](https://www.r-pkg.org/badges/version/leunbachR)](https://cran.r-project.org/package=leunbachR)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Downloads](https://cranlogs.r-pkg.org/badges/leunbachR?color=brightgreen)](https://CRAN.R-project.org/package=leunbachR)
![Downloads Status](https://cranlogs.r-pkg.org/badges/grand-total/leunbachR)
<a href="https://buymeacoffee.com/pgmj" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/default-orange.png" alt="Buy Me A Coffee" height="41" width="174"></a>
<!-- badges: end -->

This package aims to implement the Leunbach test equating procedures in the software [DIGRAM](https://biostat.ku.dk/DIGRAM/), created by Svend Kreiner and Karl Bang Christensen. Both direct and indirect equating is implemented.

The equating procedure is described in detail in the [Additional file 1](https://link.springer.com/article/10.1186/s12874-019-0768-y#additional-information) related to this article:
Adroher, N. D., Kreiner, S., Young, C., Mills, R., & Tennant, A. (2019). Test equating sleep scales: Applying the Leunbach’s model. *BMC Medical Research Methodology, 19*(1), 141. <https://doi.org/10.1186/s12874-019-0768-y>

## Installation

You can install the package from [CRAN](https://cran.r-project.org/):
```r
install.packages('leunbachR')
```

If you want to use parallel processing, which speeds up the bootstrap procedure considerably, you should also install the package `mirai`.

There is an [intro article](https://pgmj.github.io/leunbachR/articles/intro.html) under the heading Articles at the top of this page, showing the basic functionality.

## Credits 

- This work was funded by the [Swedish Defence University](https://www.fhs.se/en/swedish-defence-university.html) and the [Swedish Defence Conscription and Assessment Agency](https://pliktverket.se/om-myndigheten/in-english).
- Svend Kreiner kindly shared the source code for the DIGRAM test equating procedure and provided helpful information and guidance.
- Claude Opus 4.5 produced most of the code, based on the DIGRAM pascal code.
- Jeanette Melin and Henrik Nordström at the Swedish Defence University helped with testing and provided useful feedback.

## Author

[Magnus Johansson](https://ki.se/en/people/magnus-johansson-3) is a licensed psychologist with a PhD in behavior analysis. He works as a research specialist at [Karolinska Institutet](https://ki.se/en/cns/research/centre-for-psychiatry-research), Department of Clinical Neuroscience, Center for Psychiatry Research.

- ORCID: [0000-0003-1669-592X](https://orcid.org/0000-0003-1669-592X)
- Bluesky: [@pgmj.bsky.social](https://bsky.app/profile/pgmj.bsky.social) 