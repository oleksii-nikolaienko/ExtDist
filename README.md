# ExtDist: Extending the Range of Functions for Probability Distributions

[![](https://www.r-pkg.org/badges/version/ExtDist)](https://cran.r-project.org/package=ExtDist)
[![](https://github.com/oleksii-nikolaienko/ExtDist/workflows/R-CMD-check/badge.svg)](https://github.com/oleksii-nikolaienko/ExtDist/actions)
[![](https://badges.cranchecks.info/summary/ExtDist.svg)](https://cran.r-project.org/web/checks/check_results_ExtDist.html)
[![](https://cranlogs.r-pkg.org/badges/ExtDist?color=blue)](https://cran.r-project.org/package=ExtDist)
[![](https://cranlogs.r-pkg.org/badges/grand-total/ExtDist?color=blue)](https://cran.r-project.org/package=ExtDist)

A consistent, unified and extensible framework for estimation of parameters for probability distributions, including parameter estimation procedures that allow for weighted samples; the current set of distributions included are: the standard beta, the four-parameter beta, Burr, gamma, Gumbel, Johnson SB and SU, Laplace, logistic, normal, symmetric truncated normal, truncated normal, symmetric-reflected truncated beta, standard symmetric-reflected truncated beta, triangular, uniform, and Weibull distributions; decision criteria and selections based on these decision criteria.

Online version of documentation is available at [GitHub Pages](https://oleksii-nikolaienko.github.io/ExtDist/reference/index.html).

### Notable updates

##### 0.7-1 (2023-01-16)

Long-standing [bug](https://stackoverflow.com/questions/45208176/the-weibull-distribution-in-r-extdist) in Weibull and gamma distributions was fixed:

- for gamma distribution, `scale` parameter was passed as a `rate` in `stats::*gamma` calls

- for Weibull distribution, `stats::*gamma` was called instead of `stats::*weibull`

