\name{ExtDist-package}
\alias{ExtDist-package}
\alias{ExtDist}
\docType{package}
\title{
Extended Probability Distribution Functions
}
\description{
The package provides a consistent, unified and extensible
    framework for parameter estimation of probability distributions; it extends
    parameter estimation procedures to allow for weighted samples; moreover, it
    extends the gallery of available distributions.
}
\details{
\tabular{ll}{
Package: \tab ExtDist\cr
Version: \tab 0.6\cr
License: \tab GPL (>= 2)\cr
Imports: \tab 
numDeriv,
optimx,
VGAM,
SuppDists,
truncdist\cr
Suggests: \tab 
knitr,
ggplot2,
xtable,
PerformanceAnalytics\cr
Roxygen: \tab list(wrap = FALSE)\cr
VignetteBuilder: \tab knitr\cr
Built: \tab R 3.1.1; ; 2015-03-31 03:52:27 UTC; unix\cr
}

Index:
\preformatted{
Beta                    The standard Beta Distribution.
Beta_ab                 The four-Parameter beta Distribution.
Burr                    The Burr's Distribution.
DistSelCriteriaValues   Distribution Selection Criteria Values.
ExtDist-package         Extended Probability Distribution Functions
Gamma                   The Gamma Distribution.
Gumbel                  The Gumbel Distribution.
JohnsonSB               The Johnson SB Distribution.
JohnsonSU               The Johnson SU Distribution.
Laplace                 The Laplace Distribution.
Logistic                The Logistic Distribution.
Normal                  The Normal Distribution.
Normal_sym_trunc_ab     The symmetric truncated normal distribution.
Normal_trunc_ab         The truncated normal distribution.
SRTB_ab                 The Symmetric-Reflected Truncated Beta (SRTB)
                        Distribution.
SSRTB                   The standard Symmetric-Reflected Truncated Beta
                        (SRTB) Distribution.
Triangular              The Triangular Distribution.
Uniform                 The Uniform Distribution.
Weibull                 The Weibull Distribution.
bestDist                Best distribution for (weighted) sample.
compareDist             Compare sample and fitted distributions
eDist                   S3 methods for manipulating eDist objects.
eval.estimation         Parameter Estimation Evaluation.
wmle                    Weighted Maximum Likelihood Estimation.
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{Distributions-Beta} \tab Distributions-Beta (source)\cr
\code{Distributions-Normal} \tab Distributions-Normal (source)\cr
\code{Distributions} \tab Distributions-Index (source)\cr
\code{ParaEst-and-DistSel-by-ExtDist} \tab Parameter-Estimation-and-Distribution-Selection-by-ExtDist (source)\cr
}
}
\author{
Haizhen Wu <h.wu2@massey.ac.nz>,
    A. Jonathan R. Godfrey <A.J.Godfrey@massey.ac.nz>,
    Kondaswamy Govindaraju <k.govindaraju@massey.ac.nz>,
    Sarah Pirikahu <s.pirikahu@massey.ac.nz>

Maintainer: Oleksii Nikolaienko <oleksii.nikolaienko@gmail.com>
}
