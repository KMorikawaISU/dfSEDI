\docType{package}
\name{dfSEDI-package}
\alias{dfSEDI}
\alias{dfSEDI-package}
\title{dfSEDI: Semiparametric Efficient Data Integration for Dual-Frame Sampling}
\description{
The \pkg{dfSEDI} package implements semiparametric efficient data integration
estimators for dual-frame sampling.
}
\details{
The main user-facing estimators are:
\itemize{
  \item \code{\link{Eff}}: semiparametric efficient estimator (supports DML1/DML2 and K-fold cross-fitting)
  \item \code{\link{Eff_S}}: sub-efficient estimator
  \item \code{\link{Eff_P}}: parametric efficient estimator under a working model
}

For simulation and benchmarking, the package also includes:
\itemize{
  \item \code{\link{df_estimate_P}}: probability-sample-only estimator
  \item \code{\link{df_estimate_NP}}: non-probability-sample-only estimator (requires \code{base_fun})
  \item \code{\link{df_estimate_NP_P}}: union-sample estimator combining NP and P information (requires \code{base_fun})
}

All estimators target a population mean \eqn{\theta = E(Y)} and return a point estimate
with sandwich-type standard errors and confidence intervals.
}
\seealso{
\itemize{
  \item \code{\link{Eff}}, \code{\link{Eff_S}}, \code{\link{Eff_P}}
  \item \code{\link{df_estimate_P}}, \code{\link{df_estimate_NP}}, \code{\link{df_estimate_NP_P}}
}
}
\author{
Kosuke Morikawa and Jae Kwang Kim.
}
\keyword{package}
