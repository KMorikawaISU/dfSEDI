\name{Eff_S}
\alias{Eff_S}
\title{Sub-efficient estimator for dual-frame data integration}
\description{
Compute the sub-efficient estimator (\code{Eff_S}) for the population mean
\eqn{\theta = E(Y)}.
}
\usage{
Eff_S(dat,
      K = 2,
      logit = NULL,
      nonpara_method = "KRR",
      MAR = TRUE,
      progress = interactive(),
      x_info = TRUE)
}
\arguments{
\item{dat}{A data frame containing the required dual-frame columns (see \code{\link{Eff}}).}
\item{K}{Number of folds for cross-fitting.}
\item{logit}{Logical or \code{NULL}. If \code{TRUE} and \code{y} is binary, uses a logistic nuisance regression when applicable. If \code{NULL}, the behavior is inferred from \code{y}.}
\item{nonpara_method}{Nonparametric regression method used for nuisance estimation.}
\item{MAR}{Logical. If \code{TRUE}, uses the MAR-style augmentation. If \code{FALSE}, uses the Eq.(11) ratio/shortcut formulation.}
\item{progress}{Logical; if \code{TRUE}, prints fold-level progress.}
\item{x_info}{Logical. See \code{\link{Eff}} for the meaning of \code{x_info}.}
}
\value{
A list with components:
\describe{
  \item{theta}{Point estimate of \eqn{\theta = E(Y)}.}
  \item{var}{Estimated variance of \code{theta}.}
  \item{se}{Estimated standard error of \code{theta}.}
  \item{ci}{Confidence interval for \code{theta}.}
  \item{info}{A list of run metadata (method used, fold count, etc.).}
}
}
\examples{
\dontrun{
fit <- Eff_S(dat, K = 2, nonpara_method = "KRR")
fit$theta
}
}
