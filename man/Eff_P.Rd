\name{Eff_P}
\alias{Eff_P}
\title{Parametric efficient estimator under a working model}
\description{
Compute the parametric efficient estimator (\code{Eff_P}) for the population mean
\eqn{\theta = E(Y)} under a parametric working model for the augmentation term.
}
\usage{
Eff_P(dat,
      phi_start = NULL,
      eta4_star = 0,
      max_iter = 20,
      logit = NULL,
      progress = interactive(),
      x_info = TRUE,
      pi_x_cols = NULL)
}
\arguments{
\item{dat}{A data frame containing the required dual-frame columns (see \code{\link{Eff}}).}
\item{phi_start}{Optional initial value for the response-mechanism parameter \eqn{\phi}.}
\item{eta4_star}{(Reserved for compatibility) working value for an augmentation term parameter.}
\item{max_iter}{Maximum number of iterations used for solving the estimating equations.}
\item{logit}{Logical or \code{NULL}. If \code{TRUE}, fits a logistic working model for binary \code{y}. If \code{NULL}, the behavior is inferred from \code{y}.}
\item{progress}{Logical; if \code{TRUE}, prints iteration-level progress.}
\item{x_info}{Logical. See \code{\link{Eff}} for the meaning of \code{x_info}.}
\item{pi_x_cols}{Optional restriction on which columns of \code{X} are used in the response-mechanism model.}
}
\value{
A list with components:
\describe{
  \item{theta}{Point estimate of \eqn{\theta = E(Y)}.}
  \item{var}{Estimated variance of \code{theta}.}
  \item{se}{Estimated standard error of \code{theta}.}
  \item{ci}{Confidence interval for \code{theta}.}
  \item{phi}{Estimated response-mechanism parameter \eqn{\phi}.}
  \item{info}{A list of run metadata (options used, convergence diagnostics, etc.).}
}
}
\note{\code{Eff_P} currently returns \code{phi_var}, \code{phi_se}, and \code{phi_ci} as \code{NULL}.}
\examples{
\dontrun{
fit <- Eff_P(dat, logit = NULL)
fit$theta
}
}
