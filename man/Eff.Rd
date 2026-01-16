\name{Eff}
\alias{Eff}
\title{Semiparametric efficient estimator for dual-frame data integration}
\description{
Compute the semiparametric efficient estimator for the population mean
\eqn{\theta = E(Y)} under a dual-frame data integration setting.
}
\usage{
Eff(dat,
    K = 2,
    phi_start = NULL,
    max_restart = 10,
    type = 2,
    dml_type = 1,
    progress = interactive(),
    x_info = TRUE,
    prob_only = FALSE,
    N = NULL,
    nonpara_method = "KRR",
    pi_x_cols = NULL)
}
\arguments{
\item{dat}{
A data frame containing at least the columns \code{X}, \code{y}, \code{d_np}, \code{d_p}, and \code{pi_p}.
See the package README for the expected data format.
}
\item{K}{Number of folds for cross-fitting.}
\item{phi_start}{Optional initial value for the response-mechanism parameter \eqn{\phi}.}
\item{max_restart}{Maximum number of restarts used internally for solving estimating equations.}
\item{type}{Backward-compatible alias for \code{dml_type}. If not \code{NULL}, it overrides \code{dml_type}.}
\item{dml_type}{DML variant: \code{1}/\code{"DML1"} or \code{2}/\code{"DML2"}.}
\item{progress}{Logical; if \code{TRUE}, prints fold-level progress.}
\item{x_info}{
Logical. If \code{TRUE}, uses the full formulation (requires/uses the (d_np, d_p)=(0,0) pattern if present).
If \code{FALSE}, uses the union-sample formulation that sets augmentation terms to zero.
}
\item{prob_only}{Logical. Only relevant when \code{x_info=TRUE}; selects the probability-sample-only representation for augmentation-term estimation.}
\item{N}{
Optional population size. Required to target the full-frame mean when \code{x_info=FALSE} and the input data contain only the union sample.
}
\item{nonpara_method}{
Method used for nonparametric nuisance regressions. Currently supports kernel ridge regression (KRR) and several alternatives.
}
\item{pi_x_cols}{
Optional restriction on which columns of \code{X} are used in the response-mechanism model \eqn{\pi_{np}(X,Y)}.
Can be an integer, logical, or character vector; see README for details.
}
}
\value{
A list with components:
\describe{
  \item{theta}{Point estimate of \eqn{\theta = E(Y)}.}
  \item{var}{Estimated variance of \code{theta}.}
  \item{se}{Estimated standard error of \code{theta}.}
  \item{ci}{Confidence interval for \code{theta} (numeric vector of length 2).}
  \item{phi}{Estimated response-mechanism parameter \eqn{\phi}.}
  \item{phi_var}{Estimated variance-covariance matrix of \code{phi}.}
  \item{phi_se}{Standard errors for \code{phi}.}
  \item{phi_ci}{Confidence intervals for \code{phi}.}
  \item{info}{A list of run metadata (fold count, options used, convergence diagnostics, etc.).}
}
}
\references{
Morikawa and Kim (2026). Semiparametric efficient data integration for dual-frame sampling.
}
\examples{
\dontrun{
# dat must contain (at least) X, y, d_np, d_p, pi_p
fit <- Eff(dat, K = 2, dml_type = "DML2", x_info = FALSE, N = 100000)
fit$theta
fit$se
}
}
