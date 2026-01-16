\name{df_estimate_NP}
\alias{df_estimate_NP}
\title{Non-probability-sample-only estimator}
\description{
Compute a non-probability-sample-only estimator of \eqn{\theta = E(Y)}.
This function requires a user-supplied basis function \code{base_fun} and is mainly
intended for simulation and benchmarking.
}
\usage{
df_estimate_NP(dat, base_fun, phi_start = NULL, max_iter = 20, N = NULL, hajek = FALSE, pi_x_cols = NULL)
}
\arguments{
\item{dat}{A data frame containing at least \code{X}, \code{y}, and \code{d_np}.}
\item{base_fun}{A function that maps \eqn{L = (1, X_{\pi}, Y)} to a basis matrix used in the estimating equations.}
\item{phi_start}{Optional initial value for \eqn{\phi}.}
\item{max_iter}{Maximum number of optimization attempts used internally.}
\item{N}{Optional population size. Required when \code{hajek=FALSE} and the input does not contain full-frame rows.}
\item{hajek}{Logical. If \code{TRUE}, computes a Hajek (ratio) version that does not require \code{N}.}
\item{pi_x_cols}{Optional restriction on which columns of \code{X} are used in the response-mechanism model.}
}
\value{
A list with components \code{phi}, \code{theta}, \code{var}, \code{se}, and \code{ci}.
}
\seealso{\code{\link{df_estimate_NP_P}}, \code{\link{Eff}}}
\examples{
\dontrun{
# Example base function that uses an intercept and all columns of L
base_fun <- function(L) cbind(1, L)
fit <- df_estimate_NP(dat, base_fun = base_fun, hajek = TRUE)
fit$theta
}
}
