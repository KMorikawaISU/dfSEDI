\name{df_estimate_P}
\alias{df_estimate_P}
\title{Probability-sample-only estimator}
\description{
Compute a probability-sample-only estimator of \eqn{\theta = E(Y)}.
This function is primarily intended for simulation and benchmarking.
}
\usage{
df_estimate_P(dat, N = NULL, hajek = FALSE)
}
\arguments{
\item{dat}{A data frame containing at least \code{y}, \code{d_p}, and \code{pi_p}.}
\item{N}{Optional population size. Required when \code{hajek=FALSE} and the input does not contain full-frame rows.}
\item{hajek}{Logical. If \code{TRUE}, computes a Hajek (ratio) version that does not require \code{N}.}
}
\value{
A list with components \code{theta}, \code{var}, \code{se}, and \code{ci}.
}
\examples{
\dontrun{
fit <- df_estimate_P(dat, hajek = TRUE)
fit$theta
}
}
