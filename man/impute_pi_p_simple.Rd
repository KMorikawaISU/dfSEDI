\name{impute_pi_p_simple}
\alias{impute_pi_p_simple}
\title{Impute missing probability-sample inclusion probabilities}
\description{
Impute missing values of \code{pi_p} for units in the union sample.
This helper is mainly used internally by estimators that require \code{pi_p} on the union sample.
}
\usage{
impute_pi_p_simple(dat, sigma = NULL, nonpara_method = "KRR")
}
\arguments{
\item{dat}{A data frame containing at least \code{X}, \code{y}, and \code{pi_p}.}
\item{sigma}{Optional bandwidth/scale parameter used by the underlying kernel ridge regression. If \code{NULL}, it is estimated from the data.}
\item{nonpara_method}{Reserved for compatibility. \code{pi_p} imputation is currently performed via kernel ridge regression.}
}
\value{
The input data frame with \code{pi_p} replaced by imputed values where missing.
}
\examples{
\dontrun{
# Impute missing pi_p (if any)
dat2 <- impute_pi_p_simple(dat)
}
}
