\name{impute_pi_p_crossfit}
\alias{impute_pi_p_crossfit}
\title{Cross-fitted imputation of probability-sample inclusion probabilities}
\description{
Impute missing values of \code{pi_p} using K-fold cross-fitting.
This is useful when downstream estimators also use cross-fitting and you want to avoid reusing the same data for imputation and estimation.
}
\usage{
impute_pi_p_crossfit(dat, folds, sigma = NULL, progress = FALSE, nonpara_method = "KRR")
}
\arguments{
\item{dat}{A data frame containing at least \code{X}, \code{y}, and \code{pi_p}.}
\item{folds}{A list of integer index vectors. Each element indicates the held-out indices for a fold.}
\item{sigma}{Optional bandwidth/scale parameter used by the underlying kernel ridge regression. If \code{NULL}, it is estimated from the training data for each fold.}
\item{progress}{Logical; if \code{TRUE}, prints fold-level progress.}
\item{nonpara_method}{Reserved for compatibility. \code{pi_p} imputation is currently performed via kernel ridge regression.}
}
\value{
The input data frame with \code{pi_p} replaced by cross-fitted imputed values where missing.
}
\examples{
\dontrun{
# Example fold definition (2-fold):
set.seed(1)
K <- 2
fold_id <- sample(rep(1:K, length.out = nrow(dat)))
folds <- split(seq_len(nrow(dat)), fold_id)

dat2 <- impute_pi_p_crossfit(dat, folds = folds, progress = TRUE)
}
}
