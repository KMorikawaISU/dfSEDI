\name{df_get_X}
\alias{df_get_X}
\title{Extract the covariate matrix from a dfSEDI-style data frame}
\description{
Extract a numeric design matrix \code{X} from a data frame in the format expected by \pkg{dfSEDI}.
}
\usage{
df_get_X(dat)
}
\arguments{
\item{dat}{A data frame with a column \code{X}. The \code{X} field can be a numeric matrix stored as a single column, a data.frame (possibly with categorical variables), or a numeric vector.}
}
\value{
A numeric matrix used internally by the estimators.
}
\details{
If \code{dat$X} is a data.frame, \code{df_get_X} uses \code{model.matrix(~ ., data = X)} (without an intercept) to expand categorical variables into dummy columns.
}
\examples{
\dontrun{
X <- df_get_X(dat)
str(X)
}
}
