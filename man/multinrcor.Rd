\name{multinrcor}
\alias{multinrcor}
\title{Multiple Correlations}
\author{
Andrea Sanchez, David Barrios and Carlos Prieto. Bioinformatics, University of
Salamanca.
}
\value{
This function returns a \code{multinrcor} object.
}
\usage{
multinrcor(data, cor.method = c("pearson", "spearman"), mads = FALSE,
    sample.size = 10, iter = 10000, threads = NULL)
}
\arguments{
\item{data}{a matrix.}
\item{cor.method}{a character string indicating which correlation coefficient
(or covariance) is to be computed. One of "pearson" (default) or "spearman".}
\item{mads}{Should p value be calculated with 'mad' instead of 'sd'?
(default: FALSE)}
\item{sample.size}{a number indicating the sample size (default: 10)}
\item{iter}{The number of iterations (default: 10000)}
\item{threads}{The number of CPU threads for calculating the correlation.
Default value is the amount of CPU cores available on the system.}
}
\description{
\code{multinrcor} compute the correlation of \code{data}.
}
\examples{
obj <- multinrcor(t(iris[,1:4]))
}
