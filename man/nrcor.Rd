\name{nrcor}
\alias{nrcor}
\title{Correlation}
\author{
Andrea Sanchez, David Barrios and Carlos Prieto. Bioinformatics, University of
Salamanca.
}
\value{
This function returns a \code{nrcor} object.
}
\usage{
nrcor(x, y, r.value = NA, cor.method = c("pearson", "spearman"),
    mads = FALSE, sample.size = 10, iter = 10000, samples_results = TRUE)
}
\arguments{
\item{x}{a vector.}
\item{y}{a matrix.}
\item{r.value}{a number.}
\item{cor.method}{a character string indicating which correlation coefficient
(or covariance) is to be computed. One of "pearson" (default) or "spearman".}
\item{mads}{Should p value be calculated with 'mad' instead of 'sd'?
(default: FALSE)}
\item{sample.size}{a number.}
\item{iter}{The number of iterations (default: 10000)}
\item{samples_results}{Should 'samplesmat' and 'samplesnum' be calculated?
(Default: TRUE)}
}
\description{
\code{nrcor} compute the correlation of \code{x} and \code{y}.
}
\examples{
res <- nrcor(iris[,5],t(iris[,1:4]))
}
