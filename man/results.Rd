\name{results}
\alias{results}
\title{Create a table of results}
\author{
Andrea Sanchez, David Barrios and Carlos Prieto. Bioinformatics, University of
Salamanca.
}
\value{
This function returns a data frame.
}
\usage{
results(x, cutoff = NULL, cutoffvar = NULL)
}
\arguments{
\item{x}{a \code{multinrcor} or \code{nrcor} object.}
\item{cutoff}{filter links by the column especified by \code{cutoffvar}.}
\item{cutoffvar}{Which column should be filtered by cutoff?}
}
\description{
\code{results} creates a table from a \code{multinrcor} or \code{nrcor} object.
The table will be ordered by \code{cutoffvar} or \code{N} by default.
}
\examples{
obj <- multinrcor(t(iris[,1:4]))
net <- results(obj)
}
