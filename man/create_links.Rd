\name{create_links}
\alias{create_links}
\title{Create Links}
\author{
Andrea Sanchez, David Barrios and Carlos Prieto. Bioinformatics, University of
Salamanca.
}
\value{
This function returns a data frame of links.
}
\usage{
create_links(x, cutoff = NULL, cutoffvar = "padjust")
}
\arguments{
\item{x}{a \code{multinrcor} object.}
\item{cutoff}{filter links by the column especified by \code{cutoffvar}.}
\item{cutoffvar}{Which column should be filtered by cutoff?
(Default: padjust)}
}
\description{
\code{create_links} extracts the links from a \code{multinrcor} object.
}
\examples{
obj <- multinrcor(t(iris[,1:4]))
links <- create_links(obj)
}
