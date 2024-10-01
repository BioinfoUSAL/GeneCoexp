\name{create_network}
\alias{create_network}
\title{Create a Network}
\author{
Andrea Sanchez, David Barrios and Carlos Prieto. Bioinformatics, University of
Salamanca.
}
\value{
This function returns a \code{rD3plot} object.
}
\usage{
create_network(x, cutoff = NULL, cutoffvar = "padjust",
    nodes = NULL, name = NULL, label = NULL)
}
\arguments{
\item{x}{a \code{multinrcor} object.}
\item{cutoff}{filter links by the column especified by \code{cutoffvar}.}
\item{cutoffvar}{Which column should be filtered by cutoff?
(Default: padjust)}
\item{nodes}{a data frame with at least one column of node names. This is
useful for assigning attributes to each node.}
\item{name}{name of the column with names in the nodes data frame. It will take
the first nodes column by default.}
\item{label}{name of the column with labels in the nodes data frame.}
}
\description{
\code{create_network} creates a network from a \code{multinrcor} object.
}
\examples{
obj <- multinrcor(t(iris[,1:4]))
net <- create_network(obj)
}
