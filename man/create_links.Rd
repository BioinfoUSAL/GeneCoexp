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
create_links(x, cutoff = NULL)
}
\arguments{
\item{x}{a \code{multinrcor} object.}
\item{cutoff}{filter links by \code{padjust} value.}
}
\description{
\code{create_links} extracts the links from a \code{multinrcor} object.
}
\examples{
obj <- multinrcor(t(iris[,1:4]))
links <- create_links(obj)
}
