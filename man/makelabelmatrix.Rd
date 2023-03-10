% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makelabelmaatrix.R
\name{makelabelmatrix}
\alias{makelabelmatrix}
\title{Label matrix}
\usage{
makelabelmatrix(clust)
}
\arguments{
\item{clust}{An n-dimensional vector of class partition.}
}
\value{
\item{Z}{ A matrix of class indicator.}
}
\description{
Convert class indicator into a label maxtrix.
}
\examples{
cluster<-c(1,1,2,2,3,3)
label_maxtrix<-makelabelmatrix(cluster)
}
