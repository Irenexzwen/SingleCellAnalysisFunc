% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_Summary.R
\name{TopGenes}
\alias{TopGenes}
\title{Title Returen Top expressed genes with highest sum(expr)}
\usage{
TopGenes(COUNT, top = 10)
}
\arguments{
\item{COUNT}{a dataframe with each row represents a gene and each column represents a cell}

\item{top}{int(default 10) number of top genes to inspect}
}
\value{
DataFrame with genenames mean(expr),variance and dropout rate *expr==0/cell_number*
}
\description{
Title Returen Top expressed genes with highest sum(expr)
}
\examples{
tg <- TopGenes(ExprMatx)
}
