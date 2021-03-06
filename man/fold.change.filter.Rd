% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/filters.R
\name{fold.change.filter}
\alias{fold.change.filter}
\title{Fold change filter}
\usage{
fold.change.filter(data, theta, class1.inxs, class2.inxs)
}
\arguments{
\item{data}{- transposed gene data}

\item{theta}{- if > 1, theta is the target number of ranked variables;
otherwise, it is the top theta fraction of ranked variables.}

\item{class1.inxs}{- the row indices of class 1}

\item{class2.inxs}{- the row indices of class 2}
}
\value{
filtered dataset
}
\description{
Applies a simple fold change filter on the data.
}
\examples{
download.example.data()
webg.pathways = load.WebGestalt("Marra_0_wiki_protein_enrichment.tsv",'Wiki')
gene.ids = get.genes.wiki(webg.pathways)
Cattaneo.rna = load.gene.data("Cattaneo_array.csv",5)
prioritized.data = list.filter(Cattaneo.rna$transposed.data,gene.ids)
prioritized.filtered.data = fold.change.filter(prioritized.data,100,
                                               which(prioritized.data$labels=='adenoma'),
                                               which(prioritized.data$labels=='normal'))
}

