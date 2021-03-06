% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/get.genes.R
\name{get.genes.kegg}
\alias{get.genes.kegg}
\title{Downloading Kegg Pathways to find the gene ids}
\usage{
get.genes.kegg(webg.pathways, sleep = 10)
}
\arguments{
\item{webg.pathways}{dataframe where the second column contains the pathway ids}

\item{sleep}{Number of seconds to sleep between Kegg requests}
}
\value{
an array of gene ids
}
\description{
Returns the gene ids
}
\examples{
download.example.data()
webg.pathways = load.WebGestalt("Marra_0_kegg_protein_enrichment.tsv",'Kegg')
gene.ids = get.genes.kegg(webg.pathways)
}

