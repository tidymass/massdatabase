% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/19_REACTOME.R
\name{request_reactome_reaction_info}
\alias{request_reactome_reaction_info}
\title{Retrieve Reactome Reactions for a Given Organism}
\usage{
request_reactome_reaction_info(organism = "Homo sapiens")
}
\arguments{
\item{organism}{Character. The name of the organism. Default is "Homo sapiens".}
}
\value{
A data frame containing:
  \item{reaction_id}{Character. The Reactome reaction ID.}
  \item{reaction_link}{Character. The URL link to the reaction in Reactome.}
  \item{reaction_name}{Character. The name of the reaction.}
  \item{organism}{Character. The organism associated with the reaction.}
}
\description{
This function retrieves reaction information for a specified organism from Reactome.
}
\examples{
# Retrieve reaction information for Homo sapiens
reactions <- request_reactome_reaction_info("Homo sapiens")
head(reactions)
}
\author{
Xiaotao Shen (\email{xiaotao.shen@outlook.com})
}
