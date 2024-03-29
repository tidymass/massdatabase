% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/16_KEGG.R
\name{request_kegg_drug}
\alias{request_kegg_drug}
\title{Request one specific drug information in KEGG}
\usage{
request_kegg_drug(drug_id = "D00001", return_form = c("list", "data.frame"))
}
\arguments{
\item{drug_id}{drug id. For example, D00001}

\item{return_form}{data.frame or list.}
}
\value{
A data frame or list.
}
\description{
Request one specific drug information in KEGG
}
\examples{
x = request_kegg_drug(drug_id = "D00001", return_form = "list")
x[1:2]
y = request_kegg_drug(drug_id = "D00001", return_form = "data.frame")
head(y)
}
\author{
Xiaotao Shen
\email{shenxt1990@outlook.com}
}
