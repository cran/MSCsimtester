% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MSCsimtester.R
\docType{package}
\name{MSCsimtester}
\alias{MSCsimtester}
\title{Validity tests of simulators of the multspecies coalescent model in phylogenomics.}
\description{
The package performs comparisons of certain summary statistics for simulated gene tree samples to
theoretical predictions under the multispecies coalescent model.
The primary functions are \code{rootedTriple} for comparison of frequencies of topological rooted triples
on gene trees, and \code{pairwiseDist} and \code{ADtest} for comparison of the distributions of pairwise
distances between taxa on gene trees.
}
\details{
Required input is a collection of gene trees, stored as a \code{multiPhylo} object by the \code{ape} package,
and a rooted species tree, as a \code{Phylo} object,
with edge lengths in generations, together with constant population sizes for each edge.

\code{MSCsimtester} builds on the packages \code{ape} and \code{kSamples}.

For further examples of use and citation purposes, see \insertCite{testing19}{MSCsimtester}.
}
\references{
\insertRef{testing19}{MSCsimtester}
}
