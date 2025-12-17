#' ExpressionCellNet: Network-centric modeling of gene expression perturbations
#'
#' ExpressionCellNet provides an integrative workflow for constructing biologically
#' constrained gene–gene interaction networks by combining curated interaction evidence
#' with statistically supported expression correlations. Built networks can be analyzed
#' using centrality metrics, hub-gene identification, shortest-path queries, and
#' enrichment analysis. The framework also supports regression-based prediction to
#' simulate in silico perturbations (e.g., gene silencing/overexpression) and to
#' estimate how expression changes propagate through connected genes.
#'
#' @section Core workflow:
#' \enumerate{
#'   \item Create object: \code{createExpCellNetObj()}
#'   \item Normalize/filter: \code{CountMatrixNormalization()}
#'   \item Build network: \code{BuildNetwork()}
#'   \item Network analysis: \code{AnalyzeNetwork()}, \code{ShowHubGenes()}, \code{FindPathway()}
#'   \item Prediction and visualization: \code{NetworkPrediction()}, \code{PlotNetworkPrediction()}
#'   \item Functional interpretation: \code{GeneEnrichment()}, \code{PlotEnrichment()}
#'   \item Dimension reduction: \code{PCAforGenes()}, \code{UMAPforGenes()}
#' }
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{createExpCellNetObj}}{Create the main analysis object.}
#'   \item{\code{CountMatrixNormalization}}{Filter and normalize count matrices (CPM/log2).}
#'   \item{\code{BuildNetwork}}{Construct a multi-generation network from a seed gene.}
#'   \item{\code{AnalyzeNetwork}}{Compute centrality metrics and extract hub genes (top 10\%).}
#'   \item{\code{FindPathway}}{Shortest path between two genes in the constructed network.}
#'   \item{\code{NetworkPrediction}}{Predict expression response to gene perturbation.}
#'   \item{\code{GeneEnrichment}}{Pathway enrichment on network genes (or subsets).}
#' }
#'
#' @section Input data:
#' Required inputs typically include:
#' \itemize{
#'   \item A gene-by-sample expression matrix (counts or expression values)
#'   \item A gene annotation table with consistent gene identifiers
#'   \item A curated gene–gene interaction database (two-column edge list)
#' }
#'
#' @section Notes:
#' This package is intended for decision support and hypothesis generation and does not
#' claim causal inference.
#'
#' @docType package
#' @name ExpressionCellNet
#' @keywords internal
NULL
