% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PerturbationClustering.R
\name{PINSPerturbationClustering}
\alias{PINSPerturbationClustering}
\title{PerturbationClustering: Perturbation clustering}
\usage{
PINSPerturbationClustering(
  data,
  Kmax = 10,
  noisePercent = "med",
  iter = 200,
  kmIter = 20,
  PCAFunction = NULL
)
}
\arguments{
\item{data}{input matrix or data frame. The rows represent samples while the columns represent features.}

\item{Kmax}{the maximum number of clusters. Default value is 10.}

\item{noisePercent}{the parameter to determine the noise standard deviation. Default is "med", i.e. the noise standard deviation is the medium standard deviation of the features. If noisePercent is numeric, then the noise standard deviation is noisePercent * sd(data).}

\item{iter}{the number of perturbed datasets. Default value is 200.}

\item{kmIter}{the number of initial centers used in k-means clustering.}

\item{PCAFunction}{Custom PCA function for dimension reduction.}
}
\value{
\emph{PerturbationClustering} returns a list with at least the following components:
\item{k}{The optimal number of clusters}
\item{groups}{A vector of labels indicating the cluster to which each sample is allocated}
\item{origS}{A list of original connectivity matrices}
\item{pertS}{A list of perturbed connectivity matrices}
}
\description{
Perform subtyping using one type of high-dimensional data
}
\details{
The data are first clustered using k-means. For each value of \emph{k} in the range \emph{[2:Kmax]}, the algorithm buils an original connectivity matrix using the partitioning obtained from k-means. The algorithm then adds Gaussian noise to the data and rebuild the connectivity between samples. For each value of \emph{k}, the algorithm builds \emph{iter} connectivity matrices and then average them to provide one perturbed connectivity matrix.

For each value of \emph{k}, the algorithm then constructs a difference matrix, which is the absolute difference between the original and perturbed connectivity matrices for the given \emph{k}. It then calculates the empirical cumulative distribution functions (CDF) for the entries of the difference matrix. The area under the CDF curve (AUC) is used to assess the stability of the clustering. The algorithm chooses the optimal value of \emph{k} for which the AUC value is maximized.

It is well known that the k-means algorithm may converge to a local minimum depending on the initialization. To overcome this, the k-means algorithm is run multiple times (using \emph{kmIter} parameter) with randomly chosen seeds and the partitioing with the least residual sum of squares (RSS) is returned.
}
\references{
Tin Nguyen, Rebecca Tagett, Diana Diaz, and Sorin Draghici. A novel method for data integration and disease subtyping. Genome Research, 27(12):2025-2039, 2017.
}
\seealso{
\code{\link{kmeans}}
}
\author{
Tin Nguyen and Sorin Draghici
}
