#' @title PerturbationClustering: Perturbation clustering
#' @description Perform subtyping using one type of high-dimensional data
#'
#' @param data input matrix or data frame. The rows represent samples while the columns represent features.
#' @param Kmax the maximum number of clusters. Default value is 10.
#' @param noisePercent the parameter to determine the noise standard deviation. Default is "med", i.e. the noise standard deviation is the medium standard deviation of the features. If noisePercent is numeric, then the noise standard deviation is noisePercent * sd(data).
#' @param iter the number of perturbed datasets. Default value is 200.
#' @param kmIter the number of initial centers used in k-means clustering.
#' @param PCAFunction Custom PCA function for dimension reduction.
#'
#' @details
#'
#' The data are first clustered using k-means. For each value of \emph{k} in the range \emph{[2:Kmax]}, the algorithm buils an original connectivity matrix using the partitioning obtained from k-means. The algorithm then adds Gaussian noise to the data and rebuild the connectivity between samples. For each value of \emph{k}, the algorithm builds \emph{iter} connectivity matrices and then average them to provide one perturbed connectivity matrix.
#'
#' For each value of \emph{k}, the algorithm then constructs a difference matrix, which is the absolute difference between the original and perturbed connectivity matrices for the given \emph{k}. It then calculates the empirical cumulative distribution functions (CDF) for the entries of the difference matrix. The area under the CDF curve (AUC) is used to assess the stability of the clustering. The algorithm chooses the optimal value of \emph{k} for which the AUC value is maximized.
#'
#' It is well known that the k-means algorithm may converge to a local minimum depending on the initialization. To overcome this, the k-means algorithm is run multiple times (using \emph{kmIter} parameter) with randomly chosen seeds and the partitioing with the least residual sum of squares (RSS) is returned.
#' @return
#'
#' \emph{PerturbationClustering} returns a list with at least the following components:
#' \item{k}{The optimal number of clusters}
#' \item{groups}{A vector of labels indicating the cluster to which each sample is allocated}
#' \item{origS}{A list of original connectivity matrices}
#' \item{pertS}{A list of perturbed connectivity matrices}
#'
#' @author
#'
#' Tin Nguyen and Sorin Draghici
#'
#' @references
#'
#' Tin Nguyen, Rebecca Tagett, Diana Diaz, and Sorin Draghici. A novel method for data integration and disease subtyping. Genome Research, 27(12):2025-2039, 2017.
#'
#' @seealso \code{\link{kmeans}}
#'

PINSPerturbationClustering <- function (data, Kmax=10, noisePercent="med", iter=200, kmIter=20, PCAFunction = NULL) {
  if (is.null(PCAFunction)){
    pca <- prcomp(data)
  } else {
    pca <- list(
      x = PCAFunction(data)
    )
  }

  # get the partitioning from simply clustering the real data, for consecutive k
  message("Building original connectivity matrices");flush.console()
  origPartition <- getOriginalSimilarity(data=pca$x, clusRange = 2:Kmax)
  origS <- origPartition$origS

  #noise
  noise = getNoise(data, noisePercent)
  message ("Noise set to ", noise)

  # get perturbed similarity
  message("Building perturbed connectivity matrices");flush.console()
  pertS <- getPerturbedSimilarity(data = pca$x, clusRange=2:Kmax, iter = iter, noiseSd = noise, kmIter = kmIter)

  # get discrepancy
  #message("Calculate discrepancy between original and perturbed connectivity matrices")
  Discrepancy <- getPerturbedDiscrepancy (origS = origS, pertS = pertS, clusRange = 2:Kmax)


  ret <- NULL
  ret$k=min(which(Discrepancy$AUC==max(Discrepancy$AUC[2:Kmax])))
  ret$Kmax=10
  ret$groups <- origPartition$groupings[[ret$k]]
  ret$kmRes <- origPartition$kmRes[[ret$k]] #scPINS
  ret$origS <- origS
  ret$pertS <- pertS
  ret$Discrepancy <- Discrepancy

  message("Done. \n");flush.console()

  ret
}
