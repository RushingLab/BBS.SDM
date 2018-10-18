#' RemoveOutliers
#'
#' Removes routes that are likely outliers (far outside "likely" breeding range)
#' 
#' Estimates for each route \code{i} the mean distance of each route to its \code{k} nearest neighbors (knn-distance). Any route with knn-distance greater than \code{tresh}*sd above the mean(knn-distance) is considered an outlier and removed from the data 
#' @param counts Species count data frame obtained from GetSppCounts
#' @param thresh Number of standard deviations beyond mean knn distance to be considered outlier (default = 5)
#' @param k Number k on nearest neighbors
#' @return Data frame with same format as \code{counts} but with outlier routes removed
#' @export


RemoveOutliers <- function(raw.counts = NULL, thresh = 15, k = 1, Write = TRUE, path = NULL){

  if(is.null(raw.counts)){
    counts <- read.csv(paste0(path, "/raw_counts.csv"))
  }else{
    counts <- raw.counts
  }
  
  route_xy <- as.matrix(counts[!duplicated(counts$routeID), c("routeID", "Longitude", "Latitude")])
  dist.mat <- geosphere::distm(route_xy[,c("Longitude", "Latitude")])
  nn <- apply(dist.mat, 1, function(x) sort(x)[2:(k + 1)])
  if(k == 1){
    nn.dist <- nn
  }else{
    nn.dist <- colMeans(nn)
  }
  mu.nn <- mean(nn.dist)
  sd.nn <- sd(nn.dist)
  
  cutoff <- mu.nn + thresh * sd.nn
  keep <- route_xy[which(nn.dist < cutoff), "routeID"]
  counts2 <- dplyr::filter(counts, routeID %in% keep)
  
  if(Write){
    if(is.null(path)){
      write.csv(counts2,
                "no_outlier_counts.csv",
                row.names = FALSE)
    }else{
      write.csv(counts2,
                paste(path, "no_outlier_counts.csv", sep = "/"),
                row.names = FALSE)
    }
    
  }else{
    return(counts2)
  }
}
