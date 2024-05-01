#' @title Find the maximum Bernoulli log-likelihood parameter
#'
#' @description This function calculates the log-likelihood of a Bernoulli distribution for a range of probability values (p) and then finds the p that maximizes the log-likelihood for the given data.
#'
#' @param data A numeric vector containing 0s and 1s representing Bernoulli trials.
#'
#' @return The value of p that maximizes the log-likelihood of the Bernoulli distribution for the given data.
#'
#' @details 
#' This function leverages the idea that a Bernoulli distribution is a special case of the Binomial distribution where the number of trials (n_trials) is 1. It calculates the log-likelihood using `dbinom` for a range of p values and then identifies the p that leads to the highest log-likelihood.
#'
#' @examples
#' 
#' # Simulate some Bernoulli trial data
#' data <- rbinom(n = 100, size = 1, prob = 0.6)
#' 
#' # Find the maximum log-likelihood parameter
#' estimated_p <- logLikBernoulli(data)
#' 
#' print(estimated_p)
#'
#' @export
#' 
logLikBernoulli = function(data){
  BernouilliLogLik = function(p, x){
    # Using the idea that a Bernoulli distribution is a special case of the Binomial distribution where n_trials = 1
    loglik = sum(log(dbinom(x, size=1 ,prob=p)))
    return (loglik)
  }
  
  p = seq(from=0, to=1, by=0.1)
  start = 2
  end = length(p) - 1
  
  df = tibble(p=p[start:end])
  df$loglik = apply(matrix(data_p$p), 1, BernouilliLogLik, x=data)
  
  max_row = df |> arrange(desc(loglik)) |> head(1)
  
  return (max_row$p)
}


#' @title Unscale a vector that was previously scaled with `scale()`
#'
#' @description This function takes a numeric vector that was previously scaled using the `scale()` function and attempts to unscale it back to the original scale. It retrieves the centering and scaling factors stored as attributes by `scale()` and applies the inverse transformation.
#'
#' @param x A numeric vector that was previously scaled with `scale()`.
#'
#' @return The unscaled numeric vector.
#'
#' @details This function assumes that the input vector was scaled using the `scale()` function. It retrieves the centering (`scaled:center`) and scaling (`scaled:scale`) factors stored as attributes by `scale()`. If neither attribute is found, it raises an error indicating the vector cannot be unscaled. Otherwise, it applies the inverse transformation based on the presence of both or only one of the attributes.
#'
#' @note This function only works for vectors that were scaled with `scale()`. It might not be suitable for other scaling methods.
#'
#' @examples
#' 
#' # Simulate some data
#' data <- rnorm(n = 10)
#' 
#' # Scale the data
#' scaled_data <- scale(data)
#' 
#' # Unscale the data
#' unscaled_data <- unscale(scaled_data)
#' 
#' @export
unscale = function(x){
  x_center = attributes(x)$`scaled:center`
  x_scale = attributes(x)$`scaled:scale`
  
  if (is.null(x_center) &  is.null(x_scale)){
    stop("Error: The vector does not contain a center or scale attribute, hence it cannot be unscaled")
  } 
  else{
    if(is.null(x_scale)){
      x_raw = (x) + x_center 
    }
    else{
      if (is.null(x_center)){
        x_raw = (x * x_scale)
      }
      else {
        x_raw = (x * x_scale) + x_center      
      }
    }
  }
  
  x_unscaled = as.vector(x_raw)
  return (x_unscaled)
}

#' @title Perform Principal Component Approximation (PCA)
#'
#' @description This function performs Principal Component Analysis (PCA) on a data matrix and returns an approximation of the data using a specified number of principal components (PCs). PCA is a dimensionality reduction technique that identifies the directions of greatest variance in the data. This function centers the data and then computes the covariance matrix to find the principal components. It retains the specified number of PCs and uses them to reconstruct an approximation of the original data.
#'
#' @param x A numeric data matrix.
#'
#' @param npc An integer specifying the number of principal components to retain for approximation (default is 1).
#'
#' @return A numeric matrix representing the approximation of the original data using the specified number of principal components.
#'
#' @details This function assumes the data in `x` is in its original scale. It centers the data using `scale(center=TRUE, scale=FALSE)` before performing PCA. The function retrieves the eigenvectors corresponding to the largest eigenvalues from the covariance matrix. It then uses the selected eigenvectors to reconstruct an approximation of the centered data. Finally, it uncenters the reconstructed data to obtain the approximation in the original scale (assuming attributes are present).
#'
#' @note This function performs PCA centering but not scaling. Scaling might be beneficial depending on the data characteristics.
#'
#' @examples
#' 
#' # Simulate some data
#' data <- matrix(rnorm(n = 20, sd = 2), nrow = 10)
#' 
#' # Perform PCA approximation with 2 principal components
#' pca_approx <- pcApprox(data, npc = 2)
#' 
#' @export
pcApprox = function(x, npc){
  # Need to determine whether x is a vector or matrix
  # Center and or scale the data matrix. I think PCA only requires you to center the data matrix
  x_centered = scale(x, center=TRUE, scale=FALSE)
  
  # Compute the covariance matrix from the data
  cov_x = cov(x_centered)
  
  # Get the eigen vectors of the covariance matrix
  pcs = eigen(cov_x)
  
  pc_vectors = pcs$vectors
  print(pc_vectors)
  # To get the approximation of the data, select the n columns of the eigen vector matrix. I think R gives 
  # them ordered in descending order
  pca_vectors = pc_vectors[,1:npc]
  
  # Compute the reconstruction of the centered data. The computation is done on the centered datas
  x_reconstruct = x_centered %*% pca_vectors %*% t(pca_vectors) 
  
  # Uncenter and/or unscale the data to return the approximation of the data.
  x_result = x_reconstruct + attributes(x_centered)$`scaled:center`
  
  return (x_result)
}






