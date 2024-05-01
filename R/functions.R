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
  
  df = dplyr::tibble(p=p[start:end])
  df$loglik = apply(matrix(df$p), 1, BernouilliLogLik, x=data)
  
  max_row = df |> dplyr::arrange(desc(loglik)) |> utils::head(1)
  
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


#' @title Plot a survival curve S(t)
#'
#' @description This function takes vectors of event times (`time`) and event indicators (`status`) and calculates the Kaplan-Meier (KM) survival curve. The KM curve estimates the probability of surviving beyond a specific time point for a population experiencing events over time.
#'
#' @param status A numeric vector indicating event status (1 = event, 0 = censored).
#'
#' @param time A numeric vector of event times corresponding to the `status` vector.
#'
#' @return A ggplot object representing the Kaplan-Meier survival curve.
#'
#' @details This function utilizes dplyr verbs to manipulate the data and calculate the Kaplan-Meier curve elements. It first arranges the data by time. Then, it groups by time, calculates the number of events, censored observations, total at risk, hazard rate, and survival probability for each time point. Missing hazard rates (at the beginning) are replaced with 0. Finally, it accumulates survival probabilities, ungroups the data, adds an initial point (time=0, survival=1), rearranges by time, and selects the desired columns (time and survival) for the plot. The function then uses ggplot2 to create the survival curve visualization.
#'
#' @references Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from incomplete observations. Journal of the American Statistical Association, 53(282), 457-481.
#'
#' @examples
#' 
#' # Simulate some survival data
#' set.seed(123)
#' time <- rexp(n = 20, rate = 0.1)
#' status <- rbinom(n = 20, size = 1, prob = 0.7)
#' 
#' # Generate Kaplan-Meier survival curve
#' km_curve <- survCurv(status, time)
#' 
#' # Print the ggplot object (use ggplot2::print() for detailed view)
#' km_curve
#'
#' @export

survCurv = function(status, time){
  df = dplyr::tibble(time=time, status=status) |> dplyr::arrange(time)
  
  df.curv = df |>
    dplyr::group_by(time) |>
    dplyr::summarize(n.events = sum(status == 1), n.censored = sum(status == 0)) |>
    dplyr::mutate(n.tot = n.events + n.censored) |>
    dplyr::mutate(n.at.risk = lag(nrow(df)-cumsum(n.tot), default=nrow(df))) |>
    dplyr::mutate(Hazard = n.events/n.at.risk) |>
    dplyr::mutate(Hazard = tidyr::replace_na(Hazard,0)) |>
    dplyr::mutate(Survival = cumprod(1-Hazard)) |>
    dplyr::ungroup() |>
    dplyr::bind_rows(dplyr::tibble(time=0,Survival=1)) |>
    dplyr::arrange(time) |>
    dplyr::select(time, Survival)
  
  fig = df.curv |> ggplot2::ggplot(ggplot2::aes(time, Survival)) + ggplot2::geom_line() + ggplot2::labs(x = "Time", y = "Probability of survival")
  
  return (fig)
}


#' @title Downloads a REDCap report as a data frame
#'
#' @description This function downloads a REDCap report as a data frame using an access token retrieved from an environment variable. It constructs a form data object with the necessary parameters for the REDCap API and sends a POST request to the REDCap URL. The function then parses the response content and converts it to a data frame using `dplyr::tibble`.
#'
#' @param redcapTokenName A character string representing the environment variable name that stores the REDCap access token.
#'
#' @param redcapUrl A character string representing the base URL of the REDCap instance.
#'
#' @param redcapReportId A character string representing the unique identifier of the REDCap report to download.
#'
#' @return A data frame containing the downloaded report data (if successful). The structure and content of the data frame will depend on the specific REDCap report design.
#'
#' @details This function assumes the REDCap access token is stored in an environment variable. It retrieves the token using `Sys.getenv` and constructs a form data object with the required parameters for the REDCap API's download functionality. The function utilizes the `httr` package for making the HTTP POST request and parsing the response. It's important to ensure you have the necessary permissions and configurations set up on the REDCap server for accessing reports through the API.
#'
#' @note Downloading reports via the API might have limitations compared to the REDCap user interface. It's recommended to consult the REDCap API documentation for specific details and potential restrictions.
#'
#' @references REDCap Consortium. (2023). REDCap API Documentation. https://redcap.vanderbilt.edu/
#'
#' @examples
#' 
#' # Assuming REDCap access token is stored in an environment variable named "REDCAP_TOKEN"
#' report_data <- downloadRedcapReport(redcapTokenName = "REDCAP_TOKEN", 
#'                                     redcapUrl = "https://your_redcap_url.com/api/", 
#'                                     redcapReportId = "123")
#'
#' # Explore the downloaded data frame (structure may vary depending on the report)
#' head(report_data)
#'
#' @export
downloadRedcapReport = function(redcapTokenName, redcapUrl, redcapReportId){
  token <- Sys.getenv(paste(redcapTokenName))
  url <- redcapUrl

  formData <- list("token"=token,
                   content='report',
                   format='csv',
                   report_id=redcapReportId,
                   csvDelimiter='',
                   rawOrLabel='raw',
                   rawOrLabelHeaders='raw',
                   exportCheckboxLabel='false',
                   returnFormat='csv'
  )
  
  response <- httr::POST(url, body = formData, encode = "form")
  result <- httr::content(response)
  df_result = dplyr::tibble(result)
  
  return (df_result)
}

#' Calculate Minimum Sample Size for T-Test
#'
#' This function calculates the minimum sample size needed for a t-test with specified power and significance level.
#'
#' @param x1 Numeric vector. One or two samples of preliminary data.
#' @param x2 Numeric vector (optional). Second sample of preliminary data for a two-sample t-test.
#' @param power Numeric scalar. Desired statistical power (default is 0.8).
#' @param alpha Numeric scalar. Significance level (default is 0.05).
#'
#' @return Numeric scalar. Minimum sample size needed for the t-test.
#'
#' @details This function calculates the minimum sample size needed for a t-test, either one-sample or two-sample, with specified power and significance level. If only one sample is provided (\code{x1}), a one-sample t-test is assumed; otherwise, a two-sample t-test is assumed.
#'
#' @examples
#' # One-sample t-test
#' x <- rnorm(20)
#' minimumN(x)
#'
#' # Two-sample t-test
#' x1 <- rnorm(20)
#' x2 <- rnorm(20)
#' minimumN(x1, x2)
#'
#' @export
minimumN <- function(x1, x2 = NULL) {
  if (is.null(x2)) {
    # One-sample t-test
    result <- pwr::pwr.t.test(d = mean(x1), sig.level = 0.05, power = 0.8, alternative = "two.sided")
    return(ceiling(result$n))
  } else {
    # Two-sample t-test
    result <- pwr::pwr.t2n.test(n1 = length(x1), n2 = length(x2), sig.level = 0.05, power = 0.8, alternative = "two.sided")
    return(ceiling(result$n))
  }
}


#' @title Standardize column names to small camel case
#'
#' @description This function standardizes the column names of a data frame to small camel case (e.g., "columnA" becomes "columnA"). It utilizes `dplyr::rename_with` and `janitor::make_clean_names` to achieve consistent and valid variable names.
#'
#' @param data A data frame object.
#'
#' @return A data frame with the same data but standardized column names in small camel case.
#'
#' @details This function leverages the `dplyr` and `janitor` packages. It applies `janitor::make_clean_names` to each column name, cleaning and converting it to a valid R variable name in small camel case. Then, it uses `dplyr::rename_with` to perform the renaming based on the standardized names.
#'
#' @examples
#' 
#' # Example data frame
#' data <- data.frame(column_A = 1:5, Column_B = rnorm(5))
#' 
#' # Standardize column names to small camel case
#' standardized_data <- standardizeNames(data)
#' 
#' # Print the data frame with standardized names
#' head(standardized_data)
#'
#' @export
standardizeNames = function(data){
  data = data |> dplyr::rename_with(~ janitor::make_clean_names(.x, case = "small_camel"))
  return (data)
}


