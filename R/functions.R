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