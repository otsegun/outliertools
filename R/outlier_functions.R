


# Define Comedian Sajesh ---------------------------------------------------------------------------
#require(robustbase)
comsajesh <- function(dt){
  #set number of rows and number of columns
  p <- ncol(dt)
  # get a robust estimate of comedian location and scatter from robustbase
  comsajesh_estimate <- robustbase::covComed(dt, n.iter = 2, reweight = T)
  # set cutoff
  cutoff <- 1.4826 * (qchisq(p = 0.975, df = p) * median(comsajesh_estimate$mah))/qchisq(p = 0.5, df = p)
  # find outliers
  outlier_indices <-  which(comsajesh_estimate$mah > cutoff)
  return(outlier_indices)
}




# Define adjusted mcd function -----------------------------------------------------------------------
#library(rrcov)
mcd_adjusted_quantile <- function(dt){
  #set number of rows and number of columns
  n <- nrow(dt)
  p <- ncol(dt)
  # get a robust estimate of location and scatter using the fast mcd algorithm
  # recommended h = [(n+p+1)/2] for breakdown value of  (n-h+1)/n
  mcd_estimate <- rrcov::CovMcd(dt, alpha = ((n + p + 1)/2)/n)
  # set cutoff
  cutoff <- adjusted_quantile(dt = dt, distances = mcd_estimate@mah)

  # find outliers
  outlier_indices <-  which(mcd_estimate@mah > cutoff)
  #  n_outliers <- sum(mcd_dist >= cutoff)

  return(outlier_indices)

}


## function to identify outliers based on the mcd but using
## the adjusted cutoff of


# adjusted threshold is estimated adaptively from the data, defined for robust MCD Mahalanobis distances,
#Filmozer et al. [2005]
adjusted_quantile <- function(dt, distances){

  n <- nrow(dt)
  p <- ncol(dt)
  #value of delta
  deltan <- qchisq(p = 0.975, df = p)

  # order the distances
  ord_dist <- sort(distances)

  # calculate the empirical distribution of the values of the squared distances
  pchisq_ord_dist <- pchisq(ord_dist,p)

  # find the difference between the empirical distribution and theoretical distribution
  dif <- pchisq_ord_dist - ((.5:n)/n)

  # first condition: if ordered distance more than the value of deltan
  # Second condition: if the difference values more than zero
  cond <- (ord_dist >= deltan) & (dif > 0)
  alphan <-  ifelse(sum(cond) == 0, 0, max(dif[cond]) )
  return(max(ord_dist[n - ceiling(n*alphan)], deltan))
}


# Define mcd function ----------------------------------------------------

#require(rrcov) # for fast mcd algortithm
mcd <- function(dt){
  #set number of rows and number of columns
  p <- ncol(dt)
  # get a robust estimate of location and scatter using the fast mcd algorithm
  # recommended h = [(n+p+1)/2] for breakdown value of  (n-h+1)/n
  mcd_estimate <- rrcov::CovMcd(dt, alpha = 0.5)
  # set cutoff
  cutoff <- sqrt(qchisq(p = 0.975, df = p))

  # find outliers
  outlier_indices <-  which(sqrt(mcd_estimate@mah) > cutoff)
  return(outlier_indices)

}


# Define ogk function --------------------------------------------------

#require(rrcov) # for ogk estimate algortithm
ogk <- function(dt){
  #set number of rows and number of columns
  n <- nrow(dt)
  p <- ncol(dt)
  # get a robust estimate of location and scatter using the fast mcd algorithm
  # recommended h = [(n+p+1)/2] for breakdown value of  (n-h+1)/n
  ogk_estimate <- rrcov::CovOgk(dt)

  # compute the robust distance
  ogk_dist <- sqrt(mahalanobis(dt, ogk_estimate@center, ogk_estimate@cov))

  # set cutoff
  cutoff <- sqrt(qchisq(p = 0.975, df = p))

  # find outliers
  outlier_indices <-  which(ogk_dist > cutoff)
  #  n_outliers <- sum(mcd_dist >= cutoff)
  return(outlier_indices)
}
