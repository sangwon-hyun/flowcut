# Generated from _main.Rmd: do not edit by hand

#' A helper to split the data into (1) censored and (2) non-censored ("inner") particles.
#'
#' @param datobj List containing ylist, countslist, X, and Cbox.
#'
#' @return List of particles on /and/ inside censor boundary. And the directions
#'   in which they were censored.
#'
#' @export
split_data <- function(datobj){
  ## At each time point, set aside censored points.
  censor.direction.full.list <- lapply(datobj$ylist,
                        function(x){
                          flowcut:::censorIndicator(x[,1, drop=FALSE],
                                                    datobj$Cbox[1,,drop=FALSE])})
  censor.which.list = censor.direction.full.list %>%
    lapply(function(x) apply(x,1, function(xx){ sum(abs(xx)) > 0 })) %>%
    lapply(which)
  censored_ylist <- mapply(function(y, censor_which){y[censor_which,,drop=FALSE]},
         datobj$ylist, censor.which.list, SIMPLIFY = FALSE)
  censored_direction <- mapply(function(censor_direction, censor_which){
    censor_direction[censor_which,,drop=FALSE]
  }, censor.direction.full.list, censor.which.list, SIMPLIFY = FALSE)

  ## here are data
  inner_ylist <- mapply(function(y, censor_which){
    if(length(censor_which) > 0) return(y[-censor_which,,drop=FALSE])
    if(!length(censor_which) > 0) return(y)
  }, datobj$ylist, censor.which.list, SIMPLIFY = FALSE)
  inner_countslist <- mapply(function(counts, censor_which){
    if(length(censor_which) > 0) return(counts[-censor_which])
    if(!length(censor_which) > 0) return(counts)
  }, datobj$countslist, censor.which.list, SIMPLIFY = FALSE)

  return(list(censored_ylist = censored_ylist,
              censored_direction = censored_direction,
              inner_ylist = inner_ylist,
              inner_countslist = inner_countslist))
}
