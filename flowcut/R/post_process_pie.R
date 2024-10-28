# Generated from _main.Rmd: do not edit by hand

#' Description goes here
#' 
#' @param ga document this carefully
#' @param Xp 1 vector appended to left of covariate matrix X.
#' 
#' @export
#' @return (OO x OO) matrix with each row..
post_process_pie <- function(Xp,ga){
  Xga <- Xp %*% ga 
  pie.SB <- 1/(1+exp(-Xga))
  return(t(apply(pie.SB,1, SB2MN)))
}
