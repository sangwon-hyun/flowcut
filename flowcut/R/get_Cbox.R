# Generated from _main.Rmd: do not edit by hand

#' Helper to obtain the censoring limits |Cbox|, which is a 2-colum matrix with columns |bounds.lower| and |bounds.upper|.
#'
#' @param ylist Data (list of d-column matrices) including censored particles.
#'
#' @return |dimdat| by 2 matrix.
#' @export
get_Cbox <- function(ylist){

  ## Setup
  dimdat = ncol(ylist[[1]])

  ## Define censor limits
  bounds.lower <- Rfast::rowMins(matrix(unlist(
      lapply(ylist, function(xx) Rfast::colMins(xx,value = TRUE))), nrow = dimdat), value = TRUE)
  bounds.upper <- Rfast::rowMaxs(matrix(unlist(
      lapply(ylist, function(xx) Rfast::colMaxs(xx,value = TRUE))), nrow = dimdat), value = TRUE)

  ## bounds.lower <- lapply(ylist,function(xx)
  ##     Rfast::colMins(xx,value = TRUE)) %>%
  ##     do.call(rbind,.) %>% Rfast::colMins(.,value = TRUE)
  ## bounds.upper <- lapply(ylist,function(xx)
  ##     Rfast::colMaxs(xx,value = TRUE)) %>%
  ##     do.call(rbind,.) %>% Rfast::colMaxs(.,value = TRUE)

  ## Modify censorship slightly
  Cbox <- cbind(bounds.lower, bounds.upper)
}
