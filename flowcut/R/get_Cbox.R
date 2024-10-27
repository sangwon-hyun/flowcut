# Generated from _main.Rmd: do not edit by hand

#' Helper to obtain the censoring limits |Cbox|, which is a 2-colum matrix with columns |bounds.lower| and |bounds.upper|.
#' 
#' @param ylist List of d-dimensional cytogram data.
#'
#' @export
#' @return dimdat by 2 matrix
get_Cbox <- function(ylist){
  bounds.lower <- lapply(ylist,function(xx)
      Rfast::colMins(xx,value = TRUE)) %>%
      do.call(rbind,.) %>% Rfast::colMins(.,value = TRUE)
  bounds.upper <- lapply(ylist,function(xx)
      Rfast::colMaxs(xx,value = TRUE)) %>%
      do.call(rbind,.) %>% Rfast::colMaxs(.,value = TRUE)
  Cbox <- cbind(bounds.lower,bounds.upper)
}
