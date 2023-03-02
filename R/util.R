#' @import data.table
#' @import abind
#' @useDynLib metapop, .registration=TRUE
NULL

get_lmm <- function(params){

  lmm <- matrix(, nrow=params$n*params$n_vac,ncol=params$n*params$n_vac)
                
  for(i in 1:params$n_vac){
    for(j in 1:params$n_vac){
      lmm[ ((i-1)*params$n):(i*params$n-1)+1, ((j-1)*params$n):(j*params$n-1)+1] <- params$mixing_matrix
    }
  }
  return(lmm)
}

