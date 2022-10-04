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


#'
#' @export
fix_beta_mode_params <- function(params){

  params_keys <- names(params)
  if(!"rand_beta_sd" %in% params_keys) params$rand_beta_sd <- 0.1
  if(!"rand_beta_factors" %in% params_keys) params$rand_beta_factors <-rep(1, params$n)
  if(!"log_beta_ini" %in% params_keys) params$log_beta_ini <-0
  if(!"beta_cut_peak_param" %in% params_keys) params$beta_cut_peak_param <- rep(0, 4)
  if(!"change_factor" %in% params_keys) params$change_factor <- rep(0, 4)
  if(!"threshold" %in% params_keys) params$threshold <- c(100,0)
  if(!"spont_behav_change_params" %in% params_keys) params$spont_behav_change_params <- c(1,100000, 1, 1)
  if(!"expected_health_loss" %in% params_keys) params$expected_health_loss <- array(0, dim=c(params$n,params$n_vac))
  return(params)
}
  
