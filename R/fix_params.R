
#' Combining vaccination parameters
#'
#' @export
fix_params <- function(params, N, n_vac, n_strain, vac_pars){
  params$asympt_frac <- array(params$asympt_frac, dim=c(N, n_vac, n_strain))
  params$sympt_frac <- array(params$sympt_frac, dim=c(N, n_vac, n_strain))
  params$length_hosp <- array(outer(params$length_hosp, vac_pars$rr_los_hosp), dim=c(N, n_vac, n_strain))
  params$length_icu <- array(params$length_icu, dim=c(N, n_vac, n_strain))
  params$hosp_prob <- array(outer(params$hosp_prob, vac_pars$rr_hosp), dim=c(N, n_vac, n_strain))
  params$icu_prob <- array(outer(params$icu_prob, vac_pars$rr_icu), dim=c(N, n_vac, n_strain))
  params$pre_icu <- array(params$pre_icu, dim=c(N, n_vac, n_strain))
  params$post_icu <- array(params$post_icu, dim=c(N, n_vac, n_strain))
 # params$susceptibility <- array(outer(params$susceptibility, vac_pars$rr_inf), dim=c(N, n_vac, n_strain))
  params$susceptibility_symp <- array(outer(params$susceptibility, vac_pars$rr_inf),dim=c(N, n_vac, n_strain))
  params$susceptibility_asymp <- array(outer(params$susceptibility, vac_pars$rr_inf_asymp), dim=c(N, n_vac, n_strain))
  params$transmisibility <- array(outer(params$transmisibility, vac_pars$rr_trans), dim=c(N, n_vac, n_strain))
  params$symp_trans <- array(1, dim=c(N, n_vac, n_strain))

  params$prob_death_non_hosp <- array(outer(params$prob_death_non_hosp, vac_pars$rr_death), dim=c(N, n_vac, n_strain))
  params$prob_death_hosp <- array(outer(params$prob_death_hosp, vac_pars$rr_death), dim=c(N, n_vac, n_strain))
  params$prob_death_icu <- array(outer(params$prob_death_icu, vac_pars$rr_death), dim=c(N, n_vac, n_strain))
  return(params)
}

#' Change the timestep of the parameters
#'
#' @export
change_dt <- function(params, dt){

  params$dt <- dt
  steps_per_day <- 1/dt
  params$N_steps <- params$N_steps/dt
  #Fix Beta
  rows <- 1:nrow(params$beta_day)
  params$beta_day <- params$beta_day[rep(rows, each=steps_per_day),]

  ## Fix import_vec
  imp_vec <- params$import_vec
  imp_vec <- abind(imp_vec, array(0, dim=c(1, dim(imp_vec)[2:4])), along=1)
  rows <- c()
  for(i in 1:(dim(imp_vec)[1]-1)){
    rows <- c(rows, c(i, rep(dim(imp_vec)[1], steps_per_day-1)))
  }
  params$import_vec <- imp_vec[rows,,,]
  dim(params$import_vec) <- c(dim(params$import_vec),1)
  N_vac <- params$vaccinations
  N_vac <- abind(N_vac, array(0, dim=c(1, dim(N_vac)[2:3])), along=1)
  params$vaccinations <- N_vac[rows,,]
  return(params)
}

#'
#' @export
fix_init_params <- function(params){

  
  init_params_full <- c("I_imp_ini", "Ea_ini","Es_ini", "A_ini", "P_ini",
                   "H_ini", "ICU_H_ini", "ICU_R_ini", "ICU_P_ini","B_D_ini", "B_D_H_ini",
                   "B_D_ICU_ini", "R_ini", "D_ini", "tot_infected_ini",  "I_ini",
                   "tot_hosp_ini", "tot_resp_ini")
  init_params_min <- c("S_ini","tot_vac_ini", "tot_vac_adm_ini")

  for(p in c(init_params_full, init_params_min)){
    if(! p %in% names(params)){

      if(p %in% init_params_min){
        dims <-c(params$n, params$n_vac)
      }else{
        dims <- c(params$n, params$n_vac, params$n_strain)
      }
      params[[p]] <- array(0, dim=dims)
    }

  }
  return(params)
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
  if(!"spont_behav_change_params" %in% params_keys) params$spont_behav_change_params <- c(1,100000, 1, 1,1,1)
  if(!"expected_health_loss" %in% params_keys) params$expected_health_loss <- array(0, dim=c(params$n,params$n_vac))
  return(params)
}
  
